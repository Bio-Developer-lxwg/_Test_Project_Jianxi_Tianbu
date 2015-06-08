#include "clsqualitycontrol.h"
#include "clsalgorithm.h"
#include "local_alignment.h"

ClsQualityControl::ClsQualityControl()
{
    m_strReads1Path = "";
    m_strReads2Path = "";
}

ClsQualityControl::~ClsQualityControl()
{
}

void ClsQualityControl::Init(string strReads1Path, string strReads2Path)
{
    m_strReads1Path = strReads1Path;
    m_strReads2Path = strReads2Path;
}

/*Estimate the Gap Length -->inspire by the paper
 *同时这个部分还可以通过gap的长度来判断是否需要作repeat的填充
 * some details
 * 1: Default value is that use the map all --> i think this is make sense.
 *    (1) I think the repeats should allow map all, and let you decide use them
 * 2: The result could be classified into two classes:
 *    (1) Normal: have certain overlap between two special soft clip reads
 *    (2) AbNormal: no overlap
*/
void ClsQualityControl::EstimateGapLength(St_ScaffoldFile& stScaffold,
                                          string strReads1Path,string strReads2Path)
{
    //Step 1:Align the reads to the gap
    string strBamFilePath = ClsAlgorithm::GetInstance().CreateBamFileByMultiScaff(stScaffold.vData,
                                                                                  strReads1Path, strReads2Path,
                                                                                  "RawScafSeq");
    //Step 2: Parse the bam file
    BamReader* pBamReader = new BamReader();
    pBamReader->Open(strBamFilePath);
    pBamReader->OpenIndex(strBamFilePath + ".bai");

    //step 3: record all of position of alignment seq
    BamAlignment al;
    vector<int> clipSizes;
    vector<int> readPositions;
    vector<int> genomePositions;
    //Record SoftClip for Gap Size
    //Init Data Structure
    vector<St_GapAlnSet> vScaffoldSetGapAln;
    vScaffoldSetGapAln.resize(stScaffold.vData.size());
    int iIndex = 0;
    for(vector<St_GapAlnSet>::iterator itr = vScaffoldSetGapAln.begin();
        itr != vScaffoldSetGapAln.end(); itr++)
    {
        itr->vGapReadsSet.resize(stScaffold.vData[iIndex].vGap.size());
        iIndex++;
    }
    St_SoftClipReads stSCReads;
    //Step 4: Set Value
    int iGapIndex = 0;
    while(pBamReader->GetNextAlignment(al))  //这个地方可以尝试去定义相应的quality
    {
        if(al.RefID >= 0) //这个时候证明比上了 --> 相当于比上了某一个scaffold
        {
            unsigned int iLastLen = clipSizes.size();
            al.GetSoftClips(clipSizes, readPositions, genomePositions);
            if(clipSizes.size() - iLastLen > 0) //新增了一个soft clip的 reads
            {
                int iRefPos = *(genomePositions.end()-1);
                //1: 判断Ref Pos是不是在Gap之内
                bool bHitGap = false;
                const int COFFSET = 10;
                iGapIndex = 0;
                for(vector<St_Gap>::iterator itr = stScaffold.vData[al.RefID].vGap.begin();
                    itr != stScaffold.vData[al.RefID].vGap.end(); itr++)
                {
                    if((iRefPos-COFFSET >= itr->iStartPos && iRefPos-COFFSET <= itr->iEndPos) ||
                       (iRefPos+COFFSET >= itr->iStartPos && iRefPos+COFFSET <= itr->iEndPos))
                    {
                        bHitGap = true;
                        break;
                    }
                    iGapIndex++;
                }
                //2: 如果hit上了gap，那么我们是左边的softclip还是右边的softclip，然后对其进行记录
                if(bHitGap) //hit it
                {
                    stSCReads.iRefPos = iRefPos;
                    //判断是哪一边的
                    int icurClipLen = *(clipSizes.end()-1);
                    int icurClipPos = *(readPositions.end()-1);
                    if(icurClipPos > al.Length)
                        continue;
                    if(icurClipLen == icurClipPos)
                        icurClipPos = 0;
                    stSCReads.iClipStart = icurClipPos;
                    if(icurClipPos == 0) //this tell us : the softclip is start from the beginning
                    {
                        stSCReads.enClipPart = cpLeft;
                        stSCReads.iMissingEnd = stSCReads.iRefPos - 1;
                        stSCReads.iClipLen = icurClipLen;
                        stSCReads.strClipSeq = al.QueryBases.substr(0, stSCReads.iClipLen);
                        vScaffoldSetGapAln[al.RefID].vGapReadsSet[iGapIndex].vLeftClip.push_back(stSCReads);
                        //vScaffoldSetGapAln[al.RefID].bAddNoneMap = false;
                    }
                    else
                    {
                        stSCReads.enClipPart = cpRight;
                        stSCReads.iMissingStart = stSCReads.iRefPos + 1;
                        stSCReads.iClipLen = icurClipLen;
                        stSCReads.strClipSeq = al.QueryBases.substr(stSCReads.iClipStart, stSCReads.iClipLen);
                        //vScaffoldSetGapAln[al.RefID].bAddNoneMap = true;
                        vScaffoldSetGapAln[al.RefID].vGapReadsSet[iGapIndex].vRightClip.push_back(stSCReads);
                    }
                }
            }
        }
        /*
        else if(vScaffoldSetGapAln[al.RefID].bAddNoneMap)
        {
            if(al.IsMateMapped())
                vScaffoldSetGapAln[al.RefID].vGapReadsSet[iGapIndex].vNoneMap.push_back(St_NoneMap(al.MatePosition));
        }*/
    }
    //Step 5: After Get All those Valid Reads --> let's estimate the length of the gap
    for(vector<St_GapAlnSet>::iterator itr = vScaffoldSetGapAln.begin();
        itr != vScaffoldSetGapAln.end(); itr++)
    {
        for(vector<St_GapReads>::iterator itrGap = itr->vGapReadsSet.begin();
            itrGap != itr->vGapReadsSet.end(); itrGap++)
        {
            //Step a): Find the most right soft clip reads in vRightClip
            int iRightClipMaxLen = 0;
            string strRightClip = "";
            for(vector<St_SoftClipReads>::iterator itrRightClip = itrGap->vRightClip.begin();
                itrRightClip != itrGap->vRightClip.end(); itrRightClip++)
            {
                if(iRightClipMaxLen < itrRightClip->iClipLen)
                {
                    iRightClipMaxLen = itrRightClip->iClipLen;
                    strRightClip = itrRightClip->strClipSeq;
                }
            }
        #ifdef __Testing
            cout << "Right Clip Max: " << strRightClip << endl;
        #endif
            //Step b): Find the most left Soft clip reads in vLeftClip
            int iLeftClipMaxLen = 0;
            string strLeftClip = "";
            for(vector<St_SoftClipReads>::iterator itrLeftClip = itrGap->vLeftClip.begin();
                itrLeftClip != itrGap->vLeftClip.end(); itrLeftClip++)
            {
                if(iLeftClipMaxLen < itrLeftClip->iClipLen)
                {
                    iLeftClipMaxLen = itrLeftClip->iClipLen;
                    strLeftClip = itrLeftClip->strClipSeq;
                }
            }
        #ifdef __Testing
            cout << "Left Clip Max: " << strLeftClip << endl;
        #endif
            //Step c): try to find the number of mate with the different postion in vNoneMap
            //int iNoneClipMaxLen = 0;
            //int iTimes = 0;
            /*
            for(vector<St_NoneMap>::iterator itrNoneMap = itrGap->vNoneMap.begin();
                itrNoneMap != itrGap->vNoneMap.end(); itrNoneMap++)
            {
                iTimes++;
            }*/
            //if(iTimes > 0)
            //    iNoneClipMaxLen = 100 + iTimes;
            //Record Length
            //1： Do local alignment -->
            int iLen = 0;
            bool bAbnormalLarge = false;
            if(strLeftClip != "" && strRightClip != "")
            {
                LocalAlignment localAl;
                int iStartOrg = -1, iEndOrg = -1, iStartCmp = -1, iEndCmp = -1;
                int iScore = -1;
                string strAlnSymbol = "";
                localAl.optAlignEx(strRightClip, strLeftClip,
                                 iStartOrg, iEndOrg, iStartCmp, iEndCmp,
                                 iScore, strAlnSymbol);
                //Check if such align is reliable
                int iMinLenStr = strRightClip.length() > strLeftClip.length() ? strLeftClip.length() :
                                                                                strRightClip.length();
                if(iScore >= (iMinLenStr * .5))
                {
                    //这个case我们认为是存在了overlap，那么我们可以如下进行处理
                    iLen = CalGapLenghByOverlap(strRightClip, strLeftClip, iStartOrg, iEndOrg,
                                                    iStartCmp, iEndCmp);
                }
                else
                {
                    iLen = strRightClip.length() + strLeftClip.length();
                    bAbnormalLarge = true;
                }
            }
            else
            {
                if(strLeftClip == "")
                    iLen = strRightClip.length();
                else if(strRightClip == "")
                    iLen = strLeftClip.length();
                bAbnormalLarge = true;
            }
            itrGap->iLength = iLen;
            itrGap->bAbnormalLarge = bAbnormalLarge;
            //2： Check score to decide which kind case of length should be calculate
            cout << IntToStr(itrGap->iLength) << endl;
            cout << "Right Clip: " << strRightClip << endl;
            cout << "Left Clip: " << strLeftClip << endl;
        }
    }
    //在这里可以update的相应的gap的长度
    int iScfIndex = 0;
    for(vector<St_ScaffoldUnit>::iterator itr = stScaffold.vData.begin();
        itr != stScaffold.vData.end(); itr++)
    {
        iGapIndex = 0;
        for(vector<St_Gap>::iterator subItr = itr->vGap.begin();
            subItr != itr->vGap.end(); subItr++)
        {
            subItr->iEstimateLen = vScaffoldSetGapAln[iScfIndex].vGapReadsSet[iGapIndex].iLength;
            subItr->bAbnormalLarge = vScaffoldSetGapAln[iScfIndex].vGapReadsSet[iGapIndex].bAbnormalLarge;
        }
        iScfIndex++;
    }
    //Release memory
    delete pBamReader;
    pBamReader = NULL;
}

//注意：left的clip意味着右边的seq，right的clip意味着左边的序列
int ClsQualityControl::CalGapLenghByOverlap(string& strLeftSeq, string& strRightSeq,
                                             int iLStart, int iLEnd, int iRStart, int iREnd)
{
    //有好多个cases
    //for left seq
    En_AlnState enLeftAlnState = apMax;
    if(iLEnd >= (int)strLeftSeq.length()-6 && iLStart <= 5)
        enLeftAlnState = apWhole;
    else if(iLEnd >= (int)strLeftSeq.length()-6) // 在留下5bp的缓冲区
        enLeftAlnState = apRight;
    else if(iLStart <= 5)
        enLeftAlnState = apLeft;
    else
        enLeftAlnState = apMiddle;
    //For Right Seq
    En_AlnState enRightAlnState = apMax;
    if(iREnd >= (int)strLeftSeq.length()-6 && iRStart <= 5)
        enRightAlnState = apWhole;
    else if(iREnd >= (int)strLeftSeq.length()-6) // 在留下5bp的缓冲区
        enRightAlnState = apRight;
    else if(iRStart <= 5)
        enRightAlnState = apRight;
    else
        enRightAlnState = apMiddle;
    //根据不同的状态进行合并考虑
    int iLen = 0;
    switch((int)enLeftAlnState)
    {
        case apWhole:
            switch((int)enRightAlnState)
            {
                case apWhole:
                    iLen = strLeftSeq.length() + (strRightSeq.length() - iREnd - 1);
                    break;
                case apLeft:
                    iLen = strLeftSeq.length() + (strRightSeq.length() - iREnd - 1);
                    break;
                case apRight:
                    iLen = strRightSeq.length() - iRStart;
                    break;
                case apMiddle:
                    iLen = strRightSeq.length() - iRStart;
                    break;
            }
            break;
        case apLeft:
            switch((int)enRightAlnState)
            {
                case apWhole:
                    iLen = strRightSeq.length();
                    break;
                case apLeft:
                    iLen = strRightSeq.length();
                    break;
                case apRight:
                    iLen = iLEnd + (strRightSeq.length() - iREnd);
                    break;
                case apMiddle:
                    iLen = strRightSeq.length() - iRStart;
                    break;
            }
            break;
        case apRight:
            switch((int)enRightAlnState)
            {
                case apWhole:
                    iLen = iLEnd;
                    break;
                case apLeft:
                    iLen = iLEnd + (strRightSeq.length() - iREnd);
                    break;
                case apRight:
                    iLen = iLEnd;
                    break;
                case apMiddle:
                    iLen = iLEnd + (strRightSeq.length() - iREnd);
                    break;
            }
            break;
        case apMiddle:
            switch((int)enRightAlnState)
            {
                case apWhole:
                    iLen = iLEnd;
                    break;
                case apLeft:
                    iLen = iLEnd + strRightSeq.length() - iREnd;
                    break;
                case apRight:
                    iLen = iREnd - iRStart + iLStart;
                    break;
                case apMiddle:
                    iLen = iLStart + iREnd - iRStart + strRightSeq.length() - iREnd;
                    break;
            }
            break;
    }
    return iLen;
}

bool ClsQualityControl::EstimateFillGapByRepeat(St_RepeatUnit& stRepatUnit, St_Gap& stGap)
{
    //Left Flank------>
    int iKmerMatch = 0;
    bool bLeftKmerMatch = false;
    for(vector<KmerTypeShort>::iterator itr = stGap.vLeftKmer.begin(); itr != stGap.vLeftKmer.end(); itr++)
    {
        //Normal
        bool bGet = false;
        for(vector<KmerTypeShort>::iterator itrRep = stRepatUnit.vNormalKmerValue.begin();
            itrRep != stRepatUnit.vNormalKmerValue.end(); itrRep++)
        {
            if(*itr == *itrRep)
            {
                bGet = true;
                stGap.enLFMatchType = mtNormal;
                break;
            }
        }
        if(!bGet) //Complementary Reverse
        {
            for(vector<KmerTypeShort>::iterator itrRep = stRepatUnit.vCompleRevKmerValue.begin();
                itrRep != stRepatUnit.vCompleRevKmerValue.end(); itrRep++)
            {
                if(*itr == *itrRep)
                {
                    bGet = true;
                    stGap.enLFMatchType = mtCompleRev;
                    break;
                }
            }
        }
        if(bGet)
            iKmerMatch++;
        if(iKmerMatch > 3)
        {
            bLeftKmerMatch = true;
            break;
        }
    }

    //Right Flank-------->
    iKmerMatch = 0;
    bool bRightKmerMatch = false;
    for(vector<KmerTypeShort>::iterator itr = stGap.vRightKmer.begin(); itr != stGap.vRightKmer.end(); itr++)
    {
        //Normal
        bool bGet = false;
        for(vector<KmerTypeShort>::iterator itrRep = stRepatUnit.vNormalKmerValue.begin();
            itrRep != stRepatUnit.vNormalKmerValue.end(); itrRep++)
        {
            if(*itr == *itrRep)
            {
                bGet = true;
                stGap.enRFMatchType = mtNormal;
                break;
            }
        }
        if(!bGet) //Complementary Reverse
        {
            for(vector<KmerTypeShort>::iterator itrRep = stRepatUnit.vCompleRevKmerValue.begin();
                itrRep != stRepatUnit.vCompleRevKmerValue.end(); itrRep++)
            {
                if(*itr == *itrRep)
                {
                    bGet = true;
                    stGap.enRFMatchType = mtCompleRev;
                    break;
                }
            }
        }
        if(bGet)
            iKmerMatch++;
        if(iKmerMatch > 3)
        {
            bRightKmerMatch = true;
            break;
        }
    }

    if(bLeftKmerMatch && bRightKmerMatch)
        stGap.enGapType = gtCenter;
    else if(bLeftKmerMatch)
        stGap.enGapType = gtRight;
    else if(bRightKmerMatch)
        stGap.enGapType = gtLeft;
    else{}

    return bLeftKmerMatch || bRightKmerMatch;
}

//It has been abandoned
void ClsQualityControl::CheckFillQuality(St_ScaffoldUnit& stScaffoldUnit)
{
    //No gap need to be filled
    if(stScaffoldUnit.vGapRefinedByRept.empty())
    {
        cout << "There is no gap need to be filled!";
        return;
    }
    //=========================>
    //Step 1: Generate the final refined scaffolds -->Where all of gaps has been fixed
    //Step 2: Map the reads back to refined scaffold
    //Step 3: Check the gap area and identify its coverage
    //step 4: Get the whole coverage value and find out if the coverage of gap is reasonable
    //step 5: Update the final result
    //<=========================Go!!!!!

    St_FinalScaffoldUnit stFinalScaffoldUnit; //Final Scaffold Unit
    St_FillUnit stFillUnit;

    stFinalScaffoldUnit.strName = stScaffoldUnit.strName;
    //Step 1: fill the data for stFinalScaffoldUnit and stFillUnit
    //1: build the final scaffold by the sequence before gap
    stFinalScaffoldUnit.strFillByNorm += stScaffoldUnit.strSeq.substr(0,
                                                   stScaffoldUnit.vGapRefinedByRept[0].pOrgGap->iStartPos);
    stFinalScaffoldUnit.strFillByRevComp += stScaffoldUnit.strSeq.substr(0,
                                                   stScaffoldUnit.vGapRefinedByRept[0].pOrgGap->iStartPos);

    for(vector<St_GapRefinedByRept> ::iterator itr = stScaffoldUnit.vGapRefinedByRept.begin();
        itr != stScaffoldUnit.vGapRefinedByRept.end(); itr++) // traverse all of gaps(those gaps are refined gaps)
    {
        //Add the repeat
        stFillUnit.iStart = stFinalScaffoldUnit.strFillByNorm.length() - 1; //norm仅仅是意味着是正常的方向，并不是意味着修补的方式
        stFinalScaffoldUnit.strFillByNorm += itr->pOrgGap->strFillSeq; //itrSCR->pGap->pOrgGap->strFillSeq; //Fill by Repeat
        stFinalScaffoldUnit.strFillByNorm += itr->stBamGap.strFillSeq; //itrSCR->stBamGap.strFillSeq; //Fill by Paired-End reads
        //if(itrSCR->pGap->pOrgGap->strFillSeq == "" && itrSCR->strFillSeq == "")
        if(itr->pOrgGap->strFillSeq == "" && itr->stBamGap.strFillSeq == "")
        {
            //if both filling method are failed to fill any sequences--->Add the gap back to it.
            //for(int i=0; i<itrSCR->pGap->pOrgGap->iLen; i++)
            for(int i=0; i<itr->pOrgGap->iLen; i++)
            {
                stFinalScaffoldUnit.strFillByNorm += "N";
            }
        }
        //因为暂时还没有发现repeats是reverse complementary的情况，因此我们在此虽然仅仅是
        stFinalScaffoldUnit.strFillByRevComp += itr->pOrgGap->strFillSeq;//Fill by Repeat
        stFinalScaffoldUnit.strFillByRevComp += itr->stBamGap.strRevCompFillSeq;
        if(itr->pOrgGap->strFillSeq == "" && itr->stBamGap.strRevCompFillSeq == "")
        {
            //if both filling method are failed to fill any sequences--->Add the gap back to it.
            for(int i=0; i<itr->pOrgGap->iLen; i++)
            {
                stFinalScaffoldUnit.strFillByRevComp += "N";
            }
        }
        stFillUnit.iEnd = stFinalScaffoldUnit.strFillByNorm.length() - 1;
        stFinalScaffoldUnit.vNormGapRecord.push_back(stFillUnit);//only to cases will affect the result of gap filling: repeats and pair-end reads
        stFinalScaffoldUnit.vRevCompGapRecord.push_back(stFillUnit);//only to cases will affect the result of gap filling: repeats and pair-end reads
        //Add the sequence between current gap and the next gap
        int iLen = -1;
        int iStartPos = -1;

        //将中间的部分填补起来
        if(itr + 1 < stScaffoldUnit.vGapRefinedByRept.end())
        {
            iLen = (itr+1)->pOrgGap->iStartPos - itr->pOrgGap->iEndPos - 1;
            iStartPos = itr->pOrgGap->iEndPos+1;
        }
        else // the last one then
        {
            iLen = stScaffoldUnit.strSeq.length() - itr->pOrgGap->iEndPos - 1;
            iStartPos = itr->pOrgGap->iEndPos+1;
        }
        stFinalScaffoldUnit.strFillByNorm += stScaffoldUnit.strSeq.substr(iStartPos, iLen);
        stFinalScaffoldUnit.strFillByRevComp += stScaffoldUnit.strSeq.substr(iStartPos, iLen);
        break;
    }
    //Step2: Get Map Read Back to Refined Scaffold to check the coverage in gao area
    //1: Clear the folder rm -rf path/*  do not needed, since this function has been added into bam file creation
    //2: For the Norm Sequence
    //  (1) Make the bam file
    //  (2) Analysis the result
    string strBamFilePath = ClsAlgorithm::GetInstance().CreateBamFile(stFinalScaffoldUnit.strName,
                                                                      stFinalScaffoldUnit.strFillByNorm,
                                                                      m_strReads1Path, m_strReads2Path);
    QualityEstimate(strBamFilePath, stFinalScaffoldUnit.vNormGapRecord,
                    stFinalScaffoldUnit.strFillByNorm);
    //3: For the reverse complementary Sequence
    //  (1) Make the bam file
    //  (2) Analysis the result
    strBamFilePath = ClsAlgorithm::GetInstance().CreateBamFile(stFinalScaffoldUnit.strName,
                                                               stFinalScaffoldUnit.strFillByRevComp,
                                                               m_strReads1Path, m_strReads2Path);
    QualityEstimate(strBamFilePath, stFinalScaffoldUnit.vRevCompGapRecord,
                    stFinalScaffoldUnit.strFillByRevComp);

    //Step3: Collect the result and confirm if the gap shold be filled by normal sequence or the reverse complementary
    //Empty --> Return directly
    if(stFinalScaffoldUnit.vNormGapRecord.empty())
    {
        stFinalScaffoldUnit.strFill = stFinalScaffoldUnit.strFillByNorm;
        return;
    }
    //Not empty-->Refine it
    stFinalScaffoldUnit.strFill += stFinalScaffoldUnit.strFillByNorm.substr(0,
                                           stFinalScaffoldUnit.vNormGapRecord[0].iStart);

    vector<St_FillUnit>::iterator itrRevComp = stFinalScaffoldUnit.vRevCompGapRecord.begin();
    for(vector<St_FillUnit>::iterator itrNorm = stFinalScaffoldUnit.vNormGapRecord.begin();
        itrNorm < stFinalScaffoldUnit.vNormGapRecord.end(); itrNorm++, itrRevComp++)
    {
        //We mainly depends on the sequence of Norm and then try to refine it based on RecComp
        //Add the gap into final sequence
        stFillUnit.iStart = stFinalScaffoldUnit.strFill.length()-1;
        if(itrNorm->enQuality > itrRevComp->enQuality)
        {
            stFinalScaffoldUnit.strFill += stFinalScaffoldUnit.strFillByNorm.substr(itrNorm->iStart,
                                                                 itrNorm->iEnd - itrNorm->iStart + 1);
            stFillUnit.enQuality = itrNorm->enQuality;
            stFillUnit.iCoverage = itrNorm->iCoverage;
        }
        else
        {
            stFinalScaffoldUnit.strFill += stFinalScaffoldUnit.strFillByRevComp.substr(itrRevComp->iStart,
                                                                 itrRevComp->iEnd - itrRevComp->iStart + 1);
            stFillUnit.enQuality = itrRevComp->enQuality;
            stFillUnit.iCoverage = itrRevComp->iCoverage;
        }
        stFillUnit.iEnd = stFinalScaffoldUnit.strFill.length()-1;
        stFinalScaffoldUnit.vFinalFillGapRecord.push_back(stFillUnit);

        //Add the normal seq into final sequence
        int iStart = -1;
        int iLen = -1;
        if(itrNorm + 1 < stFinalScaffoldUnit.vNormGapRecord.end()) //Not the last item
        {
            iStart = itrNorm->iEnd + 1;
            iLen = (itrNorm + 1)->iStart - itrNorm->iEnd - 1;
        }
        else //The last item
        {
            iStart = itrNorm->iEnd + 1;
            iLen = stFinalScaffoldUnit.strFillByNorm.length() - itrNorm->iEnd - 1;
        }
        stFinalScaffoldUnit.strFill += stFinalScaffoldUnit.strFillByNorm.substr(iStart, iLen);
    }
    //return stFinalScaffoldUnit;
    //m_vFinalResult.push_back(stFinalScaffoldUnit);
}

void ClsQualityControl::QualityEstimate(string strBamFilePath,
                                        vector<St_FillUnit>& vFillUnit, string& strRef)
{
    BamReader* pBamReader = new BamReader();
    pBamReader->Open(strBamFilePath);
    pBamReader->OpenIndex(strBamFilePath + ".bai");

    BamAlignment al;
    vector<int> clipSizes;
    vector<int> readPositions;
    vector<int> genomePositions;
    /*Target:
     *1:Collect the reads which contains the gap
     *2:Calculate the coverage
     *3:Estimate the total coverage
     *4:Set the quality based on current coverage and the total average coverage
     */
    while(pBamReader->GetNextAlignment(al))
    {
        if(al.AlignedBases == "")
            continue;
        int iStart1 = al.Position;
        int iEnd1  = al.GetEndPosition();
        if(iStart1 > iEnd1)
            continue;
        for(vector<St_FillUnit>::iterator itr = vFillUnit.begin(); itr != vFillUnit.end(); itr++)
        {
            int iLen = itr->iEnd - itr->iStart + 1;
            //int iStart2 = itr->iStart; // this is the repeats
            //int iEnd2 = itr->iEnd;

            //int iStart2 = itr->iStart - iLen / 2;
            //int iEnd2 = itr->iStart + iLen / 2;

            int iStart2 = itr->iEnd - iLen / 2;
            int iEnd2 = itr->iEnd + iLen / 2;

            //cout << strRef.substr(iStart2, iEnd2 - iStart2 + 1) << endl;
            itr->iCoverage += ClsAlgorithm::GetInstance().CheckInterSection(iStart1, iEnd1,
                                                                            iStart2, iEnd2);
        }
    }
    //Set Quality
    //this is a problem!!!! --> we need to use the break point to adjust the result!!!!
    //rather than just use the pair-end read to map the current sequence
    for(vector<St_FillUnit>::iterator itr = vFillUnit.begin(); itr != vFillUnit.end(); itr++)
    {
        int iLen = itr->iEnd - itr->iStart + 1;
        //int iStart2 = itr->iStart - iLen / 2;
        //int iEnd2 = itr->iStart + iLen / 2;
        int iStart2 = itr->iEnd - iLen / 2;
        int iEnd2 = itr->iEnd + iLen / 2;
        cout << strRef.substr(iStart2, iEnd2 - iStart2 + 1) << endl;

        itr->iCoverage = itr->iCoverage / (itr->iEnd - itr->iStart + 1); //Since there are two strands
        //cout << strRef.substr(itr->iStart, itr->iEnd - itr->iStart + 1) << endl;
        cout << itr->iCoverage << endl;
        if(itr->iCoverage > 30)
            itr->enQuality = fqHigh;
        else if(itr->iCoverage <=30 && itr->iCoverage > 10)
            itr->enQuality = fqNormal;
        else if(itr->iCoverage <=10 && itr->iCoverage >= 5)
            itr->enQuality = fqLow;
    }
    cout << "===========End==========" << endl;
    delete pBamReader;
    pBamReader = NULL;
}

