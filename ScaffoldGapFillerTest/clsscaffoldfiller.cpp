#include "clsscaffoldfiller.h"
#include "local_alignment.h"
#include "smith-waterman.h"
#include "stdlib.h"
#include <unistd.h>
#include "clsalgorithm.h"
#include "algorithm"

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#define __Testing

ClsScaffoldFiller::ClsScaffoldFiller(): m_pFastaParse(NULL)
{
}

ClsScaffoldFiller::~ClsScaffoldFiller()
{
    if(m_pFastaParse)
    {
        delete m_pFastaParse;
        m_pFastaParse = NULL;
    }
}

void ClsScaffoldFiller::Init(char **argv)
{
    if(m_pFastaParse == NULL)
    {
        m_pFastaParse = new FastaParser();
    }
    //Parameter 1: Repeat File
    //Parameter 2: scaffoldFile
    //Parameter 3: reads1
    //Parameter 4: reads2
    m_stRepeat.Init(argv[1], m_pFastaParse); //Init Repeat File
    m_stScaffold.Init(argv[2], m_pFastaParse, m_stRepeat.GetMaxLen()); //Init Scaffold File
    m_strReads1Path = argv[3]; //Set Reads1
    m_strReads2Path = argv[4]; //Set Reads2
    cout << "m_stScaffold Init finished" << endl;
}

void ClsScaffoldFiller::FillScaffold()
{
    if(m_stScaffold.vData.empty())
        return;
    m_vFinalResult.clear();
    //move the file -->Which record the result of gap filling
    string strResultPath = ClsAlgorithm::GetInstance().GetCurExeFolderPath() + "/GapFillerInfo";
    ::remove(strResultPath.c_str());
    //move all of files in temperary folder
    string strCmd = "exec rm -r " +
                    ClsAlgorithm::GetInstance().GetHigherFolderPath(ClsAlgorithm::GetInstance().GetCurExeFolderPath()) +
                    "TempFile/*";
    system(strCmd.c_str());    
    //第一步，首先进行GapLength的评估---->这个是受到了那个paper的启示
    EstimateGapLength();
    //将结果输出来
#ifdef __Testing
    int iScaffoldIndex = 1;
    for( vector<St_ScaffoldUnit>::iterator itr = m_stScaffold.vData.begin();
         itr != m_stScaffold.vData.end(); itr++)
    {
        int iGapIndex = 0;
        for(vector<St_Gap>::iterator itrGap = itr->vGap.begin(); itrGap != itr->vGap.end(); itrGap++)
        {
            if(itrGap->bAbnormalLarge)
            {
                cout << "Scaffold[" << IntToStr(iScaffoldIndex) << "]: " << endl;
                cout <<"______Gap[" << IntToStr(iGapIndex) << "] Abmormal: " << "Yes" << endl;
            }
            iGapIndex++;
        }
        iScaffoldIndex++;
    }
#endif

    //在此处我想在服务器上做一个实验，也就是

    //然后我们需要做的是对剩下的相应的scaffold进行repeat的填补
    //Step 2: Fill the gap for each scaffold unit
    //整体逻辑应该是
    //(1)一次性先把scaffold 用repeat填补完
    for(vector<St_ScaffoldUnit>::iterator itr = m_stScaffold.vData.begin();
        itr != m_stScaffold.vData.end(); itr++)
    {
        //Step 1: try to use repeat to fill the gap
        FillScaffoldUnitByRepeats(*itr);
        //Step 2: 用其他scaffold里面的去填补完成
        FillScaffoldUnitByDraftGeno(*itr);
    }
    //(3)然后用pair-end reads第三次的填补: 分开的原因是因为，这一步需要作bwa，
    //   而input的数据是需要通过上面两个算法进行fixed后的数据
    FillScaffoldByPairEndReads(m_stScaffold.vData);
    //(4)最后再是用pair-end reads去进行相应的结果验证

    //(5) 将填充结果展现出来
    ofstream ofs;
    string strRootPath = ClsAlgorithm::GetInstance().GetCurExeFolderPath();
    ofs.open((strRootPath + "/GapFillerInfo").c_str(), _S_app); //Do not need to use the "a+", since we write such file just one time

    for(vector<St_ScaffoldUnit>::iterator itr = m_stScaffold.vData.begin();
        itr != m_stScaffold.vData.end(); itr++)
    {
        ofs << ">" << itr->strName << endl;
        int i = 0;
        for(vector<St_GapRefinedByRept>::iterator subItr = itr->vGapRefinedByRept.begin();
            subItr != itr->vGapRefinedByRept.end(); subItr++)
        {
            ofs << "Gap Index [" << IntToStr(i) << "]:" << endl;
            ofs << "Fill by Repeats ==>" << subItr->pOrgGap->strFillSeq << endl;
            ofs << "Final_PE: Fill by PE Reads: Norm Ore ==> " << subItr->stBamGap.strFillSeq << endl;
            ofs << "Final_PE: Fill by PE Reads: Reverse Complementary ==> " << subItr->stBamGap.strRevCompFillSeq << endl;
            ofs << "Temp_PE: left Extend: " << subItr->stBamGap.strLeftExt << endl;
            ofs << "Temp_PE: right Extend: " << subItr->stBamGap.strRightExt << endl;
            ofs << "=========================================" << endl;
            i++;
        }
    }
    ofs.close();
}

//Estimate the Gap Length -->inspire by the paper
//同时这个部分还可以通过gap的长度来判断是否需要作repeat的填充
void ClsScaffoldFiller::EstimateGapLength()
{
    //Step 1:Align the reads to the gap
    string strBamFilePath = ClsAlgorithm::GetInstance().CreateBamFileByMultiScaff(m_stScaffold.vData,
                                                                                  m_strReads1Path, m_strReads2Path,
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
    vScaffoldSetGapAln.resize(m_stScaffold.vData.size());
    int iIndex = 0;
    for(vector<St_GapAlnSet>::iterator itr = vScaffoldSetGapAln.begin();
        itr != vScaffoldSetGapAln.end(); itr++)
    {
        itr->vGapReadsSet.resize(m_stScaffold.vData[iIndex].vGap.size());
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
                for(vector<St_Gap>::iterator itr = m_stScaffold.vData[al.RefID].vGap.begin();
                    itr != m_stScaffold.vData[al.RefID].vGap.end(); itr++)
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
            int iTimes = 0;
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
    for(vector<St_ScaffoldUnit>::iterator itr = m_stScaffold.vData.begin();
        itr != m_stScaffold.vData.end(); itr++)
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
int ClsScaffoldFiller::CalGapLenghByOverlap(string& strLeftSeq, string& strRightSeq,
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

/*
void ClsScaffoldFiller::FillScaffoldUnit(St_ScaffoldUnit& stScaffoldUnit)
{
    //Step 1: try to use repeat to fill the gap
    //FillScaffoldUnitByRepeats(stScaffoldUnit);
    //Step 2: Use Pair end Reads to fix the gap:
    //FillScaffoldUnitByPairEndReads(stScaffoldUnit);

    //--------->Just for temperay output: we try to print out those result
    //In this phase: we just want to list all the gaps which could be filled by repeats and pair-end reads
    ofstream ofs;
    string strRootPath = ClsAlgorithm::GetInstance().GetCurExeFolderPath();
    ofs.open((strRootPath + "/GapFillerInfo").c_str(), _S_app); //Do not need to use the "a+", since we write such file just one time
    {
        ofs << ">" << stScaffoldUnit.strName << endl;
        int i = 0;
        for(vector<St_GapRefinedByRept>::iterator itr = stScaffoldUnit.vGapRefinedByRept.begin();
            itr != stScaffoldUnit.vGapRefinedByRept.end(); itr++)
        {
            ofs << "Gap Index [" << IntToStr(i) << "]:" << endl;
            ofs << "Fill by Repeats ==>" << itr->pOrgGap->strFillSeq << endl;
            ofs << "Final_PE: Fill by PE Reads: Norm Ore ==> " << itr->stBamGap.strFillSeq << endl;
            ofs << "Final_PE: Fill by PE Reads: Reverse Complementary ==> " << itr->stBamGap.strRevCompFillSeq << endl;
            ofs << "Temp_PE: left Extend: " << itr->stBamGap.strLeftExt << endl;
            ofs << "Temp_PE: right Extend: " << itr->stBamGap.strRightExt << endl;
            ofs << "=========================================" << endl;
            i++;
        }
    }
    //<----------------

    //Step 3: Check the Quality of fill result    
    //CheckFillQuality(stScaffoldUnit);
}*/

void ClsScaffoldFiller::FillScaffoldUnitByRepeats(St_ScaffoldUnit& stScaffoldUnit)
{
    //Go!!!!-->
    for(vector<St_Gap>::iterator itrGap = stScaffoldUnit.vGap.begin();
        itrGap != stScaffoldUnit.vGap.end(); itrGap++)
    {
        if(!itrGap->bAbnormalLarge) //do not fill by repeats if the length is not abnormal large
            continue;

        if(!itrGap->CouldFill())
            continue;

        //Possible fill then ---->
        //Use Normal Repeat and complement reverse repeat to fill the Gap (both those two type has been contained by St_RepeatUnit)
        for(vector<St_RepeatUnit>::iterator itr = m_stRepeat.vData.begin();
            itr != m_stRepeat.vData.end(); itr++)
        {
            //For each repeat --> We need to check if it could be used to fix the gap
            if(itrGap->bFilled)
                continue;
            if(EstimateFillGapByRepeat(*itr, *itrGap))
            {
                FillGapByRepeat(*itr, *itrGap);
                if(itrGap->bFilled)
                {
                    cout<<"=================Gap Filler===============" << endl;
                    cout<<"Scaffold Name: " << stScaffoldUnit.strName << endl;
                    cout<<"Filled String: " << itrGap->strFillSeq << endl;
                    cout<<"From the \'Start\' Pos of Repeat: " << IntToStr(itrGap->iReptStart) << endl;
                    cout<<"From the \'End\' Pos of Repeat: " << IntToStr(itrGap->iReptEnd) << endl;
                    cout<<"Filled \'Length\': " << IntToStr(itrGap->strFillSeq.length()) << endl;
                }
            }
        }
    }
    //We need to update the refined sequence
    stScaffoldUnit.UpdateRefinedSeqByRepeats();
}

//使用draft geno去填补当前的gap
void ClsScaffoldFiller::FillScaffoldUnitByDraftGeno(St_ScaffoldUnit& stScaffoldUnit)
{
    return;
}

bool ClsScaffoldFiller::EstimateFillGapByRepeat(St_RepeatUnit& stRepatUnit, St_Gap& stGap)
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

void ClsScaffoldFiller::FillGapByRepeat(St_RepeatUnit& stRepatUnit, St_Gap& stGap)
{
    //Get the type of match:
    En_MatchType enMatchType = mtNone;
    if(stGap.enLFMatchType == mtNone)
    {
        if(stGap.enRFMatchType != mtNone)
            enMatchType = stGap.enRFMatchType;
    }
    else if(stGap.enRFMatchType == mtNone)
    {
        if(stGap.enLFMatchType != mtNone)
            enMatchType = stGap.enLFMatchType;
    }
    else // both are not none
    {
        if(stGap.enLFMatchType == stGap.enRFMatchType)
            enMatchType = stGap.enLFMatchType;
    }
    if(enMatchType == mtNone) // 不允许左右作match后出现不一致结果的case继续进行
        return;
    string& strOrg = (enMatchType == mtNormal) ? stRepatUnit.strNormalSeq : stRepatUnit.strCompleRevSeq;
    //实际上Kmer感觉是用来判定是否可行的，具体怎么弄还是要用动态规划来做
    // evalaute whether the area can be filled by the passed-in repeat
    // use dynamic programming; use a simple score as follows
    // match: 1, mismatch/indel/: -1
    // middle gap: no penalty if length is within some threshold (+/- say 20%)
    // otherwise: -inf
    // tbl[x,y]: x is coordinate in repeats, y: scalfold (relative to the start)
    // define scoring scheme (simple for now)
    int iLeftFillPos = -1;
    int iRightFillPos = -1;
    switch((int)stGap.enGapType)
    {
        case gtLeft:
            iLeftFillPos = 0;
            iRightFillPos = GetPosFromLocalAlignment(strOrg, stGap.strRefinedRightFlank, gtLeft);;
            break;
        case gtRight:
            iLeftFillPos = GetPosFromLocalAlignment(strOrg, stGap.strRefinedLeftFlank, gtRight);
            iRightFillPos = strOrg.length()-1;
            break;
        case gtCenter:
            iLeftFillPos = GetPosFromLocalAlignment(strOrg, stGap.strRefinedLeftFlank, gtCenter, 0);
            iRightFillPos = GetPosFromLocalAlignment(strOrg, stGap.strRefinedRightFlank, gtCenter, 1);
            break;
    }
    if(iLeftFillPos < 0 || iRightFillPos < 0 || iLeftFillPos > iRightFillPos)
        return;
    stGap.iReptStart = iLeftFillPos;
    stGap.iReptEnd = iRightFillPos;
    stGap.strFillSeq = strOrg.substr(iLeftFillPos, iRightFillPos - iLeftFillPos + 1);
    stGap.bFilled = true;
}

int ClsScaffoldFiller::GetPosFromLocalAlignment(string& OrgStr, string& strCmp, En_GapType enType, int iTag)
{
    LocalAlignment al;//--> 使用现有的库来进行Local Alignment
    int iStartOrg = -1, iEndOrg = -1, iStartCmp = -1, iEndCmp = -1;
    al.optAlign(OrgStr, strCmp, iStartOrg, iEndOrg, iStartCmp, iEndCmp); //如果都对上应该是那么分数就应该是整个的长度
    int iValue = -1;
    switch((int)enType)
    {
        case gtLeft:
            iValue = iStartOrg - 2; // start from 0 and should be use the character before such value
            break;
        case gtRight:
            iValue = iEndOrg;
            break;
        case gtCenter:
        {
            if(0 == iTag) // 左边的起点
                iValue = iEndOrg; //start from 1 and use the character after such value
            else // if iTag = 1
                iValue = iStartOrg -2; // start from 0 and should be use the character before such value

        }
    }
    return iValue;
}

//Step 2: Notice: in this case we will use the "Refined Scaffoled Seq" for gap filling!!!!!!!!!!!!!!!!
//这里是需要修改的，因为我们需要对整体统一进行alignment，然后通过名字，进行相应的解析。
void ClsScaffoldFiller::FillScaffoldByPairEndReads(vector<St_ScaffoldUnit>& vScaffoldSet)
{
    //We need to fill the remained gap based on strRefinedSeq
    /*Work flow:
     * 1: Map the reads to such refined seq
     * 2: try to find:
     *   (1) if some repeats mapped to the new seq which be added as gap filling
     * 3: try to extend the gap until the terminate condition happend
     *   (1) if a reads has
     *      a) one mate in pair full mapped
     *      b) another mate has soft clipped by some sequence before the gap then extend it.
     *   (2) Terminate conditon: if some parts of sequence after the gap could be mapping to certain reads,
     *                           then we stop extension
     */
    //Let's Go!!!!!!!!!!!
    //1: Build fasta file
    ofstream ofs;
    string strFilePath = ClsAlgorithm::GetInstance().GetHigherFolderPath(ClsAlgorithm::GetInstance().GetCurExeFolderPath()) +
                        "TempFile/ScaffoldRefined.fa";
    ofs.open(strFilePath.c_str());
    for(vector<St_ScaffoldUnit>::iterator itr = vScaffoldSet.begin();
        itr != vScaffoldSet.end(); itr++)
    {
        ofs << ">" << itr->strName << endl;
        ofs << itr->strRefinedSeq << endl;
    }
    //build bam file
    string strBamFilePath = ClsAlgorithm::GetInstance().CreateBamFileByMultiScaff(vScaffoldSet,
                                                          m_strReads1Path, m_strReads2Path,
                                                          "RefinedScafSeq", sdtByRepeat);
    //Step 2:Analysis Bam Files for gap filling
    ParseBamFileForGapFill(strBamFilePath, vScaffoldSet);
}

const char* czSCFQuality[] = {"Unknown", "Low", "Normal", "High", "Max"};
void ClsScaffoldFiller::ParseBamFileForGapFill(string strBamFilePath,
                                               vector<St_ScaffoldUnit>& vScaffoldSet)
{
    /*1: 我需要弄清楚，怎样去查看pair-end reads的前后两个mate的比对情况
    //找到下面的两种类型的reads
    //(1) 一个是全匹配，但一个全丢失
          在这个case中： 如果这个全匹配的是
    //(2) 一个全匹配，但另一个半丢失
    //(3) 一个半丢失，另一个也是半丢失
    //(4) 一个全匹配，另一个也全匹配

    虽然我们暂时想不到什么方式去处理这些case，但是我隐约觉得这些case的合理利用，是能够解决gap filling 的问题的

    //2：找到一个方式可以搜集所有的soft -clip的reads
    //3：查看那些被soft-clip的reads，然后找到我们的目标reads
      4: 通过检测depth信息，来看是否新填补的gap是合理的 --->
    */
    //Step 1: 我们先搜集这样的四种类型的reads，以及他们的具体的匹配信息=========>Go!!!
    //我决定试一下chu chong封的两个类，我觉得应该是能够一定作用的 -->Go!!!!!!!!!!!
    //不得不说，还是自己写靠谱！！！！！！！！！！

    BamReader* pBamReader = new BamReader();
    pBamReader->Open(strBamFilePath);
    pBamReader->OpenIndex(strBamFilePath + ".bai");

    //step 1: record all of position of alignment seq
    BamAlignment al;
    vector<int> clipSizes;
    vector<int> readPositions;
    vector<int> genomePositions;
    //----------->使用match的reads数进行该scaffold的质量的好坏，--->这里不需要使用unique map的逻辑，因为scaffold并不是repeat
    int iMatchReadsNum = 0;
    int iSoftClipReadsNum = 0;
    int iScafUnitIndex = 0;
    vector<St_SoftClipReads> vSCReads; //soft clip reads --> they are very useful!!!
    for(vector<St_ScaffoldUnit>::iterator itr = vScaffoldSet.begin(); itr != vScaffoldSet.end(); itr++)
    {
        iMatchReadsNum = 0;
        iSoftClipReadsNum = 0;
        //Clear the old record
        vSCReads.clear();
        clipSizes.clear();
        readPositions.clear();
        genomePositions.clear();
        while(pBamReader->GetNextAlignment(al))  //这个地方可以尝试去定义相应的quality  -->现在这个里面有所有的bam alignment的结果
        {
            if(al.RefID < iScafUnitIndex)
                continue;
            if(al.RefID > iScafUnitIndex)
                break;

            //int iStart = al.Position;
            //if(iStart > 0) //这个时候应该是证明比上了
            //{
                iMatchReadsNum++;
                unsigned int iLastLen = clipSizes.size();
                al.GetSoftClips(clipSizes, readPositions, genomePositions);
                if(clipSizes.size() - iLastLen > 0)
                {
                    int icurClipLen = *(clipSizes.end()-1);
                    int icurClipPos = *(readPositions.end()-1);
                    if(icurClipPos > al.Length)
                        continue;
                    if(icurClipLen == icurClipPos)
                        icurClipPos = 0;
                    int iRefPos = *(genomePositions.end()-1);
                    UpdateBamReads(al, *itr, icurClipLen, icurClipPos, iRefPos, vSCReads);
                    iSoftClipReadsNum++;
                }
            //}
        }
        itr->UpdateQuality(iMatchReadsNum, iSoftClipReadsNum);
        //---->
        cout << "Part Match Num is: " << IntToStr(iMatchReadsNum) << endl;
        cout << "SoftClip Num is: " << IntToStr(iSoftClipReadsNum) << endl;
        cout << "The scaffold Quality is: "  << czSCFQuality[itr->enQuality] << endl;
        //<----
        //use those valid softclip data for gap filling ----->
        FillGapBySoftClipReads(vSCReads);
        iScafUnitIndex++;
    }
    //手动释放vector
    vSCReads.clear();
    clipSizes.clear();
    readPositions.clear();
    genomePositions.clear();
    //release memory
    delete pBamReader;
    pBamReader = NULL;
}

void ClsScaffoldFiller::UpdateBamReads( BamAlignment& al, St_ScaffoldUnit&  stScaffoldUnit,
                                        int icurClipLen, int icurClipPos, int iRefPos,
                                        vector<St_SoftClipReads>& vSCReads)
{
    //Step1: we should check if this this clip is valid
    const int COFFSET = 5;
    for(vector<St_GapRefinedByRept>::iterator itr = stScaffoldUnit.vGapRefinedByRept.begin();
        itr != stScaffoldUnit.vGapRefinedByRept.end(); itr++)
    {
        //use the breakpoint for confirm if such reads could be used
        if((iRefPos-COFFSET > itr->stBamGap.iStartPos && iRefPos-COFFSET < itr->stBamGap.iEndPos) ||
           (iRefPos+COFFSET > itr->stBamGap.iStartPos && iRefPos+COFFSET < itr->stBamGap.iEndPos) )
        {
            //Could be used then-->
            AddValidSCReads(al, *itr, icurClipLen, icurClipPos, iRefPos, vSCReads);
        }
    }
}

//新增soft clip reads
void ClsScaffoldFiller::AddValidSCReads(BamAlignment& al, St_GapRefinedByRept& stGapRefinedByRept,
                                        int icurClipLen, int icurClipPos, int iRefPos,
                                        vector<St_SoftClipReads>& vSCReads)
{
    //Type1: for softclip reads
    St_SoftClipReads stSCReads;
    stSCReads.iRefPos = iRefPos;
    stSCReads.iClipLen = icurClipLen;
    stSCReads.strSeq = al.AlignedBases;//QueryBases; this is the seq which could be algined!!!!!!!!!
    stSCReads.iClipStart = icurClipPos;
    stSCReads.bRevsStrand = al.IsReverseStrand();// IsMateReverseStrand();
    stSCReads.bSecondMate = al.IsSecondMate();
    //this means that which part be cut off for the current reads
    if(icurClipPos == 0) //this tell us : the softclip is start from the beginning
    {
        stSCReads.enClipPart = cpLeft;
        stSCReads.iMissingEnd = stSCReads.iRefPos - 1;
        stSCReads.strClipSeq = al.QueryBases.substr(0, stSCReads.iClipLen);
    }
    else
    {
        stSCReads.enClipPart = cpRight;
        stSCReads.iMissingStart = stSCReads.iRefPos + 1;
        stSCReads.strClipSeq = al.QueryBases.substr(stSCReads.iClipStart, stSCReads.iClipLen);
    }
    //Consider how to update this info -->
    //1: Update SoftClipReads -->
    stSCReads.pGap = &stGapRefinedByRept;
    //2-->do not need to calculate kmer, we could try to get the result by local alignment directly
    //<--
    //Use local Alignment
    LocalAlignment localAl;//--> 使用现有的库来进行Local Alignment
    int iStartOrg = -1, iEndOrg = -1, iStartCmp = -1, iEndCmp = -1;
    switch((int)stSCReads.enClipPart)
    {
        case cpLeft://left part be clipped, means that the right aligned part should be matched with the relevant flank
        {
            // the target is that clip the part that beyond the left breakpoint
            localAl.optAlign(stSCReads.pGap->stBamGap.strRefinedLeftFlank,
                        stSCReads.strClipSeq, iStartOrg, iEndOrg, iStartCmp, iEndCmp);
            //首先应看看长度，只有长度make sense的时候，才进行下面的进一步判断
            //这个if是为了防止clip的片段过短，导致跟相应ref上的flank能够偶然的成功比对上
            if( stSCReads.strClipSeq.length() - iEndCmp <=
                stSCReads.pGap->stBamGap.strRefinedLeftFlank.length() - iEndOrg)
            {
                stSCReads.strExtendSeq = stSCReads.strClipSeq;
            }
            else
            {
                if(iStartCmp <= 1)
                {
                    int iLen = stSCReads.strClipSeq.length() - (iEndCmp + 1);
                    if(iLen > 0) // that means it could be extend
                        stSCReads.strExtendSeq = stSCReads.strClipSeq.substr(iEndCmp+1, iLen);
                }
                else
                    stSCReads.strExtendSeq = stSCReads.strClipSeq;
            }
            break;
        }
        case cpRight:
        {
            // the target is that clip the part that beyond the Right breakpoint
            localAl.optAlign(stSCReads.pGap->stBamGap.strRefinedRightFlank,
                    stSCReads.strClipSeq, iStartOrg, iEndOrg, iStartCmp, iEndCmp);

            //首先应看看长度，只有长度make sense的时候，才进行下面的进一步判断
            //这个if是为了防止clip的片段过短，导致跟相应ref上的flank能够偶然的成功比对上
            if( iStartCmp <= iStartOrg)
            {
                stSCReads.strExtendSeq = stSCReads.strClipSeq;
            }
            else
            {
                if(iEndCmp >= stSCReads.strClipSeq.length()  - 2)
                {
                    int iLen = iStartCmp;
                    if(iLen > 0)
                        stSCReads.strExtendSeq = stSCReads.strClipSeq.substr(0, iLen);
                }
                else
                    stSCReads.strExtendSeq = stSCReads.strClipSeq;
            }
            break;
        }
    }
    if(stSCReads.strExtendSeq == "")
    {
        int i = 0;
    }
    vSCReads.push_back(stSCReads);
}

bool sortfunction (St_SoftClipReads* pStI, St_SoftClipReads* pStJ)
{ return pStI->iMissingStart < pStJ->iMissingStart; } // from small to large

void ClsScaffoldFiller::FillGapBySoftClipReads(vector<St_SoftClipReads>& vSCReads)
{
    if(vSCReads.empty())
        return;    

    map<St_GapRefinedByRept*, vector<St_SoftClipReads*> > mpGapReads;
    for(vector<St_SoftClipReads>::iterator itr = vSCReads.begin();
        itr != vSCReads.end(); itr++)
    {
        mpGapReads[itr->pGap].push_back(&(*itr));
    }
    //Parse the gap
    for(map<St_GapRefinedByRept*, vector<St_SoftClipReads*> >::iterator itr = mpGapReads.begin();
        itr != mpGapReads.end(); itr++)
    {
        //For current gap -->Try to fill it
        const char* Type[2] = {"Left", "Right"};

        //-->For determine if use the reverse complementary as the final result
        //int iNormOre = 0;
        //int iRevComplOre = 0;
        //<--
        vector<string> vExtendFromLeft;
        vector<string> vExtendFromRight;
        int iMaxExtendNumLeft = -1;
        int iMaxExtendNumRight = -1;
        for(vector<St_SoftClipReads*>::iterator itrSCReads = itr->second.begin();
            itrSCReads != itr->second.end(); itrSCReads++)
        {
            //if((*itrSCReads)->bRevsStrand)//我们在这个地方应该允许所有的符合条件的soft clip reads都进来，然后在最后再去定夺是否需要reverse complementary
            //    continue;
            cout << Type[(*itrSCReads)->enClipPart] << ": " << (*itrSCReads)->strExtendSeq << endl;
            switch((int)(*itrSCReads)->enClipPart)
            {
                case cpLeft:
                {
                    //Notice: left clip should be added into right extension
                    vExtendFromRight.push_back((*itrSCReads)->strExtendSeq);
                    int iLen = (*itrSCReads)->strExtendSeq.length();
                    if(iMaxExtendNumRight < iLen)
                        iMaxExtendNumRight = iLen;
                    //--->Statistic the orientation of mapping result
                    //if((*itrSCReads)->bSecondMate)
                    //    iRevComplOre++;
                    //else
                    //    iNormOre++;
                    //<---
                    break;
                }
                case cpRight:
                {
                    //Notice: right clip should be added into left extension
                    vExtendFromLeft.push_back((*itrSCReads)->strExtendSeq);
                    int iLen = (*itrSCReads)->strExtendSeq.length();
                    if(iMaxExtendNumLeft < iLen)
                        iMaxExtendNumLeft = iLen;

                    //--->Statistic the orientation of mapping result
                    //if((*itrSCReads)->bSecondMate)
                    //    iRevComplOre++;
                    //else
                    //    iNormOre++;
                    //<---
                    break;
                }
            }
        }
        //Decide if we need to make the reverse complementay to express the final result
        //bool bRevComplementary = false; //not good
        //if(iRevComplOre > iNormOre)
        //    bRevComplementary = true;

        //Try to Merge the "left" Extend sequence ==========================>
        //relevant with the soft clip reads with RIGHT part
        string strLeftExtend, strRightExtend;
        if(iMaxExtendNumLeft > 0)
        {
            vector<St_CompEtdSeq> vCompLeftSeq;
            vCompLeftSeq.resize(iMaxExtendNumLeft);
            for(vector<string>::iterator itr = vExtendFromLeft.begin();
                itr != vExtendFromLeft.end(); itr++)
            {
                for(int i = 0; i<(int)itr->length(); i++)
                {
                    vCompLeftSeq[i].mpNdCount[itr->at(i)]++;
                }
            }
            //First choose the most frequency cases
            int iSumChar = 0;
            for(vector<St_CompEtdSeq>::iterator itr = vCompLeftSeq.begin();
                itr != vCompLeftSeq.end(); itr++)
            {
                for(map<char, int>::iterator itrmp = itr->mpNdCount.begin();
                    itrmp != itr->mpNdCount.end(); itrmp++)
                {
                    if(itr->iValidCount < itrmp->second)
                    {
                        itr->cValidNd = itrmp->first;
                        itr->iValidCount = itrmp->second;
                    }
                }
                iSumChar += itr->iValidCount;
            }
            int iAverageChar = iSumChar / iMaxExtendNumLeft;
            //Get the final left extend string
            for(vector<St_CompEtdSeq>::iterator itr = vCompLeftSeq.begin();
                itr != vCompLeftSeq.end(); itr++)
            {
                if(itr->iValidCount < iAverageChar)
                    continue;
                strLeftExtend += itr->cValidNd;
            }
        }

        //try to Extend the "right" sequence ====================================>
        if(iMaxExtendNumRight > 0)
        {
            vector<St_CompEtdSeq> vCompRightSeq;
            vCompRightSeq.resize(iMaxExtendNumRight);
            for(vector<string>::iterator itr = vExtendFromRight.begin();
                itr != vExtendFromRight.end(); itr++)
            {
                for(int i = 0; i<(int)itr->length(); i++)
                {                    
                    int iIndex = iMaxExtendNumRight -(itr->length() - i);
                    vCompRightSeq[iIndex].mpNdCount[itr->at(i)]++;
                }
            }
            //First choose the most frequency cases
            int iSumChar = 0;
            for(vector<St_CompEtdSeq>::iterator itr = vCompRightSeq.begin();
                itr != vCompRightSeq.end(); itr++)
            {
                for(map<char, int>::iterator itrmp = itr->mpNdCount.begin();
                    itrmp != itr->mpNdCount.end(); itrmp++)
                {
                    if(itr->iValidCount < itrmp->second)
                    {
                        itr->cValidNd = itrmp->first;
                        itr->iValidCount = itrmp->second;
                    }
                }
                iSumChar += itr->iValidCount;
            }
            int iAverageChar = iSumChar / iMaxExtendNumRight;
            //Get the final left extend string
            for(vector<St_CompEtdSeq>::iterator itr = vCompRightSeq.begin();
                itr != vCompRightSeq.end(); itr++)
            {
                if(itr->iValidCount < iAverageChar)
                    continue;
                strRightExtend += itr->cValidNd;
            }
        }

        //--->新增逻辑: 左右外延的seq存在相应的包含关系的话，那我们需要截断掉边界的seq，保证结果的准确性
        //1：首先选取较长的，然后看较长的是否包含较短的, 防止较长的不正确的延伸
        if(strLeftExtend.length() > strRightExtend.length())
        {
            int iPos = strLeftExtend.find(strRightExtend);
            if(iPos != string::npos)
                strLeftExtend = strLeftExtend.substr(0, iPos) + strRightExtend;
        }
        else if(strLeftExtend.length() < strRightExtend.length())
        {
            int iPos = strRightExtend.find(strLeftExtend);
            if(iPos != string::npos)
                strRightExtend = strRightExtend.substr(iPos, strRightExtend.length() - iPos);
        }
        //<--------

        cout << "Left Extend Sum is: " << strLeftExtend  << endl;
        cout << "Right Extend Sum is: " << strRightExtend << endl;

        /*这里我们新增一个新的算法：
         * 1：根据左右extend的结果，取得原本序列中能尾部/首部重复的相应的长序列
         * 2：结合这种有特征的长序列，跟左右extend的结果综合考虑得到最终的结果
         * 3：这样的好处在于：比如说左延的序列其右端序列刚好跟右延的序列能够完全匹配，这样可以在很大程度上充分说明这个序列是非常的接近缺失序列的
         */
        //1:collect the sequence coincide with the condition (1)
        vector<string> vWholeLeftExt, vWholeRightExt;
        for(vector<St_SoftClipReads*>::iterator itrSCReads = itr->second.begin();
            itrSCReads != itr->second.end(); itrSCReads++)
        {
            switch((int)(*itrSCReads)->enClipPart)
            {
                case cpLeft:
                    if((*itrSCReads)->strExtendSeq.find(strLeftExtend) != string::npos) //证明找到了
                        vWholeRightExt.push_back((*itrSCReads)->strExtendSeq);
                    break;
                case cpRight:
                    if((*itrSCReads)->strExtendSeq.find(strRightExtend) != string::npos) //证明找到了
                        vWholeLeftExt.push_back((*itrSCReads)->strExtendSeq);
                    break;
            }
        }
        // 合并左flank和右flank
        string strExtend;
        if(!vWholeLeftExt.empty() || !vWholeRightExt.empty()) //只要存在不为空那么就一起来搞
        {
            strExtend = CombineExtByWholeExt(strLeftExtend, strRightExtend, vWholeLeftExt, vWholeRightExt);
        }
        else
        {
            //这里现在是存在bug的，也就是一边有一边没有，这种case是需要着重考虑的
            //We need to take consider the following valuable reads
            //(1) The single Left extend whcih contain the right extend
            //(2) The single right extend which contain the left extend
            strExtend = CombineExt(strLeftExtend, strRightExtend);
        }
        //---------->        
        cout << "Additional Extend: " << strExtend << endl;
        string strRevComplment = ClsAlgorithm::GetInstance().GetReverseCompelement(strExtend);
        cout << "Additional Extend(Reverse Complementary): " << strRevComplment << endl;

        itr->first->stBamGap.strFillSeq = strExtend;
        itr->first->stBamGap.strRevCompFillSeq = strRevComplment;
        //Just record for testing
        itr->first->stBamGap.strLeftExt = strLeftExtend;
        itr->first->stBamGap.strRightExt = strRightExtend;
    }
}

/*--->
 * 这里只是做一个记录：
 * 如果在进行refine后，左右有一定数目的碱基的延伸，那么我们对于下一个片段也要进行相应的延伸，
 * 使得数目可以进行相应的比较，从而选取frequency最高的Nucleotide作为最后的selection
 */
string ClsScaffoldFiller::CombineExtByWholeExt(string strLeftExt, string strRightExt,
                                               vector<string>& vWholeLeftExt, vector<string>& vWholeRightExt)
{
    //======>Statistic the candidate for the left Extent
    //我们以right ext 为 anchor 进行定位
    //我们采取的实现方法是，将其补全，然后进行概率的统计和比较
    vector<string> vRefinedLeftExt;
    vRefinedLeftExt.resize(vWholeLeftExt.size());
    int iIndex = 0;
    int iMaxLeftLen = 0;
    for(vector<string>::iterator itr = vWholeLeftExt.begin(); itr != vWholeLeftExt.end(); itr++)
    {
        int iExtIndex = itr->find(strRightExt);
        vRefinedLeftExt[iIndex] = itr->substr(0, iExtIndex);
        if(iMaxLeftLen < iExtIndex)
            iMaxLeftLen = iExtIndex;
        iIndex++;
    }
    //应该在这里更新一下该数组，将原本的left extend 加上去
    //这里有一点需要注意,还是左右包含的case，如果左边和右边的exted存在包含关系，那么我们就没有必要再将：当前的extend 作为case进行考虑了
    LocalAlignment localAl;//--> 使用现有的库来进行Local Alignment
    //两者做local alignment
    int iStart1 = -1, iStart2 = -1, iEnd1 = -1, iEnd2 = -1;
    if(strLeftExt.find(strRightExt) == string::npos) //只有当左边的延伸不包含右边的延伸的时候，我们才将左边的延伸作为一个case考虑进来
    {
        for(vector<string>::iterator itr = vRefinedLeftExt.begin(); itr != vRefinedLeftExt.end(); itr++)
        {
            iStart1 = -1, iStart2 = -1, iEnd1 = -1, iEnd2 = -1;
            localAl.optAlign(*itr, strLeftExt, iStart1, iEnd1, iStart2, iEnd2);
            if(iEnd1 == (int)itr->length())
            {
                string strTemp = strLeftExt.substr(0, iEnd2);
                vRefinedLeftExt.push_back(strTemp);
                if(iMaxLeftLen < (int)strTemp.length())
                    iMaxLeftLen = (int)strTemp.length();
                break;
            }
        }
    }
    //进行相应的合并
    vector<char> vRefinedWholeLeft;
    vector<int> vLeftExtAddCount;
    for(int i=0; i <iMaxLeftLen; i++)
    {
        char cDNA[vRefinedLeftExt.size()];
        int iIndex = 0;
        for(vector<string>::iterator itr = vRefinedLeftExt.begin(); itr != vRefinedLeftExt.end(); itr++)
        {
            if((int)itr->length() > i)
                cDNA[iIndex] = itr->at(itr->length()-i-1);
            else
                cDNA[iIndex] = 'N';
            iIndex++;
        }
        //找到出现频率最高的字符
        int aryCount[4] = {0, 0, 0, 0};
        char aryDNA[4] = {'A', 'T', 'G', 'C'};
        for(int iChar=0; iChar<(int)vRefinedLeftExt.size(); iChar++)
        {
            switch(cDNA[iChar])
            {
                case 'A':
                    aryCount[0]++;
                    break;
                case 'T':
                    aryCount[1]++;
                    break;
                case 'G':
                    aryCount[2]++;
                    break;
                case 'C':
                    aryCount[3]++;
                    break;
                case 'N':
                    break;
            }
        }
        int iMaxIndex = 0;
        int iMaxCount = aryCount[0];
        for(int i=1; i<4; i++)
        {
            if(iMaxCount < aryCount[i])
            {
                iMaxIndex = i;
                iMaxCount = aryCount[i];
            }
        }
        if(iMaxCount != 0) //防止盲目的加入，因为考虑到都为"空"的情况
        {
            vRefinedWholeLeft.insert(vRefinedWholeLeft.begin(), aryDNA[iMaxIndex]);
            vLeftExtAddCount.insert(vLeftExtAddCount.begin(), iMaxCount);
        }
        else
        {
            vRefinedWholeLeft.insert(vRefinedWholeLeft.begin(), NULL);
            vLeftExtAddCount.insert(vLeftExtAddCount.begin(), -1);
        }
    }
    //&&&&&&&&&&--> 额外代码获取左边的延伸
    int iLeftExtAddNum = (int)vRefinedWholeLeft.size(); // additional left extend number
    vector<string> vRightExtAddLeftSeq;
    vRightExtAddLeftSeq.resize(vWholeRightExt.size());
    iIndex = 0;
    int iMaxRightExtAddLeftLen = 0;
    for(vector<string>::iterator itr = vWholeRightExt.begin(); itr != vWholeRightExt.end(); itr++)
    {
        int iExtIndex = itr->find(strLeftExt);
        int iRefinedLen = iLeftExtAddNum > iExtIndex ? iExtIndex : iLeftExtAddNum;
        vRightExtAddLeftSeq[iIndex] = itr->substr(iExtIndex - iRefinedLen, iRefinedLen);
        if(iMaxRightExtAddLeftLen < iRefinedLen)
            iMaxRightExtAddLeftLen = iRefinedLen;
        iIndex++;
    }
    //<--
    //得到之后，我们计算每个位置的词频率，然后选取最大的并记录下次数。==>
    vector<int> vRightExtAddLeftCount;
    vector<char> vRightExtAddLeft;
    for(int i=0; i <iMaxRightExtAddLeftLen; i++)
    {
        char cDNA[vRightExtAddLeftSeq.size()];
        int iIndex = 0;
        for(vector<string>::iterator itr = vRightExtAddLeftSeq.begin();
            itr != vRightExtAddLeftSeq.end(); itr++)
        {
            if((int)itr->length() > i)
                cDNA[iIndex] = itr->at(itr->length()-i-1);
            else
                cDNA[iIndex] = 'N';
            iIndex++;
        }
        //找到出现频率最高的字符
        int aryCount[4] = {0, 0, 0, 0};
        char aryDNA[4] = {'A', 'T', 'G', 'C'};
        for(int iChar=0; iChar<(int)vRightExtAddLeftSeq.size(); iChar++)
        {
            switch(cDNA[iChar])
            {
                case 'A':
                    aryCount[0]++;
                    break;
                case 'T':
                    aryCount[1]++;
                    break;
                case 'G':
                    aryCount[2]++;
                    break;
                case 'C':
                    aryCount[3]++;
                    break;
                case 'N':
                    break;
            }
        }
        int iMaxIndex = 0;
        int iMaxCount = aryCount[0];
        for(int i=1; i<4; i++)
        {
            if(iMaxCount < aryCount[i])
            {
                iMaxIndex = i;
                iMaxCount = aryCount[i];
            }
        }
        if(iMaxCount != 0) //防止盲目的加入，因为考虑到都为"空"的情况
        {
            vRightExtAddLeft.insert(vRightExtAddLeft.begin(), aryDNA[iMaxIndex]);
            vRightExtAddLeftCount.insert(vRightExtAddLeftCount.begin(), iMaxCount);
        }
        else
        {
            vRightExtAddLeft.insert(vRightExtAddLeft.begin(), NULL);
            vRightExtAddLeftCount.insert(vRightExtAddLeftCount.begin(), -1);
        }
    }
    //然后跟之前的 strRefinedWholeLeft进行比较
    string strRefinedWholeLeft;
    int iOffSet = iLeftExtAddNum - iMaxRightExtAddLeftLen;
    if(iOffSet > 0)
    {
        for(int i=0; i<iOffSet; i++)
            strRefinedWholeLeft += vRefinedWholeLeft[i];
    }
    for(int i=0; i<iMaxRightExtAddLeftLen; i++)
    {
        if(vRightExtAddLeftCount[i] > vLeftExtAddCount[i])
            strRefinedWholeLeft += vRightExtAddLeft[i];
        else
            strRefinedWholeLeft += vRefinedWholeLeft[iOffSet + i];
    }
    //<==
    strRefinedWholeLeft = strRefinedWholeLeft + strRightExt;

    //=========>Statistic the candidate for the right Extent
    //--> *右边的扩展(由右向左的扩展): 这种多余的ext的比较只针对右边的extend进行补充
    //核心思想就是，除了正常的外，还要多扩展几个
    vector<string> vRefinedRightExt;
    vRefinedRightExt.resize(vWholeRightExt.size());
    iIndex = 0;
    int iMaxRightLen = 0;
    for(vector<string>::iterator itr = vWholeRightExt.begin(); itr != vWholeRightExt.end(); itr++)
    {
        int iExtIndex = itr->find(strLeftExt);
        int iRefinedLen = itr->length()-(iExtIndex + strLeftExt.length()-1);
        vRefinedRightExt[iIndex] = itr->substr(iExtIndex + strLeftExt.length(), iRefinedLen);
        if(iMaxRightLen < iRefinedLen)
            iMaxRightLen = iRefinedLen;
        iIndex++;
    }

    //应该在这里更新一下该数组，将原本的left extend 加上去
    //两者做local alignment
    //只有当右边的延伸不包左边的延伸的时候，我们才将右边的延伸作为一个case考虑进来
    if(strRightExt.find(strLeftExt) == string::npos)
    {
        for(vector<string>::iterator itr = vRefinedRightExt.begin(); itr != vRefinedRightExt.end(); itr++)
        {
            iStart1 = -1, iStart2 = -1, iEnd1 = -1, iEnd2 = -1;
            localAl.optAlign(*itr, strRightExt, iStart1, iEnd1, iStart2, iEnd2);
            if(iStart1 == 1)//这个返回值是从1开始，因为不要最后续的计算，所以在这里不需要将它们都统一“--”
            {
                string strTemp = strRightExt.substr(iStart2 - 1, strRightExt.length() - iStart2 + 1);
                vRefinedRightExt.push_back(strTemp);
                if(iMaxRightLen < (int)strTemp.length())
                    iMaxRightLen = (int)strTemp.length();
                break;
            }
        }
    }

    string strRefinedWholeRight = "";
    for(int i=0; i <iMaxRightLen; i++)
    {
        char cDNA[vRefinedRightExt.size()];
        int iIndex = 0;
        for(vector<string>::iterator itr = vRefinedRightExt.begin(); itr != vRefinedRightExt.end(); itr++)
        {
            if(i < (int)itr->length())
                cDNA[iIndex] = itr->at(i);
            else
                cDNA[iIndex] = 'N';
            iIndex++;
        }
        //找到出现频率最高的字符
        int aryCount[4] = {0, 0, 0, 0};
        char aryDNA[4] = {'A', 'T', 'G', 'C'};
        for(int iChar=0; iChar<(int)vRefinedRightExt.size(); iChar++)
        {
            switch(cDNA[iChar])
            {
                case 'A':
                    aryCount[0]++;
                    break;
                case 'T':
                    aryCount[1]++;
                    break;
                case 'G':
                    aryCount[2]++;
                    break;
                case 'C':
                    aryCount[3]++;
                    break;
                case 'N':
                    break;
            }
        }
        int iMaxIndex = 0;
        int iMaxCount = aryCount[0];
        for(int i=1; i<4; i++)
        {
            if(iMaxCount < aryCount[i])
            {
                iMaxIndex = i;
                iMaxCount = aryCount[i];
            }
        }
        if(iMaxCount != 0)
            strRefinedWholeRight = strRefinedWholeRight + aryDNA[iMaxIndex];
    }
    strRefinedWholeRight = strLeftExt + strRefinedWholeRight;

    //最后还要跟原来的Right Extend进行合并矫正一下
    //Make local alignment
    //两者做local alignment 然后去交集
    iStart1 = -1, iStart2 = -1, iEnd1 = -1, iEnd2 = -1;
    string strTempLeft, strTempRight;
    if(vWholeLeftExt.empty()) // 只有右边的可以延伸
    {
        strTempLeft = strLeftExt;
        strTempRight = strRefinedWholeRight;
    }
    else if(vWholeRightExt.empty()) // 只有左边的可以延伸
    {
        strTempLeft = strRefinedWholeLeft;
        strTempRight = strRightExt;
    }
    else //两边都可以延伸
    {
        strTempLeft = strRefinedWholeLeft;
        strTempRight = strRefinedWholeRight;
    }
    //local alignment
    localAl.optAlign(strTempLeft, strTempRight, iStart1, iEnd1, iStart2, iEnd2);
    //因为起始位置是1，所以上的所有start 和 end的值需要依次去递减一个：
    iStart1--, iEnd1--, iStart2--, iEnd2--;
    //分析结果，然后取得交集
    string strResult = "";
    if(iStart2 == 0 &&
       iStart1 >=0 && iEnd1 >=0 && iEnd2 >=0)
    {
        if(iStart1 != 0)
            strResult = strTempLeft.substr(0, iStart1) + strTempRight;
        else if(iEnd1 == (int)strTempLeft.length() - 1)
            strResult = strTempLeft + strTempRight.substr(iEnd2 + 1, strTempRight.length() - iEnd2);
    }
    else
        return CombineExt(strLeftExt, strRightExt);
    return strResult;
}

//this is a recursive
string ClsScaffoldFiller::CombineExt(string strExt1, string strExt2/*, int iStart1, int iStart2, int iEnd1, int iEnd2*/)
{
    //Step1: 首先判断两边有一者为空的情况，因为这种情况是需要去作local alignment的（强行做会使得得到的start和end的值异常大，从而导致后面的计算出错）
    // 实际上如果存在一遍的数目为空，实际上这个case是不应该常发生的。因此，这个可以在一定程度上认定，这个是属于unreliable的
    if(strExt1 == "" || strExt2 == "")
        return strExt1 + strExt2;

    //我们可以考虑一下-->global alignment --> 这样说不定会有比较好的结果 -->今天到此为止吧 --->一会儿去看看python去，完成一下那个老师的作业
    LocalAlignment localAl;//--> 使用现有的库来进行Local Alignment
    int iStart1 = -1, iStart2 = -1, iEnd1 = -1, iEnd2 = -1;
    localAl.optAlign(strExt1, strExt2, iStart1, iEnd1, iStart2, iEnd2);

    En_AlnState enExt1Side, enExt2Side;
    //case 0: 有一个是全匹配，那么纠结，直接最长的那个就是结果
    //For enExt1Side  // 这个的返回值是按照1为start的，并不是0 这个要注意！！！！！， 我么可以把他们搞成0；
    iStart1--;
    iEnd1--;
    iStart2--;
    iEnd2--;
    if( (iStart1 == 0 && iEnd1 == (int)strExt1.length() - 1) ||
        (iStart2 == 0 && iEnd2 == (int)strExt2.length() - 1))
    {
        return strExt1.length() >= strExt2.length() ? strExt1 : strExt2;
    }

    //Make sure the statement of cut side
    //For enExt1Side    
    if(iStart1 > 0 && iEnd1 < (int)strExt1.length() - 1)
        enExt1Side = apMiddle;
    else if(iStart1 > 0 && iEnd1 >= (int)strExt1.length() - 1)
        enExt1Side = apRight;
    else if(iStart1 <= 0 && iEnd1 <= (int)strExt1.length() - 1)
        enExt1Side = apLeft;
    else if(iStart1 == 0 && iEnd1 == (int)strExt1.length() - 1)
        enExt1Side =apWhole;
    else {}

    //For enExt2Side
    if(iStart2 > 0 && iEnd2 < (int)strExt2.length() - 1)
        enExt2Side = apMiddle;
    else if(iStart2 > 0 && iEnd2 == (int)strExt2.length() - 1)
        enExt2Side = apRight;
    else if(iStart2 == 0 && iEnd2 < (int)strExt2.length() - 1)
        enExt2Side = apLeft;
    else if(iStart2 == 0 && iEnd2 == (int)strExt2.length() - 1)
        enExt2Side = apWhole;
    else {}

    //Let's Compare
    string strSubExt1, strCmbExt1, strCmbExt1_1, strCmbExt1_2,
           strSubExt2, strCmbExt2, strCmbExt2_1, strCmbExt2_2;
    switch((int)enExt1Side)
    {
        case apLeft:
        {
            strSubExt1 = strExt1.substr(iStart1, iEnd1 - iStart1 + 1);
            strCmbExt1 = strExt1.substr(iEnd1 + 1, strExt1.length() - iEnd1);
            switch((int)enExt2Side)
            {
                case apLeft:
                {
                    //sub：意味着比对成功的 ， cmb：意味着需要继续考虑的
                    strSubExt2 = strExt2.substr(iStart2, iEnd2 - iStart2 + 1);
                    strCmbExt2 = strExt2.substr(iEnd2+1, strExt2.length() - iEnd2);
                    return strExt2;
                }
                case apRight:
                {
                    strSubExt2 = strExt2.substr(iStart2, iEnd2 - iStart2 + 1);
                    strCmbExt2 = strExt2.substr(0, iStart2);
                    return strExt2;
                }
                case apMiddle:
                {
                    strSubExt2 = strExt2.substr(iStart2, iEnd2 - iStart2 + 1);
                    strCmbExt2_1 = strExt2.substr(0, iStart2);
                    strCmbExt2_2 = strExt2.substr(iEnd2+1, strExt2.length() - iEnd2);
                    return strSubExt1 + strCmbExt2_2;
                    //do not consider                    
                }
                case apWhole:
                {
                    return strExt2;
                }
            }
            break;
        }
        case apRight:
        {
            strSubExt1 = strExt1.substr(iStart1, iEnd1 - iStart1 + 1);
            strCmbExt1 = strExt1.substr(0, iStart1);
            switch((int)enExt2Side)
            {
                case apLeft:
                {
                    strSubExt2 = strExt2.substr(iStart2, iEnd2 - iStart2 + 1);
                    strCmbExt2 = strExt2.substr(iEnd2+1, strExt2.length() - iEnd2);
                    return strCmbExt1 + strSubExt1 + strCmbExt2; // 也就是 strExt1
                }
                case apRight:
                {
                    strSubExt2 = strExt2.substr(iStart2, iEnd2 - iStart2 + 1);
                    strCmbExt2 = strExt2.substr(0, iStart2);
                    return strExt1;
                }
                case apMiddle:
                {
                    strSubExt2 = strExt2.substr(iStart2, iEnd2 - iStart2 + 1);
                    strCmbExt2_1 = strExt2.substr(0, iStart2);
                    strCmbExt2_2 = strExt2.substr(iEnd2+1, strExt2.length() - iEnd2);
                    return strExt1 + strCmbExt2_2;
                }
                case apWhole:
                {
                    return strExt1;
                }
            }
            break;
        }
        case apMiddle:
        {
            strSubExt1 = strExt1.substr(iStart1, iEnd1 - iStart1 + 1);
            strCmbExt1_1 = strExt1.substr(0, iStart1);
            strCmbExt1_2 = strExt1.substr(iEnd1+1, strExt1.length() - iEnd2);
            switch((int)enExt2Side)
            {
                case apLeft:
                {
                    strSubExt2 = strExt2.substr(iStart2, iEnd2 - iStart2 + 1);
                    strCmbExt2 = strExt2.substr(iEnd2+1, strExt2.length() - iEnd2);
                    return strCmbExt1_1 + strSubExt1 + strCmbExt2; // 也就是 strExt1
                }
                case apRight:
                {
                    strSubExt2 = strExt2.substr(iStart2, iEnd2 - iStart2 + 1);
                    strCmbExt2 = strExt2.substr(0, iStart2);
                    return strCmbExt1_1 + strSubExt1;
                }
                case apMiddle:
                {
                    strSubExt2 = strExt2.substr(iStart2, iEnd2 - iStart2 + 1);
                    strCmbExt2_1 = strExt2.substr(0, iStart2);
                    strCmbExt2_2 = strExt2.substr(iEnd2+1, strExt2.length() - iEnd2);
                    return strCmbExt1_1 + strSubExt1 + strCmbExt2_2;
                }
                case apWhole:
                {
                    return strCmbExt1_1 + strCmbExt1_1;
                }
            }
            break;
        }
        case apWhole:
        {            
            switch((int)enExt2Side)
            {
                case apLeft:
                {
                    strSubExt2 = strExt2.substr(iStart2, iEnd2 - iStart2 + 1);
                    strCmbExt2 = strExt2.substr(iEnd2+1, strExt2.length() - iEnd2);
                    return strExt2; // 也就是 strExt1
                }
                case apRight:
                {
                    strSubExt2 = strExt2.substr(iStart2, iEnd2 - iStart2 + 1);
                    strCmbExt2 = strExt2.substr(0, iStart2);
                    return strExt1;
                }
                case apMiddle:
                {
                    strSubExt2 = strExt2.substr(iStart2, iEnd2 - iStart2 + 1);
                    strCmbExt2_1 = strExt2.substr(0, iStart2);
                    strCmbExt2_2 = strExt2.substr(iEnd2+1, strExt2.length() - iEnd2);
                    return strExt1 + strCmbExt2_2;
                }
                case apWhole:
                {
                    return strExt1;
                }
            }
            break;
        }
    }
    return "";
    /*
    strExtend = strLeftExtend; //Init Value;

    //这里的几个extend还需要跟上面的original的去比较才合适
    if(iStartCmp > 1) //== Expand from left
        strExtend = strRightExtend.substr(0, iStartCmp-1) + strExtend;
    if(iEndCmp < strRightExtend.length())
        strExtend = strExtend + strRightExtend.substr(iEndCmp, strRightExtend.length() - iEndCmp);
    */
}

void ClsScaffoldFiller::CheckFillQuality(St_ScaffoldUnit& stScaffoldUnit)
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
    m_vFinalResult.push_back(stFinalScaffoldUnit);
}

void ClsScaffoldFiller::QualityEstimate(string strBamFilePath,
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
            int iEnd2 = itr->iEnd + iLen /
                    2;

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
