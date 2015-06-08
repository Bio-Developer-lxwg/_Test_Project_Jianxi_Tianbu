#include "clspairendfill.h"
#include "clsalgorithm.h"
#include "local_alignment.h"
#include "smith-waterman.h"
#include "algorithm" //This is the system algorithm such as: sort, etc
#include "math.h"

ClsPairEndFill::ClsPairEndFill()
{
    m_strReads1Path = "";
    m_strReads2Path = "";
}

ClsPairEndFill::~ClsPairEndFill()
{
}

void ClsPairEndFill::BaseInit()
{}

void ClsPairEndFill::Init(string strReads1Path, string strReads2Path)
{
    BaseInit();
    m_strReads1Path = strReads1Path;
    m_strReads2Path = strReads2Path;
}

const char* czSCFQuality[] = {"Unknown", "Low", "Normal", "High", "Max"};
void ClsPairEndFill::FillScafUnit(St_ScaffoldUnit& stScafUnit)
{
    //build bam file by CURRENT SCAFFOLD UNIT !
    string strBamFilePath = ClsAlgorithm::GetInstance().CreateBamFile(stScafUnit.strName, stScafUnit.strRefinedSeq,
                                                                      m_strReads1Path, m_strReads2Path);
    //Parse bam file to fill the gap
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
    vector<St_SoftClipReads> vSCReads; //soft clip reads --> they are very useful!!!
    while(pBamReader->GetNextAlignment(al))  //这个地方可以尝试去定义相应的quality  -->现在这个里面有所有的bam alignment的结果
    {
        if(!al.IsMapped()) // if do not map --> break out
            continue;
        iMatchReadsNum++;

        if(al.Position > 1900)
        {
            int i = 0;
        }

        unsigned int iLastLen = clipSizes.size();
        al.GetSoftClips(clipSizes, readPositions, genomePositions);
        if(clipSizes.size() - iLastLen > 0)
        {
            int icurClipLen = *(clipSizes.end()-1);
            int icurClipPos = *(readPositions.end()-1);
            if(icurClipPos > al.Length) //In most of cases, this case will not happen
                continue;
            if(icurClipLen == icurClipPos)
                icurClipPos = 0;
            int iRefPos = *(genomePositions.end()-1);
            UpdateBamReads(al, stScafUnit, icurClipLen, icurClipPos, iRefPos, vSCReads);
            iSoftClipReadsNum++;
        }
        /*
        else
        {
            int iRefStartPos = al.Position;
            int iRefEndPos = al.GetEndPosition();
            int iReadlLen = 100;
            enum En_CigarType{ctM = 0, ctI, ctD, ctN, ctS, ctMax};
            int arrCount[ctMax];
            memset(arrCount, 0, sizeof(int)*ctMax);
            //There will be no soft clip if just small number of character be matched
            //So let's check the mapping position
            for(vector<St_GapRefinedByRept>::iterator itr = stScafUnit.vGapRefinedByRept.begin();
                itr != stScafUnit.vGapRefinedByRept.end(); itr++)
            {
                //use the breakpoint for confirm if such reads could be used
                if((iRefStartPos < itr->stBamGap.iStartPos && iRefStartPos + iReadlLen > itr->stBamGap.iStartPos)||
                   (iRefEndPos > itr->stBamGap.iEndPos && iRefStartPos + iReadlLen < itr->stBamGap.iEndPos))
                {
                    //Check Cigar
                    for(vector<CigarOp>::iterator itr = al.CigarData.begin(); itr != al.CigarData.end(); itr++)
                    {
                        switch(itr->Type)
                        {
                            case 'M':
                                arrCount[ctM]++;
                                break;
                            case 'I':
                                arrCount[ctI]++;
                                break;
                            case 'D':
                                arrCount[ctD]++;
                                break;
                            case 'N':
                                arrCount[ctN]++;
                                break;
                            case 'S':
                                arrCount[ctS]++;
                                break;
                        }
                    }
                    //Update Data:
                    //AddValidSCReads(al, *itr, icurClipLen, icurClipPos, iRefPos, vSCReads);
                }
            }
        }*/
    }
    stScafUnit.UpdateQuality(iMatchReadsNum, iSoftClipReadsNum);
    //---->
    cout << "Part Match Num is: " << IntToStr(iMatchReadsNum) << endl;
    cout << "SoftClip Num is: " << IntToStr(iSoftClipReadsNum) << endl;
    cout << "The scaffold Quality is: "  << czSCFQuality[stScafUnit.enQuality] << endl;
    //<----
    //use those valid softclip data for gap filling ----->
    FillGapBySoftClipReads(vSCReads);
    //release memory
    delete pBamReader;
    pBamReader = NULL;
}

//Step 2: Notice: in this case we will use the "Refined Scaffoled Seq" for gap filling!!!!!!!!!!!!!!!!
//这里是需要修改的，因为我们需要对整体统一进行alignment，然后通过名字，进行相应的解析。
//I got a big mistake here! --->More than one gap will be abtained in Current Scaffold.
//We need to classify those paire end redas whihc relvant with current gap !!!!!!!!!!
void ClsPairEndFill::FillScafSet(vector<St_ScaffoldUnit>& vScafSet)
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
    for(vector<St_ScaffoldUnit>::iterator itr = vScafSet.begin();
        itr != vScafSet.end(); itr++)
    {
        ofs << ">" << itr->strName << endl;
        ofs << itr->strRefinedSeq << endl;
    }
    ofs.close();
    //build bam file
    string strBamFilePath = ClsAlgorithm::GetInstance().CreateBamFileByMultiScaff(vScafSet,
                                                          m_strReads1Path, m_strReads2Path,
                                                          "RefinedScafSeq", sdtByRepeat);
    //Step 2:Analysis Bam Files for gap filling
    ParseBamFileForGapFill(strBamFilePath, vScafSet);
}

void ClsPairEndFill::ParseBamFileForGapFill(string strBamFilePath,
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
                    if(icurClipPos > al.Length) //It seems that it is impossible
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

void ClsPairEndFill::UpdateBamReads(BamAlignment& al, St_ScaffoldUnit&  stScaffoldUnit,
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

void ClsPairEndFill::FillGapBySoftClipReads(vector<St_SoftClipReads>& vSCReads)
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
        //ofstream ofs;
        //ofs.open("/home/lq/Desktop/scaff_output.ini");
        for(vector<St_SoftClipReads*>::iterator itrSCReads = itr->second.begin();
            itrSCReads != itr->second.end(); itrSCReads++)
        {
            //if((*itrSCReads)->bRevsStrand)//我们在这个地方应该允许所有的符合条件的soft clip reads都进来，然后在最后再去定夺是否需要reverse complementary
            //    continue;
            cout << Type[(*itrSCReads)->enClipPart] << ": " << (*itrSCReads)->strExtendSeq << endl;
            //ofs << Type[(*itrSCReads)->enClipPart] << ": " << (*itrSCReads)->strExtendSeq << endl;
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
        //ofs.close();
        //Decide if we need to make the reverse complementay to express the final result
        //bool bRevComplementary = false; //not good
        //if(iRevComplOre > iNormOre)
        //    bRevComplementary = true;

        //Try to Merge the "left" Extend sequence ==========================>
        //relevant with the soft clip reads with RIGHT part
        string strLeftExtend, strRightExtend;
        //In this phase ,we need to consider two diff cases--->
        //再此处我们将已经被repeats成功填补的case和没有被repeat填补的case分开来进行考虑，因为他们的情况是非常不一样的。
        //已经被repeat修补过的case相对而言要简单很多（可能吧）
        if(itr->first->pOrgGap->bFilled) //For the gap wich has been filled by repeat
        {
            //Consider the type of current filled gap
            const int cKmerCmpNum = 6;
            if(itr->first->pOrgGap->enGapType == gtRight) // this is the case of current
            {
                /*在这个实现里面，有以下几点表现
                 * 1：截掉长的左边的，留下长的右边的
                 * 2：去掉剩下的过长的
                 * 3：去掉剩下的过短的
                 * 4：整体右对齐
                 * 5：取频率最高的character作为最后的结果
                 */
                //至少我今天觉得，在这个地方，我们应该首先将数据及过滤一遍，然后再按照现在的去作
                //1: get 8-mer to check if there is the relationship of containing
                string strKmer = "";
                if(itr->first->pOrgGap->strFillSeq.length() > cKmerCmpNum) //最后6个，而不是前面6个
                    strKmer = itr->first->pOrgGap->strFillSeq.substr(itr->first->pOrgGap->strFillSeq.length()-cKmerCmpNum,
                                                                     cKmerCmpNum);
                else
                    strKmer = itr->first->pOrgGap->strFillSeq;
                //2: traverse the data to refine the data set
                bool bAllPass = false; // If it could be whole filled
                for(vector<string>::iterator itr = vExtendFromRight.begin();
                    itr != vExtendFromRight.end(); itr++)
                {
                    if(itr->find(strKmer.c_str()) != string::npos)
                    {
                        int iPosStart = itr->find(strKmer.c_str()) + strKmer.length();
                        *itr = itr->substr(iPosStart, itr->length()-iPosStart);
                        bAllPass = true;
                    }
                }
                if(bAllPass)
                {
                    itr->first->pOrgGap->enFillResult = frAllPass; // 因为右边跟左边存在右边包含左边的关系，这就说明了此gap是能够被whole fill的
                }
                //erase the null items
                for(vector<string>::iterator itr = vExtendFromRight.end() - 1;
                    itr >= vExtendFromRight.begin(); itr--)
                {
                    if(*itr == "")
                        vExtendFromRight.erase(itr);
                }
                //sort from small to large
                sort(vExtendFromRight.begin(), vExtendFromRight.end(), sortfunctionLength);
                int iDelNum= vExtendFromRight.size() * .1;
                //filter1: the first minium
                for(int i=0; i<iDelNum; i++) //del the first several
                {
                    vExtendFromRight.erase(vExtendFromRight.begin());
                }
                //filter2: the last maxiumn
                for(int i=0; i<iDelNum; i++) //del the last several
                {
                    vExtendFromRight.erase(vExtendFromRight.end()-1);
                }
                //Filter 3: calc the average length and filter the one which longer than 2 times of average
                int iSumLen = 0;
                for(vector<string>::iterator itr = vExtendFromRight.begin();
                    itr != vExtendFromRight.end(); itr++)
                {
                    iSumLen += itr->length();
                }
                size_t iLenThreshold = floor((float)iSumLen / vExtendFromRight.size()) * 2;
                //因为这里已经排过序了，因此我们发现有小于iLenThreshold的就可以跳出循环了
                for(vector<string>::iterator itr = vExtendFromRight.end() - 1;
                    itr >= vExtendFromRight.begin(); itr--)
                {
                    if(itr->length() > iLenThreshold)
                        vExtendFromRight.erase(itr);
                    else
                        break;
                }
                //Function: Get the most frequency Char --->
                //Get the maxiumn length
                iMaxExtendNumRight = 0;
                for(vector<string>::iterator itr = vExtendFromRight.begin();
                    itr != vExtendFromRight.end(); itr++)
                {
                    if(iMaxExtendNumRight < (int)itr->length())
                        iMaxExtendNumRight = itr->length();
                }
                //get the most frequent char
                vector<char> vMostFreqt;
                vector<int> vCount;
                ClsAlgorithm::GetInstance().GetMostFreqChar(vExtendFromRight, vMostFreqt, vCount,
                                                            iMaxExtendNumRight, atRight); //右对齐
                //<---
                //Set value
                for(vector<char>::iterator itr = vMostFreqt.begin(); itr != vMostFreqt.end(); itr++)
                    strRightExtend += *itr;
                strLeftExtend = "";
            }
            else if(itr->first->pOrgGap->enGapType == gtLeft)
            {
                //just consider right softclip
                /*相对于gtright而言，在这个实现里面，有以下几点表现
                 * 1：截掉长的右边的，留下长的左边的
                 * 2：去掉剩下的过长的
                 * 3：去掉剩下的过短的
                 * 4：整体左对齐
                 * 5：取频率最高的character作为最后的结果
                 */
                string strKmer = "";
                if(itr->first->pOrgGap->strFillSeq.length() > cKmerCmpNum) //最后8个，而不是前面8个
                    strKmer = itr->first->pOrgGap->strFillSeq.substr(0, cKmerCmpNum);
                else
                    strKmer = itr->first->pOrgGap->strFillSeq;
                //2: traverse the data to refine the data set
                bool bAllPass = false; // If it could be whole filled
                for(vector<string>::iterator itr = vExtendFromLeft.begin();
                    itr != vExtendFromLeft.end(); itr++)
                {
                    if(itr->find(strKmer.c_str()) != string::npos)
                    {
                        int iPosStart = itr->find(strKmer.c_str()) + strKmer.length();
                        *itr = itr->substr(0, iPosStart); //取前面的，也就是从begin到该kemer的区间内的字段
                        bAllPass = true;
                    }
                }
                if(bAllPass)
                {
                    itr->first->pOrgGap->enFillResult = frAllPass;
                }
                //erase the null items
                for(vector<string>::iterator itr = vExtendFromLeft.end() - 1;
                    itr >= vExtendFromLeft.begin(); itr--)
                {
                    if(*itr == "")
                        vExtendFromLeft.erase(itr);
                }
                //sort from small to large
                sort(vExtendFromLeft.begin(), vExtendFromLeft.end(), sortfunctionLength);
                int iDelNum= vExtendFromLeft.size() * .1;
                //filter1: the first minium
                for(int i=0; i<iDelNum; i++) //del the first several
                {
                    vExtendFromLeft.erase(vExtendFromLeft.begin());
                }
                //filter2: the last maxiumn
                for(int i=0; i<iDelNum; i++) //del the last several
                {
                    vExtendFromLeft.erase(vExtendFromLeft.end()-1);
                }
                //Filter 3: calc the average length and filter the one which longer than 2 times of average
                int iSumLen = 0;
                for(vector<string>::iterator itr = vExtendFromLeft.begin();
                    itr != vExtendFromLeft.end(); itr++)
                {
                    iSumLen += itr->length();
                }
                int iLenThreshold = floor((float)iSumLen / vExtendFromLeft.size()) * 2;
                //因为这里已经排过序了，因此我们发现有小于iLenThreshold的就可以跳出循环了
                for(vector<string>::iterator itr = vExtendFromLeft.end() - 1;
                    itr >= vExtendFromLeft.begin(); itr--)
                {
                    if((int)itr->length() > iLenThreshold)
                        vExtendFromRight.erase(itr);
                    else
                        break;
                }
                //Get the maxiumn length
                iMaxExtendNumRight = 0;
                for(vector<string>::iterator itr = vExtendFromLeft.begin();
                    itr != vExtendFromLeft.end(); itr++)
                {
                    if(iMaxExtendNumRight < (int)itr->length())
                        iMaxExtendNumRight = (int)itr->length();
                }
                //get the most frequent char
                vector<char> vMostFreqt;
                vector<int> vCount;
                ClsAlgorithm::GetInstance().GetMostFreqChar(vExtendFromLeft, vMostFreqt, vCount,
                                                            iMaxExtendNumRight, atLeft); //左对齐
                //Set value
                for(vector<char>::iterator itr = vMostFreqt.begin(); itr != vMostFreqt.end(); itr++)
                   strLeftExtend += *itr;
                strRightExtend = "";
            }
            else // The case of center
            {
                strLeftExtend = "";
                strRightExtend = "";
            }
        }
        else //Do not repaired by repeats
        {
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
                int iAverageChar = iSumChar / iMaxExtendNumLeft * .6;
                //Get the final left extend string
                for(vector<St_CompEtdSeq>::iterator itr = vCompLeftSeq.begin();
                    itr != vCompLeftSeq.end(); itr++)
                {
                    if(itr->iValidCount < iAverageChar)
                        break;//continue;
                    strLeftExtend += itr->cValidNd;
                }
            }

            //try to Extend the "right" sequence ====================================>实际上是位于右边的延伸，实际上是从右往左进行延伸
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
                for(vector<St_CompEtdSeq>::iterator itr = vCompRightSeq.end() - 1;
                    itr >= vCompRightSeq.begin(); itr--)
                {
                    if(itr->iValidCount < iAverageChar)
                        break;
                    strRightExtend = itr->cValidNd + strRightExtend;
                }
                /*
                for(vector<St_CompEtdSeq>::iterator itr = vCompRightSeq.begin();
                    itr != vCompRightSeq.end(); itr++)
                {
                    if(itr->iValidCount < iAverageChar)
                        continue;
                    strRightExtend += itr->cValidNd;
                }*/
            }

            //--->新增逻辑: 左右外延的seq存在相应的包含关系的话，那我们需要截断掉边界的seq，保证结果的准确性
            //1：首先选取较长的，然后看较长的是否包含较短的, 防止较长的不正确的延伸
            bool bAllPass = false;
            if(strLeftExtend.length() > strRightExtend.length())
            {
                size_t iPos = strLeftExtend.find(strRightExtend);
                if(iPos != string::npos)
                {
                    strLeftExtend = strLeftExtend.substr(0, iPos) + strRightExtend;
                    bAllPass = true;
                }
            }
            else if(strLeftExtend.length() < strRightExtend.length())
            {
                size_t iPos = strRightExtend.find(strLeftExtend);
                if(iPos != string::npos)
                {
                    strRightExtend = strRightExtend.substr(iPos, strRightExtend.length() - iPos);
                    bAllPass = true;
                }
            }
            if(bAllPass)
            {
                itr->first->pOrgGap->enFillResult = frAllPass;
            }
            //<--------
        }

        cout << "Left Extend Sum is: " << strLeftExtend  << endl;
        cout << "Right Extend Sum is: " << strRightExtend << endl;

        vector<string> vWholeLeftExt, vWholeRightExt;
        if(!itr->first->pOrgGap->bFilled) // 当没有被repeat填充的时候，我们才触发这个补充算法
        {
            //下面的这个补充算法，应该是对于不没有被repeat填充的case才有作用

            /*这里我们新增一个新的算法：
             * 1：根据左右extend的结果，取得原本序列中能尾部/首部重复的相应的长序列
             * 2：结合这种有特征的长序列，跟左右extend的结果综合考虑得到最终的结果
             * 3：这样的好处在于：比如说左延的序列其右端序列刚好跟右延的序列能够完全匹配，这样可以在很大程度上充分说明这个序列是非常的接近缺失序列的
             */
            //1:collect the sequence coincide with the condition (1)
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
        }
        // 合并左flank和右flank
        string strExtend;
        //如果存在这样的soft clip reads，那么说明，这个gap是能够被准确的全部填充的
        if(!vWholeLeftExt.empty() || !vWholeRightExt.empty()) //只要存在不为空那么就一起来搞
        {
            strExtend = CombineExtByWholeExt(strLeftExtend, strRightExtend, vWholeLeftExt, vWholeRightExt);
            itr->first->pOrgGap->enFillResult = frAllPass;
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

//新增soft clip reads
void ClsPairEndFill::AddValidSCReads(BamAlignment& al, St_GapRefinedByRept& stGapRefinedByRept,
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
                if(iEndCmp >= (int)stSCReads.strClipSeq.length()-2)
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
    }
    vSCReads.push_back(stSCReads);
}

/*--->
 * 这里只是做一个记录：
 * 如果在进行refine后，左右有一定数目的碱基的延伸，那么我们对于下一个片段也要进行相应的延伸，
 * 使得数目可以进行相应的比较，从而选取frequency最高的Nucleotide作为最后的selection
 */
string ClsPairEndFill::CombineExtByWholeExt(string strLeftExt, string strRightExt,
                                            vector<string>& vWholeLeftExt, vector<string>& vWholeRightExt)
{
    //======>Statistic the candidate for the left Extent
    //我们以right ext 为 anchor 进行定位
    //我们采取的实现方法是，将其补全，然后进行概率的统计和比较
    vector<string> vRefinedLeftExt;
    int iIndex = 0;
    int iMaxLeftLen = 0;
    if(!vWholeLeftExt.empty())
    {
        vRefinedLeftExt.resize(vWholeLeftExt.size());
        for(vector<string>::iterator itr = vWholeLeftExt.begin(); itr != vWholeLeftExt.end(); itr++)
        {
            int iExtIndex = itr->find(strRightExt);
            if(iExtIndex > 0)
                vRefinedLeftExt[iIndex] = itr->substr(0, iExtIndex);
            else
                vRefinedLeftExt[iIndex] = "";
            if(iMaxLeftLen < iExtIndex)
                iMaxLeftLen = iExtIndex;
            iIndex++;
        }
        //将里面是空(“”)的元素给删除掉
        for(vector<string>::iterator itr = vRefinedLeftExt.end()-1; itr >= vRefinedLeftExt.begin(); itr--)
        {
            if(*itr == "")
                vRefinedLeftExt.erase(itr);
        }
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
    //这里的合并的基础也是右对其
    vector<char> vRefinedWholeLeft;
    vector<int> vLeftExtAddCount;
    ClsAlgorithm::GetInstance().GetMostFreqChar(vRefinedLeftExt, vRefinedWholeLeft, vLeftExtAddCount,
                                                iMaxLeftLen, atRight); //右对齐

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
    //这里是选择的右对齐
    vector<char> vRightExtAddLeft;
    vector<int> vRightExtAddLeftCount;
    ClsAlgorithm::GetInstance().GetMostFreqChar(vRightExtAddLeftSeq, vRightExtAddLeft,
                                                vRightExtAddLeftCount, iMaxRightExtAddLeftLen, atRight); //右对齐

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
    iIndex = 0;
    int iMaxRightLen = 0;
    if(!vWholeRightExt.empty())
    {
        vRefinedRightExt.resize(vWholeRightExt.size());
        for(vector<string>::iterator itr = vWholeRightExt.begin(); itr != vWholeRightExt.end(); itr++)
        {
            int iExtIndex = itr->find(strLeftExt);
            int iRefinedLen = itr->length()-(iExtIndex + strLeftExt.length()-1);
            vRefinedRightExt[iIndex] = itr->substr(iExtIndex + strLeftExt.length(), iRefinedLen);
            if(iMaxRightLen < iRefinedLen)
                iMaxRightLen = iRefinedLen;
            iIndex++;
        }
        //将里面是空(“”)的元素给删除掉
        for(vector<string>::iterator itr = vRefinedRightExt.end()-1; itr >= vRefinedRightExt.begin(); itr--)
        {
            if(*itr == "")
                vRefinedRightExt.erase(itr);
        }
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

    //这里是选择左对齐跟之前的是有所区别的 (之前的对齐方式是右对齐)
    string strRefinedWholeRight = "";
    vector<char> vRefinedWholeRight;
    vector<int> vRefinedWholeRightCount;
    ClsAlgorithm::GetInstance().GetMostFreqChar(vRefinedRightExt, vRefinedWholeRight,
                                                vRefinedWholeRightCount, iMaxRightLen, atLeft); //这里是左对齐
    for(vector<char>::iterator itr = vRefinedWholeRight.begin(); itr != vRefinedWholeRight.end(); itr++)
    {
        if(*itr != NULL)
            strRefinedWholeRight += *itr;
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

string ClsPairEndFill::CombineExt(string strExt1, string strExt2/*, int iStart1, int iStart2, int iEnd1, int iEnd2*/)
{
    //Step1: 首先判断两边有一者为空的情况，因为这种情况是需要去作local alignment的（强行做会使得得到的start和end的值异常大，从而导致后面的计算出错）
    // 实际上如果存在一遍的数目为空，实际上这个case是不应该常发生的。因此，这个可以在一定程度上认定，这个是属于unreliable的
    if(strExt1 == "" || strExt2 == "")
        return strExt1 + strExt2;

    //我们可以考虑一下-->global alignment --> 这样说不定会有比较好的结果 -->今天到此为止吧 --->一会儿去看看python去，完成一下那个老师的作业
    /*LocalAlignment localAl;//--> 使用现有的库来进行Local Alignment
     int iStart1 = -1, iStart2 = -1, iEnd1 = -1, iEnd2 = -1;
     localAl.optAlign(strExt1, strExt2, iStart1, iEnd1, iStart2, iEnd2);
    */
    //Use我自己改过的local alignment去进行相应的比对工作
    //Notice Ext1 is the left part, Ext2 is the right part
    SmithWaterman* pAln = new SmithWaterman(strExt1, strExt2);
    pAln->align();
    string strAlnExt1;
    string strAlnExt2;
    pAln->GetAlnSeqResult(strAlnExt1, strAlnExt2);
    delete pAln;
    pAln = NULL;
    //Refine the Ext1 and Ext2:
    //基于一个思想：
    //实际上从左边走的，左边一定是最准确的
    //然后从右边走的，其右边一定是最准确的 -->
    // 在这里我们要找一个从后面开始找的函数，或者返回最后一次找到的位置！！！！，以后再思考吧，这里是肯定需要修改和替换的
    // ***不要向上面那么去干，我们直接去找第一个就可以了 -->不行，还是要先通过这个进行序列的不起，然后再使用第一次出现的作为相应的start和end
    int iStart1 = -1, iStart2 = -1, iEnd1 = -1, iEnd2 = -1;
    int iPos = -1;
    iPos = -1;
    for(int i=0; i<(int)strExt1.size() - (int)strAlnExt1.size() + 1; i++)
    {
        int iTempPos = strExt1.find(strAlnExt1.c_str(), i);
        if(iTempPos != -1)
            iPos = iTempPos;
    }
    if(iPos > 0)
    {
        int iCount = iPos;
        strExt1 = strExt1.substr(0, iCount) + strAlnExt1;
        iStart1 = strExt1.find(strAlnExt1) + 1;
        iEnd1 = iStart1 + strAlnExt1.length() - 1;
    }
    else
    {
        strExt1 = strAlnExt1;
        iStart1 = 1;
        iEnd1 = strExt1.length();
    }

    if(strExt2.find(strAlnExt2.c_str()) + strAlnExt2.length() != strExt2.length())
    {
        int iStart = strExt2.find(strAlnExt2.c_str()) + strAlnExt2.length();
        int iCount = strExt2.length() - iStart;
        strExt2 = strAlnExt2 + strExt2.substr(iStart, iCount);
        iStart2 = 1;
        iEnd2 = strAlnExt2.length();
    }
    else
    {
        strExt2 = strAlnExt2;
        iStart2 = 1;
        iEnd2 = strExt2.length();
    }
    //<--
    //下一步思想：我觉得这个地方进行相应的更改后，实际上后面也是可以进行更改的，这样后面的逻辑会变得更加的简单-->下一步再进行吧

    En_AlnState enExt1Side, enExt2Side;
    //case 0: 有一个是全匹配，那么纠结，直接最长的那个就是结果
    //For enExt1Side  // 这个的返回值是按照1为start的，并不是0 这个要注意！！！！！， 我么可以把他们搞成0；
    iStart1--;
    iEnd1--;
    iStart2--;
    iEnd2--;

    /*
    if( (iStart1 == 0 && iEnd1 == (int)strExt1.length() - 1) ||
        (iStart2 == 0 && iEnd2 == (int)strExt2.length() - 1))
    {
        return strExt1.length() >= strExt2.length() ? strExt1 : strExt2;
    }*/
    //这后面一定需要全部的重新修改，因为这里实际上成功的local alignment 并不意味着完全匹配，其中是允许indel以及mismatch的。
    //今天把这个结了～～～～～！！！

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
                    return strExt2 + strCmbExt1;
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
                    return strCmbExt1 + strExt2; //strExt1;
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
                    return strCmbExt1_1 + strExt2 + strCmbExt1_1;
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
                    return  strExt1 + strCmbExt2; //strExt2; // 也就是 strExt1
                }
                case apRight:
                {
                    strSubExt2 = strExt2.substr(iStart2, iEnd2 - iStart2 + 1);
                    strCmbExt2 = strExt2.substr(0, iStart2);
                    return strCmbExt2 + strExt1;
                }
                case apMiddle:
                {
                    strSubExt2 = strExt2.substr(iStart2, iEnd2 - iStart2 + 1);
                    strCmbExt2_1 = strExt2.substr(0, iStart2);
                    strCmbExt2_2 = strExt2.substr(iEnd2+1, strExt2.length() - iEnd2);
                    return strCmbExt2_1 + strExt1 + strCmbExt2_2;
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




