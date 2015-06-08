#include "clsrepeatbuild.h"
#include "algorithm" //This is the system file
#include <math.h>
#include <ctype.h>

ClsRepeatBuild::ClsRepeatBuild(): m_strBamPath("")
{
}

ClsRepeatBuild::~ClsRepeatBuild()
{}

void ClsRepeatBuild::Init(string strScafPath, string strReads1Path, string strReads2Path)
{
    //Use Mapping all to get all of reasult
    m_strScafPath = strScafPath;
    m_strReads1Path = strReads1Path;
    m_strReads2Path = strReads2Path;
}

void ClsRepeatBuild::GapFillByAbnPE(St_ScaffoldUnit& stScafUnit, vector<St_Fasta>& vScafContigSet)
{
    m_strBamPath = ClsAlgorithm::GetInstance().CreateBamFileForMultiRefPEReads( m_strScafPath,
                                                                                m_strReads1Path,
                                                                                m_strReads2Path,
                                                                                "DraftGeno");

    for(vector<St_Gap>::iterator itr = stScafUnit.vGap.begin(); itr != stScafUnit.vGap.end(); itr++)
    {
        FillCurGapByPEReads(stScafUnit, *itr, vScafContigSet);
    }
}

/*
 * Idea：
 * 1：定义insert length: 2 * Reads Length + Interval Length
 * 2: 找到起点
 * 3：找到终点
 * *起点跟终点的依据在于，使得map上的reads的mate落在gap之间 (这就意味这这些reads的mat本来就不应该被map上)
 * 4：有效范围内是否存在异常的PE reads
 * 5：对异常PE的定义如下
 * (1) 如果严重大于insert size 两倍 (同一scaffold上)
 * (2) 位于另外一个不同的draft Geno上
 */

enum En_AbnType{abtSameRef, abtDiffRef, abtMax}; //abnormal type

struct St_AbnMapReads // abnormal mapping reads
{
    string strMateSeq;
    int iRefID;
    En_AbnType enType;
    int iStartPos;

    St_AbnMapReads():strMateSeq(""), iRefID(-1), enType(abtMax), iStartPos(-1)
    {}
};

/*
 * the case we need to take into consideration
 * (1)The start should be the Start(Gap) - Len(MaxInsertion)
 * (2)The End part should be the the Start(Gap) - 1
 * (3)We need to check the statement of its mates
 *   a)If it should be the same ref as current one
 * Caution:
 * 1: we should know that: the scaffold itself many be not correct enough!!!! If it is, it will be hard for us to get the real gap filling result
 * #这里还有可能scaffold很长，所以就multi map了, 在这个case里面，map到的结果都是在同一个draft geno上面
 *
 */
void ClsRepeatBuild::FillCurGapByPEReads(St_ScaffoldUnit& stScafUnit, St_Gap& stGap, vector<St_Fasta>& vScafContigSet)
{
    //Data Prepare
    int iReadsLen = 100;
    int iIntervalLen = 200;
    int iInsertLen = iIntervalLen + iReadsLen;
    int iAbnormalLen = iInsertLen * 1.5;

    vector<St_AbnMapReads> vAbnMapReads;
    //Reads Mapping
    BamReader* pBamReader = new BamReader();
    pBamReader->Open(m_strBamPath);
    pBamReader->OpenIndex(m_strBamPath + ".bai");
    BamAlignment al;
    //From left to right
    while(pBamReader->GetNextAlignment(al))  //现在最简单的方式：找到一个就可以了--->
    {        
        //只有match到当前scaffold上的才能接受
        if(!al.IsMapped())
            continue;
        if(al.RefID != stScafUnit.iID) //这个证明是配到了当前的reference
            continue;
        //if(al.Position < iInsertLen)
        //    continue;
        //到了一定的距离我们还应该停止，因为没有必要了 (对于当前的gap而言没有必要了)
        if(al.Position > stGap.iStartPos) //do it again and again until it reaches the start position of gap
            break;

        //Erase the reads which is too far away from the current gap
        if(al.Position + iReadsLen + iIntervalLen + iReadsLen < stGap.iStartPos)
            continue;

        //好了现在这个时候就是正经能够对上的时候了
        //--->
        //如果发生了这个情况，那么意味着
        //1：gap本身要小于insertion length
        //2：scaffol足够的长，能够使得reads的两个pair在其上面都能够匹配上
        if(al.IsMateMapped()) // Mat Map 上了，
        {
            //We try to check if they should be in the same draft Geno-->
            bool bShouldSameRef = false;
            int iEptMateStart = al.Position + iReadsLen + iIntervalLen + stGap.iLen + .5 * iReadsLen;
            if(iEptMateStart <= (int)stScafUnit.strSeq.length())
                bShouldSameRef = true;
            //<--
            if(bShouldSameRef)
            {
                if(al.RefID != al.MateRefID) //this case means repetitive fragment happend(missing current ref, appeared in other geno)
                {
                    //Map to another place
                    St_AbnMapReads stAbnMapReads;
                    stAbnMapReads.enType = abtDiffRef;
                    stAbnMapReads.iRefID = al.MateRefID;
                    stAbnMapReads.iStartPos = al.MatePosition;
                    int iLen = iReadsLen;
                    if(iLen + al.MatePosition > (int)vScafContigSet[stAbnMapReads.iRefID].strSeq.length())
                        iLen = vScafContigSet[stAbnMapReads.iRefID].strSeq.length() - al.MatePosition;
                    stAbnMapReads.strMateSeq = vScafContigSet[stAbnMapReads.iRefID].strSeq.substr(al.MatePosition, iLen);
                    //valid type 1.1:scaffold 很长，在本scaffold内部就存在multi mapping
                    vAbnMapReads.push_back(stAbnMapReads);
                }
                else
                {
                    if(abs(abs(al.InsertSize) - stGap.iLen) < iAbnormalLen) // 如果这个长度是正常的，我们认为可以不进行考虑
                        continue;
                    else
                    {
                        //如果长度不正常，我们就要进行相应的值存储:
                        St_AbnMapReads stAbnMapReads;
                        stAbnMapReads.enType = abtSameRef;
                        stAbnMapReads.iRefID = al.MateRefID;
                        stAbnMapReads.iStartPos = al.MatePosition;
                        int iLen = iReadsLen;
                        if(iLen + al.MatePosition > (int)vScafContigSet[stAbnMapReads.iRefID].strSeq.length())
                            iLen = vScafContigSet[stAbnMapReads.iRefID].strSeq.length() - al.MatePosition;
                        stAbnMapReads.strMateSeq = vScafContigSet[stAbnMapReads.iRefID].strSeq.substr(al.MatePosition, iLen);
                        //valid type 1.1:scaffold 很长，在本scaffold内部就存在multi mapping
                        vAbnMapReads.push_back(stAbnMapReads);
                    }
                }
                continue;
            }
            else // should Map to the different geno
            {
                //This case may also include the repetitive case, but we can not extect it out of from such case
                continue;
            }
        }
        else // this means the left part cannot be mapped
        {
            //没配上----> 于是就可以完全不考虑了
            continue;
        }
        //<---
    }
    delete pBamReader;
    pBamReader = NULL;

    //-->
    //因为我们针对的是已经sorted后的片段，因此我们不需要进行右边的比对了，左边若比不上，右边肯定是不需要考虑的
    //<--
    //Handle with the vAbnMapReads ----------->
    //现在需要找寻的是我得到了一组的reads，看看我们应该怎样去将他们combine到一起来
    for(vector<St_AbnMapReads>::iterator itr = vAbnMapReads.begin(); itr != vAbnMapReads.end(); itr++)
    {
        // 到目前为止，虽然我们得到了一些异常的（全部是配到另外一个draft geno上面的）的reads，但是根据我的手工查看，我觉得这些reads都不值得用来作片段的还原。
        // 或者你可以说这样的到的所谓的异常发的reads，其来源以及本身的值是没法估计的，他们基本上跟repetitive 片段没有任何关系
        // 但是为啥会出现这样的match不上的情况，这个还需要进一步深究。
    }
    //<----
}

/* -------->Go!!!!!!!!!!!!!!!!!!!
 * This is the most important part of my project
 * Strategy:
 * (1)Use Kmer Counting get the high frequency Kmer  (Use Jellyfish for Kmer counting)
 * (2)Extend those kmer and try to find out
 *   (a)the the maximum length of repeats which contains the current kmer (Use Blast)
 *   (b)the position of "N"
 * (3)Try to combine the relevant fragment of the same frequency kmer (Use Blast)
 * (4)Use the new created repeat to fill those Gap
 * Comments:
 * we do not need the info of pair-end reads currently
 */
void ClsRepeatBuild::BuildPseudoRepeat(St_ScaffoldUnit& stScaf, vector<St_Fasta>& vDraftSet,
                                       vector<St_ScaffoldUnit>& vScafSet)
{
    //Step1: Use build result of k-mer frequency by Jelly fish
    string strJF = KmerCounting(m_strScafPath);
    //Step2: Filter those result of k-mer freqency with the flank of gap by using current scaffold    
    for(vector<St_Gap>::iterator itr = stScaf.vGap.begin(); itr != stScaf.vGap.end(); itr++)
    {        
        FillCurGapByDraftGeno(*itr, strJF, vDraftSet, vScafSet);
    }

    //Step3: Try to extend those filtered result and get the longer sequence of repeats
    //Step4: Use those kind of repeats to fill those gaps
}

/*Use Jelly fish to make Kmer Counting
//Our strategy:
//1: 使用jellyfish组建kmer counting
//2: 收集左右的flank
//3: 使用16mer进行相应的查询
//4： 将allele frequency的查询结果存储起来
//5： 根据一组连续的kmercouting 的值，相应的flank的freqency结果
//6： record: 3
//7:  record:
*/
string ClsRepeatBuild::KmerCounting(string strDraftGeno)
{
    //Generate jf file (the output file of jellyfish)
    string strJellyfish = ClsAlgorithm::GetInstance().GetHigherFolderPath(ClsAlgorithm::GetInstance().GetCurExeFolderPath()) +
                          "ThirdPartyTools/Jellyfish/bin/jellyfish";
    string strOutputFile = ClsAlgorithm::GetInstance().GetHigherFolderPath(ClsAlgorithm::GetInstance().GetCurExeFolderPath()) +
                          "ThirdPartyTools/Jellyfish/data/mer_counts.jf";
    string strCmd = "";
    int iKmerLength = 16;
    //Option: "-C" let both regular direction and reverse complementary be caculated together
    strCmd = strJellyfish + " count -m " + IntToStr(iKmerLength) + " -s 100M -t 10 -C " + strDraftGeno +
            " -o " + strOutputFile;
    system(strCmd.c_str());
    return strOutputFile;
}

void ClsRepeatBuild::FillCurGapByDraftGeno(St_Gap& stGap, string& strJF, vector<St_Fasta>& vDraftSet,
                                           vector<St_ScaffoldUnit>& vScafSet)
{
    //Get the left flank and the right flank to check which kind of them has the high kmer counting result
    //For left flank
    string strLeftAnchor = "";
    int iLeftFreq = 0;
    GetAnchorInfo(stGap.strLeftFlank, strJF, strLeftAnchor, iLeftFreq, asLeft);
    //For Right Flank
    string strRightAnchor = "";
    int iRightFreq = 0;
    GetAnchorInfo(stGap.strRightFlank, strJF, strRightAnchor, iRightFreq, asRight);
    //Chick which side should be took into consideration --->Go!!!!
    En_AnchorSide enAnchorSide = asMax;
    int iMinFreq = 7;
    if(iLeftFreq < iMinFreq && iRightFreq < iMinFreq)
        return;
    else
    {
        if(iRightFreq > iLeftFreq)
        {
            if(iRightFreq < iMinFreq)
                return;
            else
                enAnchorSide = asRight;
        }
        else
        {
            if(iLeftFreq < iMinFreq)
                return;
            enAnchorSide = asLeft;
        }
    }

    //try to extend the anchor and fill it relevant gap
    string strFillResult = "";    
    //1: Let's try to find out if there is the gap which share the same flank with current one and has been filled by the algorithm
    //-->Go!!!!!Scaffold 79 support this assumption
    //2: Extend it
    if(enAnchorSide == asLeft) //left side
    {
        bool bFillByCopy = FillByCopyFrag(strLeftAnchor, vDraftSet, vScafSet,
                                          strFillResult, enAnchorSide);
        if(!bFillByCopy)
        {
            //consider both reversary and its reverse complementory
            AnchorLeftExtend(strLeftAnchor, iLeftFreq,
                            vDraftSet, strJF, strFillResult); //---> Do it after launch--->!!!!
        }
    }
    else if(enAnchorSide == asRight) //How to do the right side
    {
        bool bFillByCopy = FillByCopyFrag(strRightAnchor, vDraftSet, vScafSet,
                                          strFillResult, enAnchorSide);
        if(!bFillByCopy)
        {
            //consider both reversary and its reverse complementory
            AnchorRightExtend(strRightAnchor, iLeftFreq,
                              vDraftSet, strJF, strFillResult); //---> Do it after launch--->!!!!
        }
    }
    else{}
    //Set result
    if(strFillResult != "") // It means some kind of sequence could be found to fill such gap (this subsequece must relevant with repeats)
    {
        //unsigned int i = strFillResult.size();
        stGap.strFillSeq = strFillResult;
        stGap.bFilled = true;
        //在这个case中，我们坚持认为：最后一个most frequency的kmer 不会跟另一边的flank产生交集，
        //因此这也意味着我们不需要在这个阶段对enFillResult进行赋值
        if(enAnchorSide == asLeft)
            stGap.enGapType = gtRight;
        else if(enAnchorSide == asRight)
            stGap.enGapType = gtLeft;
        else{}
    }
}

/* Step 1:Mapping such part back to drft geno
 * Step 2:Make some filter
 * Step 3: Combine those connected parts
 * Notice:
 * 1: The condition of extension
 * 2: The condition of termination
 * Use Bwa to map such flank back to the draft geno --> Map back to both scaffold and contig directly
   1:Stop if meet N
   2:Stop in the end of each sequence
   Notice: the extension orientation for both normal one and reverse comolementary
*/

void ClsRepeatBuild::AnchorExtdByBWA(string strAnchor, int iFreq, vector<St_Fasta>& vDraftSet)
{
    //Step 1: Use BWA for mapping
    //(1) Build Reads
    string strRootPath = ClsAlgorithm::GetInstance().GetHigherFolderPath(ClsAlgorithm::GetInstance().GetCurExeFolderPath());
    string strReadsPathFa = strRootPath + "TempFile/RepeatsExtd.fa";
    string strReadsPathFq = strRootPath + "TempFile/RepeatsExtd.fq";
    ofstream ofs;
    ofs.open(strReadsPathFa.c_str());
    ofs << ">FlankAnchor"<< endl << strAnchor << endl;
    ofs.close();
    //(2) Generate fastq file
    string strCmd = "perl " + strRootPath + "/ThirdPartyTools/fasta_to_fastq.pl " +
             strReadsPathFa + " > " + strReadsPathFq;
    system(strCmd.c_str());
    //(3) Make BWA
    string strBamPath = ClsAlgorithm::GetInstance().CreateBamFileForMultiRefSingleReads(m_strScafPath,
                                                                       strReadsPathFq, "RepeatAnchor");
    //Step 2: Parse such bam file and try to get all of their remaining sequence
    //Notice: How to chose the remaining sequence by either original seq and reverse complementary
    BamReader* pBamReader = new BamReader();
    pBamReader->Open(strBamPath);
    pBamReader->OpenIndex(strBamPath + ".bai");
    BamAlignment al;
    vector<string> vExtSeq;
    while(pBamReader->GetNextAlignment(al))  //现在最简单的方式：找到一个就可以了--->
    {
        //All the mapping reads are in the first part, since it has been sorted
        if(!al.IsMapped())
            break;
        string strSeq = "";
        if(al.IsReverseStrand()) //Reverse complementary
        {
            strSeq = ClsAlgorithm::GetInstance().GetReverseCompelement(vDraftSet[al.RefID].strSeq.substr(0, al.Position));
        }
        else //Normal Orientation
        {
            strSeq = vDraftSet[al.RefID].strSeq.substr(al.GetEndPosition() + 1,
                                                       vDraftSet[al.RefID].strSeq.length() - al.GetEndPosition() - 1);
        }
        vExtSeq.push_back(strSeq);
    }
    delete pBamReader;
    pBamReader = NULL;
    //Step 3:Try to merge those subsequence
}

bool ClsRepeatBuild::FillByCopyFrag(string strAnchor, vector<St_Fasta>& vDraftSet,
                                    vector<St_ScaffoldUnit>& vScafSet, string& strCmb,
                                    En_AnchorSide enAnchorSide)
{
    //Algorithm:
    //1: Find out all of places and then
    //2: Check if such gap
    int iScafIndex = 0;
    bool bGot = false;
    for(vector<St_Fasta>::iterator itr = vDraftSet.begin(); itr != vDraftSet.end(); itr++)
    {
        if(iScafIndex > (int)vScafSet.size())
            break;
        int iQueryPos = 0;
        while(itr->strSeq.find(strAnchor, iQueryPos) != string::npos)
        {
            int iStart = itr->strSeq.find(strAnchor, iQueryPos);
            switch(enAnchorSide)
            {
                case asLeft:
                {
                    if(iStart + strAnchor.length() + 1 >= itr->strSeq.length())
                        break;
                    if(itr->strSeq.at(iStart + strAnchor.length() + 1) == 'N') // It means N behind it
                    {
                        //Got one Gap
                        for(vector<St_Gap>::iterator itrGap = vScafSet[iScafIndex].vGap.begin();
                            itrGap != vScafSet[iScafIndex].vGap.end(); itrGap++)
                        {
                            if(itrGap->iStartPos = iStart + strAnchor.length() + 1) //Get current gap
                            {
                                if(itrGap->strFillSeq != "")
                                {
                                    strCmb = itrGap->strFillSeq;
                                    bGot = true;
                                    break;
                                }
                            }
                        }
                        if(bGot)
                            break;
                    }
                    break;
                }
                case asRight:
                {
                    if(itr->strSeq.at(iStart - 1) == 'N') // It means N behind it
                    {
                        //Got one Gap
                        for(vector<St_Gap>::iterator itrGap = vScafSet[iScafIndex].vGap.begin();
                            itrGap != vScafSet[iScafIndex].vGap.end(); itrGap++)
                        {
                            if(itrGap->iStartPos = iStart - 1) //Get current gap
                            {
                                if(itrGap->strFillSeq != "")
                                {
                                    strCmb = itrGap->strFillSeq;
                                    bGot = true;
                                    break;
                                }
                            }
                        }
                        if(bGot)
                            break;
                    }
                    break;
                }
                default:
                    break;
            }
            if(bGot)
                break;
            iQueryPos = iStart + strAnchor.length() + 1;
        }
        if(bGot)
            break;
        iScafIndex++;
    }
    return bGot;
}

void ClsRepeatBuild::AnchorLeftExtend(string strAnchor, int iFreq,
                                      vector<St_Fasta>& vDraftSet, string& strJF, string& strCmb)
{
    //We use string finding directly to check the possible sequence
    vector<string> vExtSeq;
    for(vector<St_Fasta>::iterator itr = vDraftSet.begin(); itr != vDraftSet.end(); itr++)
    {
        //Normal Orientation: use the sequence after the anchor
        //Find the sequece after anchor
        int iQueryPos = 0;
        bool bHit = false;
        while(iQueryPos < (int)itr->strSeq.length())
        {
            if(::FindStrPosByNS(itr->strSeq, strAnchor, iQueryPos) != string::npos) // Find it
            {
                int iStartPos = ::FindStrPosByNS(itr->strSeq, strAnchor, iQueryPos) +
                                (int)strAnchor.length();
                int iEndPos = itr->strSeq.length()-1;
                if(itr->strSeq.find("N", iStartPos) != string::npos)
                {
                    iEndPos = itr->strSeq.find("N", iStartPos) - 1;
                    if(::FindStrPosByNS(itr->strSeq, strAnchor, iStartPos) != string::npos)
                    {
                        int iNextAnchorPos = ::FindStrPosByNS(itr->strSeq, strAnchor, iStartPos);
                        if(iEndPos > iNextAnchorPos)
                            iEndPos = iNextAnchorPos - 1;
                    }
                }
                //Save the result --->
                if(iEndPos - iStartPos + 1 > 0)
                    vExtSeq.push_back(itr->strSeq.substr(iStartPos, iEndPos - iStartPos + 1));
                bHit = true;
                //<---
                iQueryPos = iEndPos + 1;
            }
            else
                iQueryPos = itr->strSeq.length();
        }
        if(bHit) //It means the draft geno belongs to the ragular sequence
        {
            continue;
        }
        //If the original orientation do not match current draft geno, we try to search its reverse complementary
        //One thing should be noticed: Do not need to find the next anchor, since we use the part from start to current!!!!!
        //Reverse Complementary Orientation:
        //Notice: use the sequence before the anchor and then make the reverse complementary
        //Find the sequence before the anchor
        iQueryPos = 0;
        string strRC = ClsAlgorithm::GetInstance().GetReverseCompelement(strAnchor);
        while(iQueryPos < (int)itr->strSeq.length())
        {
            if(::FindStrPosByNS(itr->strSeq, strRC, iQueryPos) != string::npos) // Find it
            {
                int iStartPos = iQueryPos;
                int iEndPos = ::FindStrPosByNS(itr->strSeq, strRC, iQueryPos) - 1;
                if(itr->strSeq.find("N", iStartPos) != string::npos)
                {
                    int iLastNPos = itr->strSeq.find("N", iStartPos);
                    if(iLastNPos > iEndPos) //we just need to consider the N between start and end
                    {}
                    else
                    {
                        while(itr->strSeq.substr(iLastNPos, 1) == "N")
                            iLastNPos++;
                        iStartPos = iLastNPos;
                    }
                }
                //Save the result --->
                if(iEndPos - iStartPos + 1 > 0)
                    vExtSeq.push_back(ClsAlgorithm::GetInstance().GetReverseCompelement(itr->strSeq.substr(iStartPos, iEndPos - iStartPos + 1)));
                //<---
                iQueryPos = iEndPos + strRC.length() + 1;
            }
            else
                iQueryPos = itr->strSeq.length();
        }
    }
    //Combine vExtSeq --->
    //Check if all of the sequence contained by vExtSeq is correct one    
    int iPrevPos = (int)strCmb.length();
    GetCombinedSeq(strCmb, vExtSeq, strJF, iFreq, atLeft);
    //----------------->
    //Try to use the last 16 Kmer to located a new group of sequence and them combined them together.
    //Check Kmer Frequency
    //Find the condition of termination
    //1: Continue, if all of sub k-mer with the high frequency
    //2: Stop if there is the low frequency kmer
    //3: Get the breakpoint between the high freqency one and the low frequency one
    bool bTerminate = true;
    int iKmerLen = 16;
    int iCurFreq = 0;
    string strSubSeq = "";
    if((int)strCmb.length() >= iKmerLen)
    {
        //Here: we need some way to use the most useful kmer as the achor rather than just the last one --->Go!!!
        strSubSeq = strCmb.substr(strCmb.length() - iKmerLen, iKmerLen);
        iCurFreq = GetKmerFreq(strSubSeq, strJF);
        if( iCurFreq > iFreq ||
            abs(iCurFreq - iFreq) < (iFreq * .3))
            bTerminate = false;
    }
    //<-------
    // It is either the low frequency one or not long enough-->Try to find out the last high freqency kmer
    if(bTerminate)
    {
        if((int)strCmb.length() >= iKmerLen)
        {            
            if(iPrevPos == 0) //This is the first time combined sequence
            {               
                //In this case we need to find it from the end to start rather than from start to end,
                //since some mismatch will make the result be rejected
                bool bHighFreq = false; //Target find the first high freqency one from end to start !!!!!
                int iOffSet = 1;
                while(!bHighFreq &&  //Jump out when low frequency hit
                      iOffSet < (int)strCmb.length() - iKmerLen) //Jump out when the end met
                {
                    strSubSeq = strCmb.substr(strCmb.length() - iKmerLen - iOffSet, iKmerLen);
                    //Skip "AAAAAAAAAAAAAAAA"
                    if(strSubSeq.find("AAAAAAAAAAAA") != string::npos) //Do not allowed a lot of A be existed in query seed
                    {
                        iOffSet++;
                        continue;
                    }
                    iCurFreq = GetKmerFreq(strSubSeq, strJF);                    
                    if( iCurFreq > iFreq ||
                        abs(iCurFreq - iFreq) < (iFreq * .3)) // Belongs to the high freqency one
                        bHighFreq = true;
                    else
                        iOffSet++;
                }
                if(bHighFreq)
                    strCmb = strCmb.substr(0, strCmb.length() - iOffSet);
                else
                    strCmb = "";
            }
            else //Second, third, forth or .... time
            {
                bool bHighFreq = true;
                int iOffSet = 1;
                string strSubSeq = "";
                while(bHighFreq &&  //Jump out when low frequency hit
                      iOffSet + iPrevPos < (int)strCmb.length()) //Jump out when the end met
                {
                    strSubSeq = strCmb.substr(iPrevPos - iKmerLen + iOffSet, iKmerLen);
                    iCurFreq = GetKmerFreq(strSubSeq, strJF);
                    if( iCurFreq > iFreq ||
                        abs(iCurFreq - iFreq) < (iFreq * .3))
                        iOffSet++;
                    else
                        bHighFreq = false;
                }
                //In this place, we extend it from start to end. As such, the low frequency could be viewed as the terminate condition
                if(!bHighFreq)
                    strCmb = strCmb.substr(0, iPrevPos + iOffSet - 1);                                
            }
        }
        return;
    }
    else
        AnchorLeftExtend(strSubSeq, iFreq, vDraftSet, strJF, strCmb);
}

void ClsRepeatBuild::AnchorRightExtend(string strAnchor, int iFreq,
                       vector<St_Fasta>& vDraftSet, string& strJF, string& strCmb)
{
    //We use string finding directly to check the possible sequence
    vector<string> vExtSeq;
    for(vector<St_Fasta>::iterator itr = vDraftSet.begin(); itr != vDraftSet.end(); itr++)
    {
        //Normal Orientation: use the sequence after the anchor
        //Find the sequence before anchor
        int iQueryPos = 0;
        bool bHit = false;
        while(iQueryPos < (int)itr->strSeq.length())
        {
            if(::FindStrPosByNS(itr->strSeq, strAnchor, iQueryPos) != string::npos) // Find it
            {
                int iStartPos = iQueryPos;
                int iEndPos = itr->strSeq.find(strAnchor, iQueryPos) - 1;
                if(itr->strSeq.find("N", iStartPos) != string::npos)
                {
                    int iLastPos = itr->strSeq.find("N", iStartPos);
                    if(iLastPos > iEndPos) //we just need to consider the N between start and end
                    {}
                    else
                    {
                        while(itr->strSeq.substr(iLastPos, 1) == "N")
                            iLastPos++;
                        iStartPos = iLastPos;
                    }
                }
                //Save the result --->
                if(iEndPos - iStartPos + 1 > 0)
                    vExtSeq.push_back(itr->strSeq.substr(iStartPos, iEndPos - iStartPos + 1));
                bHit = true;
                //<---
                iQueryPos = iEndPos + (int)strAnchor.length() + 1;
            }
            else
                iQueryPos = itr->strSeq.length();
        }
        if(bHit) //It means the draft geno belongs to the ragular sequence
        {
            continue;
        }
        //Reverse Complementary Orientation:
        //Notice: use the sequence before the anchor and then make the reverse complementary
        //Find the sequence after anchor
        iQueryPos = 0;
        string strRC = ClsAlgorithm::GetInstance().GetReverseCompelement(strAnchor);
        while(iQueryPos < (int)itr->strSeq.length())
        {
            if(::FindStrPosByNS(itr->strSeq, strRC, iQueryPos) != string::npos) // Find it
            {
                int iStartPos = ::FindStrPosByNS(itr->strSeq, strRC, iQueryPos) + (int)strRC.length();
                int iEndPos = (int)itr->strSeq.length() - 1;
                if(itr->strSeq.find("N", iStartPos) != string::npos)
                {
                    iEndPos = itr->strSeq.find("N", iStartPos) - 1;
                    if(::FindStrPosByNS(itr->strSeq, strRC, iStartPos) != string::npos)
                    {
                        int iNextAnchorPos = ::FindStrPosByNS(itr->strSeq, strRC, iStartPos);
                        if(iEndPos > iNextAnchorPos)
                            iEndPos = iNextAnchorPos - 1;
                    }
                }
                //Save the result --->
                if(iEndPos - iStartPos + 1 > 0)
                    vExtSeq.push_back(ClsAlgorithm::GetInstance().GetReverseCompelement(itr->strSeq.substr(iStartPos, iEndPos - iStartPos + 1)));
                //<---
                iQueryPos = iEndPos + 1;
            }
            else
                iQueryPos = itr->strSeq.length();
        }
    }
    //Combine vExtSeq --->
    //Check if all of the sequence contained by vExtSeq is correct one
    int iPrevPos = (int)strCmb.length();
    GetCombinedSeq(strCmb, vExtSeq, strJF, iFreq, atRight);
    //----------------->
    //Try to use the first 16 Kmer to located a new group of sequence and them combined them together.
    //Check Kmer Frequency
    //Find the condition of termination
    //1: Continue, if all of sub k-mer with the high frequency
    //2: Stop if there is the low frequency kmer
    //3: Get the breakpoint between the high freqency one and the low frequency one
    bool bTerminate = true;
    int iKmerLen = 16;
    int iCurFreq = 0;
    string strSubSeq = "";
    if((int)strCmb.length() >= iKmerLen)
    {
        strSubSeq = strCmb.substr(0, iKmerLen);
        iCurFreq = GetKmerFreq(strSubSeq, strJF);
        if( iCurFreq > iFreq ||
            abs(iCurFreq - iFreq) < (iFreq * .3))
            bTerminate = false;
    }
    //<-------
    // It is either the low frequency one or not long enough-->Try to find out the last high freqency kmer
    if(bTerminate)
    {
        if((int)strCmb.length() >= iKmerLen) //In this case: the anchor is located in the right side of flank
        {            
            if(iPrevPos == 0) //The first time
            {
                //In this case we need to find it from the end to start rather than from start to end,
                //since some mismatch will make the result be rejected
                bool bHighFreq = false; //Target find the first high freqency one from end to start !!!!!
                int iOffSet = 0;
                while(bHighFreq &&  //Jump out when low frequency hit
                      iOffSet < (int)strCmb.length() - iKmerLen) //Jump out when the begin met
                {
                    strSubSeq = strCmb.substr(iOffSet, iKmerLen); // we need backward here!!!!
                    iCurFreq = GetKmerFreq(strSubSeq, strJF);
                    if( iCurFreq > iFreq ||
                        abs(iCurFreq - iFreq) < (iFreq * .3)) // Belongs to the high freqency one
                        bHighFreq = true;
                    else
                        iOffSet++;
                }
                if(bHighFreq)
                    strCmb = strCmb.substr(iOffSet, strCmb.length() - iOffSet);
                else
                    strCmb = "";
            }
            else // More than one time
            {
                bool bHighFreq = true;
                int iOffSet = 0;
                string strSubSeq = "";
                while(bHighFreq &&  //Jump out when low frequency hit
                      (int)strCmb.length() - iPrevPos - iOffSet >= 0) //Jump out when the begin met
                {
                    strSubSeq = strCmb.substr(iPrevPos - iOffSet, iKmerLen); // we need backward here!!!!
                    iCurFreq = GetKmerFreq(strSubSeq, strJF);
                    if( iCurFreq > iFreq ||
                        abs(iCurFreq - iFreq) < (iFreq * .3))
                        iOffSet++;
                    else
                        bHighFreq = false;
                }
                if(!bHighFreq)
                    strCmb = strCmb.substr(iPrevPos - iOffSet + 1, (int)strCmb.length() - (iPrevPos - iOffSet + 1));
                else
                    strCmb = "";
            }
        }
        return;
    }
    else
        AnchorRightExtend(strSubSeq, iFreq, vDraftSet, strJF, strCmb);
}

void ClsRepeatBuild::GetCombinedSeq(string& strCmbSeq, vector<string>& vSeqSet, string strJF,
                                    int iFreq, En_AlignType enAlignType)
{
    if(vSeqSet.empty()) //do not conside the case without candidate
        return;
    int iKmerLen =16;
    //strCmbSeq = ""; //it should keep its org value,since there is traverse for its parent function

    //1: First: Erase the too short one (Less than kmer)
    for(vector<string>::iterator itr = vSeqSet.end() - 1; itr >= vSeqSet.begin(); itr--)
    {
        if(itr->length() < iKmerLen * 0.5)
            vSeqSet.erase(itr);
    }  
    if(vSeqSet.empty())
        return;

    //2: Second: 试图找到以lower case 开头，的那个字串，从而求得相应子串较大值，将其作为一个截断的标准：
    //Motivation： 较长的由纯小写字母组成的字串实际上在很大程度上较精确的对应了某一个repeat
    //We just care about the sub sequence which start from the lower case!!!
    //Do not need to do any tolerance of string length, since some abnormal one have been erased by step 1
    int iContinuousLower = 0;
    for(vector<string>::iterator itr = vSeqSet.begin(); itr != vSeqSet.end(); itr++)
    {
        int iIndex = 0;
        while(islower(itr->at(iIndex)))
            iIndex++;
        if(iContinuousLower < iIndex)
            iContinuousLower = iIndex;
    }

    //3: Third: Make the cutoff by allele frequency
    vector<int> vLowFreqLen;
    string strKmer = "";
    for(vector<string>::iterator itr = vSeqSet.begin(); itr != vSeqSet.end(); itr++)
    {
        if((int)itr->length() < iKmerLen)
            continue;
        strKmer = itr->substr(itr->length() - iKmerLen, iKmerLen);
        int iCurFreq = GetKmerFreq(*itr, strJF);
        //record the sequence length with the low kmer frequency value
        //The tail with low allele frequecy could be viwed as the max length of candidate sequence
        if( iCurFreq < iFreq ||
            (iFreq - iCurFreq) > (iFreq * .5)) // do not use abs,since iCurSeq the large the better.
        vLowFreqLen.push_back(itr->length());
    }
    int iMaxLen = -1;
    if(!vLowFreqLen.empty()) // If it is not empty
    {
        sort(vLowFreqLen.begin(), vLowFreqLen.end());
        // We extend the accaptable length from 1X to 2X.--> I think it is reasonable
        if(iContinuousLower > iFreq * 2 && vLowFreqLen[0] > iContinuousLower * 3) // to avoid the length of lowFreqLen is too large
            iMaxLen =  iContinuousLower * 2;
        else
            iMaxLen = (vLowFreqLen[0] > iContinuousLower ? vLowFreqLen[0] : iContinuousLower) * 2;
    }
    else
    {
        if(iContinuousLower > iKmerLen * 2)
            iMaxLen = iContinuousLower;
    }
    if(iMaxLen > 0)
    {
        for(vector<string>::iterator itr = vSeqSet.begin(); itr != vSeqSet.end(); itr++)
        {
            if((int)itr->length() > iMaxLen)
                *itr = itr->substr(0, iMaxLen);
        }
    }
    //4: Make some erase of the final sequence
    //4.1: If there are several sequence which share the similar length but some of them are totally upper case.
    //Then, we need to keep the sequence with whole upper case, and erase others
    sort(vSeqSet.begin(), vSeqSet.end(), sortfunctionLength);
    if(vSeqSet.size() > 1) // use the second last one as reference --> since the last one will bring some abnormal result (abnormal large)
        iMaxLen = (vSeqSet.end() - 2)->length();
    else
        iMaxLen = (vSeqSet.end() - 1)->length(); // use the last one
    bool bWholeUpperCaseFilter = false;
    for(vector<string>::iterator itr = vSeqSet.end() - 1; itr >= vSeqSet.begin(); itr--)
    {
        if(itr->length() > (iMaxLen * .8) &&
           ClsAlgorithm::GetInstance().IsWholeUpperCase(*itr))
        {
            bWholeUpperCaseFilter = true;
            break;
        }
    }

    if(bWholeUpperCaseFilter)
    {
        for(vector<string>::iterator itr = vSeqSet.end() - 1; itr >= vSeqSet.begin(); itr--)
        {
            if(!ClsAlgorithm::GetInstance().IsWholeUpperCase(*itr))
                vSeqSet.erase(itr);
        }
    }
    else
    {
        //4.2: Forth: of all erase some sequence with the lower case exsited
        //In this place, we do not need to erase the sequence contain the lower case, since the sequence
        //which only contain the upper case may be too short to give you all of information.
        //So our strategy:
        //I need to add the upper case and the lower case adjustment in the function of "GetMostFreqChar!!!"
        //Notice: we could just try to change the lower case from acgt to " " to coincide with the later combination.
        bool bHasUpper = false;
        bool bHasLower = false;
        vector<char> vCharSet;

        for(int i=0; i<iMaxLen; i++)
        {
            bHasUpper = false;
            vCharSet.clear();
            for(vector<string>::iterator itr = vSeqSet.begin(); itr != vSeqSet.end(); itr++)
            {
                if(i < (int)itr->length())
                {
                    vCharSet.push_back(itr->at(i));
                    if(isupper(itr->at(i)))
                        bHasUpper = true;
                    else if(islower(itr->at(i)))
                        bHasLower = true;
                }
                else
                    // 我们在这里insert n的目的是：为了滥竽充数，使得后面统计index比较方便，没有其他特殊的意思，N在此仅仅是作为一种特殊字符
                    vCharSet.push_back('N');
            }
            //Erase the sequence relevant with lower case if there is uppercase.
            if(bHasUpper && bHasLower)
            {
                for(int j = (int)vCharSet.size() - 1; j >= 0; j--)
                {
                    if(vCharSet[j] == 'N')
                        continue;
                    if(islower(vCharSet[j]))
                    {
                        //vSeqSet.erase(vSeqSet.begin() + j); //Here we do not erase you
                        vSeqSet[j].at(i) = ' ';
                    }
                }
                //iMaxLen = (vSeqSet.end() - 1)->length();
            }
            //if((int)vSeqSet.size() <= 1)
            //    break;
        }
    }

    string strSeq = "";
    if(vSeqSet.empty())
    {}
    else if(vSeqSet.size() == 1)
    {
        strSeq = vSeqSet[0];
    }
    else
    {
        //2: Get the most frequency one directly, since the interupted has been filterd
        vector<char> vMostFreqChar;
        vector<int> vFreqValue;
        //We use the left alignment
        //Change All of Seq to Upper Case -->
        for(vector<string>::iterator itr = vSeqSet.begin(); itr != vSeqSet.end(); itr++)
        {
            std::transform(itr->begin(), itr->end(), itr->begin(), ::toupper);
        }
        //<--
        ClsAlgorithm::GetInstance().GetMostFreqChar(vSeqSet, vMostFreqChar, vFreqValue,
                                                    iMaxLen, enAlignType, false);
        //Set Value

        for(vector<char>::iterator itr = vMostFreqChar.begin(); itr != vMostFreqChar.end(); itr++)
        {
            if(*itr == '\0')
                break;
            strSeq += *itr;
        }
    }
    if(enAlignType == atLeft)
        strCmbSeq = strCmbSeq + strSeq;
    else // atRight
        strCmbSeq = strSeq + strCmbSeq;
}

void ClsRepeatBuild::GetAnchorInfo(string strFlank, string& strJF, string& strAnchor,
                                   int& iFreq, En_AnchorSide enSide)
{
    int iKmerNum = 16;
    int iCndNum = 5; // we need 5 different kemer for balance the kmer counting result
    //Transfer strFlank to uppercase --> Coincide with the search function later
    vector<string> vSubFlank;
    for(int i=0; i<iCndNum; i++)
    {
        if(i + iKmerNum > (int)strFlank.length())
            break;
        if(enSide == asLeft)
            vSubFlank.push_back(strFlank.substr(strFlank.length() - iKmerNum - i, iKmerNum));
        else if(enSide == asRight)
            vSubFlank.push_back(strFlank.substr(i, iKmerNum));
    }
    if(!vSubFlank.empty())
    {
        int iSum = 0;
        for(vector<string>::iterator itr = vSubFlank.begin(); itr != vSubFlank.end(); itr++)
        {
            //Get the value of kmer counting jellyfish
            iSum += GetKmerFreq(*itr, strJF);
        }
        iFreq = floor((float)iSum / (int)vSubFlank.size());
        strAnchor = vSubFlank[0]; // this time we think such flank with the high kmer frequency
    }
}

int ClsRepeatBuild::GetKmerFreq(string& strKmer, string& strJF)
{
    string strCmd = "";
    string strJellyfish = ClsAlgorithm::GetInstance().GetHigherFolderPath(ClsAlgorithm::GetInstance().GetCurExeFolderPath()) +
                          "ThirdPartyTools/Jellyfish/bin/jellyfish";
    string strOutput = ClsAlgorithm::GetInstance().GetHigherFolderPath(ClsAlgorithm::GetInstance().GetCurExeFolderPath()) +
                       "ThirdPartyTools/Jellyfish/data/log.ini";
    strCmd = strJellyfish + " query " + strJF + " " + strKmer + " > " + strOutput;
    //Option: "-C"
    //This result already contained the statistics result in both ordinary direction and its reverse comolementary
    system(strCmd.c_str());
    //Parse the number and return it
    ifstream ifs;
    ifs.open(strOutput.c_str());
    string strLine = "";
    getline(ifs, strLine);
    int iFreq = 0;
    if(strLine.find(" ") != string::npos &&       //Exsited
       strLine.find(" ") != strLine.length() - 1) //Not the last one
    {
        iFreq = atoi(strLine.substr(strLine.find(" ") + 1, strLine.length() - strLine.find(" ") - 1).c_str());
    }
    ifs.close();
    return iFreq;
}
