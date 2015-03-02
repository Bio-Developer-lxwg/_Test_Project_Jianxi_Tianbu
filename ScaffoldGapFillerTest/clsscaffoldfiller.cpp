#include "clsscaffoldfiller.h"
#include "local_alignment.h"
#include "smith-waterman.h"
#include "stdlib.h"
#include <unistd.h>
#include "clsalgorithm.h"
#include "algorithm"

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

/////////////////St_Gap///////////////////
bool St_Gap::CouldFill()
{
    if(iLen < 10)
        return false;
    else
        return true;
}

void St_Gap::RefinedByN(string& strOrg, string& strDest)
{
    if(strOrg.find('N', 0) == string::npos)
    {
        strDest = strOrg;
        return;
    }
    strDest = "";
    int iNPos = -1;
    int iNNum = 0;
    int iLastValidPos = 0;
    while( (iNPos = strOrg.find('N', iLastValidPos)) != (int)string::npos)
    {
        iNNum = 1;
        for(int i = iNPos + 1; i<strOrg.length(); i++)
        {
            if(strOrg.at(i) == 'N')
                iNNum++;
            else
                break;
        }
        if(iNNum < 10)
            strDest += strOrg.substr(iLastValidPos, (iNPos - iLastValidPos));
        else
            strDest = "";
        iLastValidPos = iNPos + iNNum;
    }
    // add the last part
    if(iLastValidPos < strOrg.length())
        strDest += strOrg.substr(iLastValidPos, strOrg.length() - iLastValidPos + 1);
}
////////////////////////////St_RepeatFile/////////////////////////////////////////////
//-------->The function of Structure
void St_RepeatFile::Init(string strFilePath, FastaParser* pFastaParse)
{
    vector<St_Fasta> vFastaData;
    pFastaParse->ReadFasta(strFilePath, vFastaData);
    //Init Repeat -->
    vData.resize(vFastaData.size());
    vector<St_Fasta>::iterator itrFasta = vFastaData.begin();
    for(vector<St_RepeatUnit>::iterator itr = vData.begin();
        itr != vData.end() && itrFasta != vFastaData.end(); itr++, itrFasta++)
    {
        itr->strName = itrFasta->strName;
        itr->strNormalSeq = itrFasta->strSeq;
        //--->Calculate the relevant Kmer-->
        int iKmerNum = itr->strNormalSeq.length() - itr->iKmerLen + 1;
        if(iKmerNum <= 0)
            continue;
        itr->vNormalKmerValue.resize(iKmerNum);
        int iKmerPosStart = 0;
        for(vector<KmerTypeShort>::iterator itrKmer = itr->vNormalKmerValue.begin();
            itrKmer != itr->vNormalKmerValue.end(); itrKmer++)
        {
            FormKmerTypeShortSeg(itr->strNormalSeq.c_str(), iKmerPosStart, itr->iKmerLen, *itrKmer);
            iKmerPosStart++;
        }
        //<---

        itr->strCompleRevSeq = pFastaParse->GetRevCompleString(itr->strNormalSeq);
        //--->Calculate the relevant Kmer-->
        iKmerNum = itr->strCompleRevSeq.length() - itr->iKmerLen + 1;
        if(iKmerNum <= 0)
            continue;
        itr->vCompleRevKmerValue.resize(iKmerNum);
        iKmerPosStart = 0;
        for(vector<KmerTypeShort>::iterator itrKmer = itr->vCompleRevKmerValue.begin();
            itrKmer != itr->vCompleRevKmerValue.end(); itrKmer++)
        {
            FormKmerTypeShortSeg(itr->strCompleRevSeq.c_str(), iKmerPosStart, itr->iKmerLen, *itrKmer);
            iKmerPosStart++;
        }
        //<---

    }
}

unsigned int St_RepeatFile::GetMaxLen()
{
    if(vData.empty())
        return 0;

    unsigned int iMaxLen = 0;
    for(vector<St_RepeatUnit>::iterator itr = vData.begin(); itr != vData.end(); itr++)
    {
        if(iMaxLen < itr->strNormalSeq.length())
            iMaxLen = itr->strNormalSeq.length();
    }
    return iMaxLen;
}

//////////////////////////St_ScaffoldUnit/////////////////////////////////////////////////

void St_ScaffoldUnit::FindGaps()
{
    if("" == strSeq)
        return;
    vGap.clear();
    St_Gap stGap;
    for(unsigned int i=0; i<strSeq.length(); i++)
    {
        //
        if(IsMissing(strSeq.at(i) ))
        {
            // this is a gap
            stGap.iStartPos = i;
            while(IsMissing(strSeq.at(i) ) && i < strSeq.length())
            {
                i++;
            }
            stGap.iEndPos = i-1;
            stGap.iLen = stGap.iEndPos - stGap.iStartPos + 1;
            vGap.push_back(stGap);
        }
    }
}

void St_ScaffoldUnit::RefineFlank(St_Gap& stCurGap)
{
    //Consider the N, and the influence from other gaps
    //1:Refined by other gap -->We do not need to consider this case, since the following
    //logic will never cause this problem

    //2: Refine left flank
    stCurGap.RefinedByN(stCurGap.strLeftFlank, stCurGap.strRefinedLeftFlank);
    //3: Refine right flank
    stCurGap.RefinedByN(stCurGap.strRightFlank, stCurGap.strRefinedRightFlank);
}

void St_ScaffoldUnit::InitKmerInfoForOneGap(St_Gap& stGap, string& strRefSeq, int iMaxReptLen)
{
    //we need to set value for the following variantns
    //思路：我们按着MaxLen的方向去找，然后以发现其中有N作为中断的条件
    int iMaxLen = iMaxReptLen * stGap.fRatioTolerent;
    //Step 1:Find the strLeftValidSeq -->我们先找左边的
    int iLeftFlankStart = stGap.iStartPos - iMaxLen;
    int iLeftFlankLen = iMaxLen;
    if(iLeftFlankStart < 0)
    {
        iLeftFlankStart = 0;
        iLeftFlankLen = stGap.iStartPos;
    }
    stGap.strLeftFlank = strRefSeq.substr(iLeftFlankStart, iLeftFlankLen);
    stGap.iLeftFlankStart = iLeftFlankStart;

    //Step 2:Find the strRightValidSeq -->我们再找右边的
    int iRightFlankLen = iMaxLen;
    int iRightFlankEnd = stGap.iEndPos + iMaxLen;
    if(iRightFlankEnd >= (int)strRefSeq.length())
    {
        iRightFlankLen = strRefSeq.length() - stGap.iEndPos;
    }
    if(iRightFlankLen <= 0)
    {
        stGap.strRightFlank = "";
        stGap.iRightFlankEnd = stGap.iEndPos;
    }
    else
    {
        strRefSeq.substr(stGap.iEndPos+1, iRightFlankLen);
        stGap.iRightFlankEnd = stGap.iEndPos + iRightFlankLen;
    }

    //Step 3: Refine both left and right flank
    RefineFlank(stGap);

    //Step 4: 计算左边的kmer值
    if(stGap.strRefinedLeftFlank != "")
    {
        int iKmerNum = stGap.strRefinedLeftFlank.length() - stGap.iKmerLen + 1;
        stGap.vLeftKmer.resize(iKmerNum);
        int iKmerPosStart = 0;
        for(vector<KmerTypeShort>::iterator itrKmer = stGap.vLeftKmer.begin();
            itrKmer != stGap.vLeftKmer.end(); itrKmer++)
        {
            FormKmerTypeShortSeg(stGap.strRefinedLeftFlank.c_str(), iKmerPosStart, stGap.iKmerLen, *itrKmer);
            iKmerPosStart++;
        }
    }
    //step 5: 计算右边的kmer值
    if(stGap.strRefinedRightFlank != "")
    {
        int iKmerNum = stGap.strRefinedRightFlank.length() - stGap.iKmerLen + 1;
        stGap.vRightKmer.resize(iKmerNum);
        int iKmerPosStart = 0;
        for(vector<KmerTypeShort>::iterator itrKmer = stGap.vRightKmer.begin();
            itrKmer != stGap.vRightKmer.end(); itrKmer++)
        {
            FormKmerTypeShortSeg(stGap.strRefinedRightFlank.c_str(), iKmerPosStart, stGap.iKmerLen, *itrKmer);
            iKmerPosStart++;
        }
     }
}

void St_ScaffoldUnit::UpdateRefinedSeqByRepeats()
{
    strRefinedSeq = "";
    vGapRefinedByRept.clear();
    St_GapRefinedByRept stNewGap;
    if(vGap.empty())
    {
        strRefinedSeq = strSeq;
        return;
    }
    //1: Append the seq before the first gap
    strRefinedSeq += strSeq.substr(0, vGap[0].iStartPos);
    //2: Update the sequence by Gap
    const char* cGAPSYMBOL = "NNNNNNNNNNNNNNNNNNNN"; //20Ns
    for(vector<St_Gap>::iterator itrGap = vGap.begin(); itrGap != vGap.end(); itrGap++)
    {
        if(itrGap != vGap.begin())
        {
            //add the valid sequence between current and the last gap into refined seq
            strRefinedSeq += strSeq.substr((itrGap-1)->iEndPos + 1, (itrGap->iStartPos - (itrGap-1)->iEndPos - 1));
        }
        if(itrGap->bFilled) //1: if such gap has been fixed
        {
            strRefinedSeq += itrGap->strFillSeq;
            stNewGap.stGap.iStartPos = strRefinedSeq.length();
            strRefinedSeq += cGAPSYMBOL;
            stNewGap.stGap.iEndPos = strRefinedSeq.length() - 1;
            stNewGap.stGap.iLen = string(cGAPSYMBOL).length();
            stNewGap.pOrgGap = &(*itrGap);
        }
        else //Do not fix by repeat
        {
            stNewGap.stGap.iStartPos = strRefinedSeq.length();
            for(int i=0; i<itrGap->iLen; i++)
            {
                strRefinedSeq += "N";
            }
            stNewGap.stGap.iEndPos = strRefinedSeq.length() - 1;
            stNewGap.stGap.iLen = itrGap->iLen;
            stNewGap.pOrgGap = &(*itrGap);
        }
        vGapRefinedByRept.push_back(stNewGap); // now, we have collected the either refined and unrefined gap
    }       
    //3: Append the valid seq after the last gap
    unsigned int iLastIndex = vGap.size() - 1;
    strRefinedSeq += strSeq.substr(vGap[iLastIndex].iEndPos + 1, strSeq.length() - vGap[iLastIndex].iEndPos);

    //4: Since the string do not build over, we need to calculate the relevant kmer in the here (the last step)
    //--->We need to calculate the new kmer
    for(vector<St_GapRefinedByRept>::iterator itr = vGapRefinedByRept.begin();
        itr != vGapRefinedByRept.end(); itr++)
    {
        InitKmerInfoForOneGap(itr->stGap, strRefinedSeq, 100); // in face we just conside the 30 characters: 1 kmer will be very good
    }
    //<---
}

//////////////////////////St_ScaffoldFile/////////////////////////////////////////////////

void St_ScaffoldFile::Init(string strFilePath, FastaParser* pFastaParse, int iMaxRepLen)
{
    vector<St_Fasta> vFastaData;
    pFastaParse->ReadFasta(strFilePath, vFastaData);
    //Transfer vData to vector<St_ScaffoldUnit>
    unsigned int iDataNum = vFastaData.size();
    cout << IntToStr(iDataNum) << endl;
    vData.resize(iDataNum);
    vector<St_Fasta>::iterator itrFasta = vFastaData.begin();
    int iIndex = 0;
    for(vector<St_ScaffoldUnit>::iterator itr = vData.begin();
        itr != vData.end() && itrFasta != vFastaData.end(); itr++, itrFasta++)
    {       
        if(itrFasta->strName.find(' ') == string::npos) // discard the name after ' '
            itr->strName = itrFasta->strName;
        else
        {
            int iValidLen = itrFasta->strName.find(' ');
            itr->strName = itrFasta->strName.substr(0, iValidLen);
        }

        if(iIndex == 91748)
        {
            cout << itr->strName << endl;
        }

        itr->strSeq = itrFasta->strSeq;
        if(iIndex == 91748)
        {
            cout << "Start Find Gap" << endl;
        }
        itr->FindGaps();
        if(iIndex == 91748)
        {
            cout << "End Find Gap" << endl;
        }
        //-->InitKmerInfoForEachGap
        if(iIndex == 91748)
        {
            cout << "Start InitKmerInfo" << endl;
        }
        for(vector<St_Gap>::iterator itrGap = itr->vGap.begin();
            itrGap != itr->vGap.end(); itrGap++)
        {
            itr->InitKmerInfoForOneGap(*itrGap, itr->strSeq, iMaxRepLen);
        }
        if(iIndex == 91748)
        {
            cout << "End Init KmerInfo" << endl;
        }
        //<--
        cout << IntToStr(iIndex) << "//" << IntToStr(iDataNum) << endl;
        iIndex++;
    }
}
//<--------

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

void ClsScaffoldFiller::Init(string strRepeatFilePath, string strScaffoldPath)
{
    m_pFastaParse = new FastaParser();
    cout << "Start to init m_stRepeats" << endl;
    m_stRepeat.Init(strRepeatFilePath, m_pFastaParse);
    cout << "m_stRepeats init finished" << endl;
    cout << "Start to init m_stScaffold" << endl;
    m_stScaffold.Init(strScaffoldPath, m_pFastaParse, m_stRepeat.GetMaxLen());
    cout << "m_stScaffold Init finished" << endl;
}

void ClsScaffoldFiller::FillScaffold()
{
    if(m_stScaffold.vData.empty())
        return;
    m_vFinalResult.clear();
    //move the file
    string strResultPath = ClsAlgorithm::GetInstance().GetCurExeFolderPath() + "/GapFillerInfo";
    ::remove(strResultPath.c_str());
    //
    for(vector<St_ScaffoldUnit>::iterator itr = m_stScaffold.vData.begin();
        itr != m_stScaffold.vData.end(); itr++)
    {
        FillScaffoldUnit(*itr);
    }  
}

void ClsScaffoldFiller::FillScaffoldUnit(St_ScaffoldUnit& stScaffoldUnit)
{
    //Step 1: try to use repeat to fill the gap
    FillScaffoldUnitByRepeats(stScaffoldUnit);
    //Step 2: Use Pair end Reads to fix the gap:
    FillScaffoldUnitByPairEndReads(stScaffoldUnit);

    //--------->Just for temperay output: we try to print out those result
    //In this phase: we just want to list all the gaps which could be filled by repeats and pair-end reads
    ofstream ofs;
    string strRootPath = ClsAlgorithm::GetInstance().GetCurExeFolderPath();
    ofs.open((strRootPath + "/GapFillerInfo").c_str(), _S_app); //Do not need to use the "a+", since we write such file just one time
    {
        ofs << ">" << stScaffoldUnit.strName << endl;
        int i = 0;
        for(vector<St_FillResultBySCReads>::iterator itr = m_stBamReads.vSCResult.begin();
            itr != m_stBamReads.vSCResult.end(); itr++)
        {
            ofs << "Gap Index [" << IntToStr(i) << "]:" << endl;
            ofs << "Fill by Repeats ==>" << itr->pGap->pOrgGap->strFillSeq << endl;
            ofs << "Fill by PE Reads: Norm Ore ==> " << itr->strFillSeq << endl;
            ofs << "Fill by PE Reads: Reverse Complementary ==> " << itr->strRevCompFillSeq << endl;
            i++;
        }
    }
    //<----------------

    //Step 3: Check the Quality of fill result    
    //CheckFillQuality(stScaffoldUnit);
}

void ClsScaffoldFiller::FillScaffoldUnitByRepeats(St_ScaffoldUnit& stScaffoldUnit)
{
    //Go!!!!-->
    for(vector<St_Gap>::iterator itrGap = stScaffoldUnit.vGap.begin();
        itrGap != stScaffoldUnit.vGap.end(); itrGap++)
    {
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
                stGap.enMatchType = mtNormal;
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
                    stGap.enMatchType = mtCompleRev;
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
                stGap.enMatchType = mtNormal;
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
                    stGap.enMatchType = mtCompleRev;
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
    string& strOrg = (stGap.enMatchType == mtNormal) ? stRepatUnit.strNormalSeq : stRepatUnit.strCompleRevSeq;
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
    switch(stGap.enGapType)
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
    switch(enType)
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
void ClsScaffoldFiller::FillScaffoldUnitByPairEndReads(St_ScaffoldUnit& stScaffoldUnit)
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
    //Step 1: View current scaffold as the a reference seq and map the reads into this sequence
    //1:Create Bam File
    string strBamFilePath = CreateBamFile(stScaffoldUnit.strName, stScaffoldUnit.strRefinedSeq);
    //2:Analysis Bam Files for gap filling
    ParseBamFileForGapFill(strBamFilePath, stScaffoldUnit);
}

string ClsScaffoldFiller::CreateBamFile(string& strFaName, string& strRefSeq)
{
    cout << "Create Bam File" << endl;
    //1: Save such scaffold as the reference file
    //Get the exe path and the the fill path  for testing   ==>Try to seach internet!!!!
    string strRootPath = ClsAlgorithm::GetInstance().GetHigherFolderPath(get_current_dir_name());
    strRootPath += "TempFile/";
    //Clear the files under this folder
    string strCmd = "";//"rm -rf " + strRootPath + "*";
    //system(strCmd.c_str());
    //Set the valuefor Reference Fasta file
    string strRefPath = strRootPath + strFaName + ".fa";

    ofstream ofs;
    ofs.open(strRefPath.c_str());
    ofs << ">" << strFaName << endl;
    ofs << strRefSeq;
    ofs.close();

    //2: Build Index for this Reference Fa File
    strCmd = "bwa index -a bwtsw " + strRefPath;
    system(strCmd.c_str());

    //3: Alignment-->Build SAI
    string strRead1Path = strRootPath + "read1.fq";
    string strRead2Path = strRootPath + "read2.fq";
    string strSai1Path = strRootPath + "read1.sai";
    strCmd = "bwa aln -1 " + strRefPath + " " + strRead1Path + " > " + strSai1Path; //Read 1
    system(strCmd.c_str());
    string strSai2Path = strRootPath + "read2.sai";
    strCmd = "bwa aln -2 " + strRefPath + " " + strRead2Path + " > " + strSai2Path; //Read 1
    system(strCmd.c_str());

    //4: Generate SAM file
    string strSamPath = strRootPath + "Read.sam";
    strCmd = "bwa sampe -f " + strSamPath + " " +
             strRefPath + " " + strSai1Path + " " + strSai2Path + " " +
             strRead1Path + " " + strRead2Path;
    system(strCmd.c_str());

    //5: transfer sam to bam
    string strBamPath = strRootPath + "Read.bam";
    strCmd = "samtools view -bS " + strSamPath + " > " + strBamPath;
    system(strCmd.c_str());

    //6: Sort bam file
    string strSortedBamPath = strRootPath + "Read.sorted.bam";
    strCmd = "bamtools sort -in " +
            strBamPath + " -out " + strSortedBamPath;
    system(strCmd.c_str());

    //7:build index file for bam file
    strCmd = "bamtools index -in " + strSortedBamPath;
    system(strCmd.c_str());

    return strBamPath;
}

void ClsScaffoldFiller::ParseBamFileForGapFill(string strBamFilePath, St_ScaffoldUnit& stScaffoldUnit)
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

    //Step 1: Clear the old record
    m_stBamReads.vSCReads.clear();

    BamReader* pBamReader = new BamReader();
    pBamReader->Open(strBamFilePath);
    pBamReader->OpenIndex(strBamFilePath + ".bai");

    //step 1: record all of position of alignment seq
    BamAlignment al;
    vector<int> clipSizes;
    vector<int> readPositions;
    vector<int> genomePositions;
    while(pBamReader->GetNextAlignment(al))
    {        
        int iStart = al.Position;                
        if(iStart > 0) //这个时候应该是证明比上了
        {
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
                UpdateBamReads(al, stScaffoldUnit, icurClipLen, icurClipPos, iRefPos);
            }
        }
    }
    //use those valid softclip data for gap filling ----->
    FillGapBySoftClipReads();

    delete pBamReader;
    pBamReader = NULL;
}

void ClsScaffoldFiller::UpdateBamReads( BamAlignment& al, St_ScaffoldUnit&  stScaffoldUnit,
                                        int icurClipLen, int icurClipPos, int iRefPos)
{
    bool bValidReads = false;

    //Step1: we should check if this this clip is valid
    const int COFFSET = 5;
    for(vector<St_GapRefinedByRept>::iterator itr = stScaffoldUnit.vGapRefinedByRept.begin();
        itr != stScaffoldUnit.vGapRefinedByRept.end(); itr++)
    {
        //use the breakpoint for confirm if such reads could be used
        if(iRefPos-COFFSET > itr->stGap.iStartPos && iRefPos-COFFSET < itr->stGap.iEndPos ||
           iRefPos+COFFSET > itr->stGap.iStartPos && iRefPos+COFFSET < itr->stGap.iEndPos)
        {
            //Could be used then-->
            AddValidBamReads(al, *itr, icurClipLen, icurClipPos, iRefPos);
        }
    }
}

void ClsScaffoldFiller::AddValidBamReads(BamAlignment& al, St_GapRefinedByRept& stGapRefinedByRept,
                                         int icurClipLen, int icurClipPos, int iRefPos)
{
    //Type1: for softclip reads
    St_SoftClipReads stSCReads;
    stSCReads.iRefPos = iRefPos;
    stSCReads.iClipLen = icurClipLen;
    stSCReads.strSeq = al.AlignedBases;//QueryBases;
    stSCReads.iClipStart = icurClipPos;
    stSCReads.bRevsStrand = al.IsReverseStrand();// IsMateReverseStrand();
    stSCReads.bSecondMate = al.IsSecondMate();
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
    switch(stSCReads.enClipPart)
    {
        case cpLeft:
        {
            localAl.optAlign(stSCReads.pGap->stGap.strRefinedLeftFlank,
                        stSCReads.strClipSeq, iStartOrg, iEndOrg, iStartCmp, iEndCmp);
            int iLen = stSCReads.strClipSeq.length() - iEndCmp -
                       (iEndOrg - stSCReads.pGap->stGap.strRefinedLeftFlank.length());
            if(iEndCmp < stSCReads.strClipSeq.length()) // that means it could be extend
                stSCReads.strExtendSeq = stSCReads.strClipSeq.substr(iEndCmp, iLen);
            break;
        }
        case cpRight:
        {
            localAl.optAlign(stSCReads.pGap->stGap.strRefinedRightFlank,
                    stSCReads.strClipSeq, iStartOrg, iEndOrg, iStartCmp, iEndCmp);
            int iLen = iStartCmp - iStartOrg;
            if(iStartCmp > 0 && iLen > 0)
                stSCReads.strExtendSeq = stSCReads.strClipSeq.substr(0, iLen);
            break;
        }
    }
    m_stBamReads.vSCReads.push_back(stSCReads);
}

string St_SoftClipReads::GetExtdSeq()
{
    if(!bRevsStrand)
        return strExtendSeq;
    return ClsAlgorithm::GetInstance().GetReverseCompelement(strExtendSeq);
}

bool sortfunction (St_SoftClipReads* pStI, St_SoftClipReads* pStJ)
{ return pStI->iMissingStart < pStJ->iMissingStart; } // from small to large

void ClsScaffoldFiller::FillGapBySoftClipReads()
{
    if(m_stBamReads.vSCReads.empty())
        return;
    m_stBamReads.vSCResult.clear();

    map<St_GapRefinedByRept*, vector<St_SoftClipReads*> > mpGapReads;
    for(vector<St_SoftClipReads>::iterator itr = m_stBamReads.vSCReads.begin();
        itr != m_stBamReads.vSCReads.end(); itr++)
    {
        mpGapReads[itr->pGap].push_back(&(*itr));
    }
    //Parse the gap
    for(map<St_GapRefinedByRept*, vector<St_SoftClipReads*> >::iterator itr = mpGapReads.begin();
        itr != mpGapReads.end(); itr++)
    {
        //For current gap -->Try to fill it
        char* Type[2] = {"Left", "Right"};

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
            if((*itrSCReads)->bRevsStrand)
                continue;
            cout << Type[(*itrSCReads)->enClipPart] << ": " << (*itrSCReads)->strExtendSeq << endl;
            switch((*itrSCReads)->enClipPart)
            {
                case cpLeft:
                {
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
                case cpRight:
                {
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
            }
        }
        //Decide if we need to make the reverse complementay to express the final result
        //bool bRevComplementary = false; //not good
        //if(iRevComplOre > iNormOre)
        //    bRevComplementary = true;

        //Try to Merge the "left" Extend sequence ==========================>        
        string strLeftExtend, strRightExtend;
        if(iMaxExtendNumLeft > 0)
        {
            vector<St_CompEtdSeq> vCompLeftSeq;
            vCompLeftSeq.resize(iMaxExtendNumLeft);
            for(vector<string>::iterator itr = vExtendFromLeft.begin();
                itr != vExtendFromLeft.end(); itr++)
            {
                for(int i = 0; i<itr->length(); i++)
                {
                    int iIndex = iMaxExtendNumLeft -(itr->length() - i);
                    vCompLeftSeq[iIndex].mpNdCount[itr->at(i)]++;
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
                for(int i = 0; i<itr->length(); i++)
                {
                    vCompRightSeq[i].mpNdCount[itr->at(i)]++;
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

        cout << "Left Extend Sum is: " << strLeftExtend << endl;
        cout << "Right Extend Sum is: " << strRightExtend << endl;

        // 合并左flank和右flank
        string strExtend;        
        strExtend = CombineExt(strLeftExtend, strRightExtend/*, iStartOrg, iStartCmp, iEndOrg, iEndCmp*/);        
        //---------->        
        cout << "Additional Extend: " << strExtend << endl;
        string strRevComplment = ClsAlgorithm::GetInstance().GetReverseCompelement(strExtend);
        cout << "Additional Extend(Reverse Complementary): " << strRevComplment << endl;
        m_stBamReads.vSCResult.push_back(St_FillResultBySCReads(itr->first, strExtend, strRevComplment));
    }
}

enum En_StringSide{ssLeft=0, ssRight, ssCenter, ssMax};
//this is a recursive
string ClsScaffoldFiller::CombineExt(string strExt1, string strExt2/*, int iStart1, int iStart2, int iEnd1, int iEnd2*/)
{
    //我们可以考虑一下-->global alignment --> 这样说不定会有比较好的结果 -->今天到此为止吧 --->一会儿去看看python去，完成一下那个老师的作业

    LocalAlignment localAl;//--> 使用现有的库来进行Local Alignment
    int iStart1 = -1, iStart2 = -1, iEnd1 = -1, iEnd2 = -1;
    localAl.optAlign(strExt1, strExt2, iStart1, iEnd1, iStart2, iEnd2);

    En_StringSide enExt1Side, enExt2Side;
    //case 0: 有一个是全匹配，那么纠结，直接最长的那个就是结果
    //For enExt1Side  // 这个的返回值是按照1为start的，并不是0 这个要注意！！！！！， 我么可以把他们搞成0；
    iStart1--;
    iEnd1--;
    iStart2--;
    iEnd2--;
    if( (iStart1 == 0 && iEnd1 == strExt1.length() - 1) ||
        (iStart2 == 0 && iEnd2 == strExt2.length() - 1))
    {
        return strExt1.length() >= strExt2.length() ? strExt1 : strExt2;
    }

    //Make sure the statement of cut side
    //For enExt1Side
    if(iStart1 > 0 && iEnd1 < strExt1.length() - 1)
        enExt1Side = ssCenter;
    else if(iStart1 > 0 && iEnd1 == strExt1.length() - 1)
        enExt1Side =ssLeft;
    else if(iStart1 == 0 && iEnd1 < strExt1.length() - 1)
        enExt1Side =ssRight;
    else {}
    //For enExt2Side
    if(iStart2 > 0 && iEnd2 < strExt2.length() - 1)
        enExt2Side = ssCenter;
    else if(iStart2 > 0 && iEnd2 == strExt2.length() - 1)
        enExt2Side =ssLeft;
    else if(iStart2 == 0 && iEnd2 < strExt2.length() - 1)
        enExt2Side =ssRight;
    else {}

    //Let's Compare
    switch(enExt1Side)
    {
    case ssLeft:
    {
        string strSubExt1 = strExt1.substr(0, iStart1);
        string strCmbExt1 = strExt1.substr(iStart1, iEnd1 - iStart1 + 1);
        switch(enExt2Side)
        {
        case ssLeft:
        {
            string strSubExt2 = strExt2.substr(0, iStart2);
            string strCmbExt2 = strExt2.substr(iStart2, iEnd2 - iStart2 + 1);
            return CombineExt(strSubExt1, strSubExt2) +
                   (strCmbExt1.length() >= strCmbExt2.length() ? strCmbExt1 : strCmbExt2);
        }
        case ssRight:
        {
            string strSubExt2 = strExt2.substr(iEnd2, strSubExt2.length() - iEnd2);
            string strCmbExt2 = strExt2.substr(iStart2, iEnd2 - iStart2 +1);
            return strSubExt1 +
                   (strCmbExt1.length() >= strCmbExt2.length() ? strCmbExt1 : strCmbExt2) +
                   strSubExt2;
        }
        case ssCenter:
        {
            return "";
            //do not consider
            int i = 0;
            break;
        }
        }
        break;
    }
    case ssRight:
    {
        string strSubExt1 = strExt1.substr(iEnd1, strExt1.length() - iEnd1);
        string strCmbExt1 = strExt1.substr(iStart1, iEnd1 - iStart1 + 1);
        switch(enExt2Side)
        {
        case ssLeft:
        {
            string strSubExt2 = strExt2.substr(0, iStart2);
            string strCmbExt2 = strExt2.substr(iStart2, iEnd2 - iStart2 + 1);
            return strSubExt2 +
                   (strCmbExt1.length() >= strCmbExt2.length() ? strCmbExt1 : strCmbExt2) +
                   strSubExt1;
        }
        case ssRight:
        {
            string strSubExt2 = strExt2.substr(iEnd2, strSubExt2.length() - iEnd2);
            string strCmbExt2 = strExt2.substr(iStart2, iEnd2 - iStart2 +1);
            return (strCmbExt1.length() >= strCmbExt2.length() ? strCmbExt1 : strCmbExt2) +
                   CombineExt(strSubExt1, strSubExt2);
            break;
        }
        case ssCenter:
        {
            return "";
            // do nothing
            break;
        }
        }
        break;
    }
    case ssCenter:
    {
        return "";
        //do nothing
        string strSubExt1Left = strExt1.substr(0, iStart1);
        string strSubExt1Right = strExt1.substr(iEnd1, strExt1.length() - iEnd1);
        switch(enExt2Side)
        {
        case ssLeft:
        {
            break;
        }
        case ssRight:
        {
            break;
        }
        case ssCenter:
        {
            break;
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
        //Notice: the gap cannot be filled will be kept
        //Check each  gaps
        for(vector<St_FillResultBySCReads>::iterator itrSCR = m_stBamReads.vSCResult.begin(); //SCR: soft clip reads
            itrSCR != m_stBamReads.vSCResult.end(); itrSCR++)
        {
            if(itrSCR->pGap == &(*itr)) //got it
            {
                //Add the repeat
                stFillUnit.iStart = stFinalScaffoldUnit.strFillByNorm.length() - 1;
                stFinalScaffoldUnit.strFillByNorm += itrSCR->pGap->pOrgGap->strFillSeq; //Fill by Repeat
                stFinalScaffoldUnit.strFillByNorm += itrSCR->strFillSeq; //Fill by Paired-End reads
                if(itrSCR->pGap->pOrgGap->strFillSeq == "" && itrSCR->strFillSeq == "")
                {
                    //if both filling method are failed to fill any sequences--->Add the gap back to it.
                    for(int i=0; i<itrSCR->pGap->pOrgGap->iLen; i++)
                    {
                        stFinalScaffoldUnit.strFillByNorm += "N";
                    }
                }
                stFinalScaffoldUnit.strFillByRevComp += itrSCR->pGap->pOrgGap->strFillSeq; //Fill by Repeat
                stFinalScaffoldUnit.strFillByRevComp += itrSCR->strRevCompFillSeq;
                if(itrSCR->pGap->pOrgGap->strFillSeq == "" && itrSCR->strRevCompFillSeq == "")
                {
                    //if both filling method are failed to fill any sequences--->Add the gap back to it.
                    for(int i=0; i<itrSCR->pGap->pOrgGap->iLen; i++)
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
                if(itrSCR + 1 < m_stBamReads.vSCResult.end())
                {
                    iLen = (itrSCR+1)->pGap->pOrgGap->iStartPos - itrSCR->pGap->pOrgGap->iEndPos - 1;
                    iStartPos = itrSCR->pGap->pOrgGap->iEndPos+1;
                }
                else // the last one then
                {
                    iLen = stScaffoldUnit.strSeq.length() - itrSCR->pGap->pOrgGap->iEndPos - 1;
                    iStartPos = itrSCR->pGap->pOrgGap->iEndPos+1;
                }
                stFinalScaffoldUnit.strFillByNorm += stScaffoldUnit.strSeq.substr(iStartPos, iLen);
                stFinalScaffoldUnit.strFillByRevComp += stScaffoldUnit.strSeq.substr(iStartPos, iLen);
                break;
            }
        }
    }
    //Step2: Get Map Read Back to Refined Scaffold to check the coverage in gao area
    //1: Clear the folder rm -rf path/*  do not needed, since this function has been added into bam file creation
    //2: For the Norm Sequence
    //  (1) Make the bam file
    //  (2) Analysis the result
    string strBamFilePath = CreateBamFile(stFinalScaffoldUnit.strName, stFinalScaffoldUnit.strFillByNorm);
    QualityEstimate(strBamFilePath, stFinalScaffoldUnit.vNormGapRecord,
                    stFinalScaffoldUnit.strFillByNorm);
    //3: For the reverse complementary Sequence
    //  (1) Make the bam file
    //  (2) Analysis the result
    strBamFilePath = CreateBamFile(stFinalScaffoldUnit.strName, stFinalScaffoldUnit.strFillByRevComp);
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

