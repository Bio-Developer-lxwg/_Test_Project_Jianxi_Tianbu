//#include "clscorestructure.h"
#include "clsalgorithm.h"
#include <algorithm>
/*
 * The very basical structure of gap filling: it defines what the gap it is.
 */
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
        for(int i = iNPos + 1; i<(int)strOrg.length(); i++)
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
    if(iLastValidPos < (int)strOrg.length())
        strDest += strOrg.substr(iLastValidPos, strOrg.length() - iLastValidPos + 1);
}
//////////////////////////////////////////////////////////////////////////////////////

/* This structure focus on record all of Repeats info from the "repeat file"
 * We use the "St_RepeatUnit" as each Unit, and "St_RepeatFile" to manage the whole repeats
 */
//-------->The function of Structure
//NOTICE: we need to read the repeats at one time, for this is the very basical element for gap filling
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
        //--->Calculate the relevant Kmer--> For normal sequence
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
        //--->Calculate the relevant Kmer--> for reverse complementary
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
///////////////////////////////////////////////////////////////////////////////////////

/*The section of Scaffold structure identification:
 * 1: St_GapRefinedByRept: record the updated gap after filled by repeat
 * 2: St_ScaffoldUnit: this is the scaffold uint, designing for the single scaffold
 * 3: St_ScaffoldFile: manage all of scaffold
 */
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
            stNewGap.stBamGap.iStartPos = strRefinedSeq.length();
            strRefinedSeq += cGAPSYMBOL;
            stNewGap.stBamGap.iEndPos = strRefinedSeq.length() - 1;
            stNewGap.stBamGap.iLen = string(cGAPSYMBOL).length();
            stNewGap.pOrgGap = &(*itrGap);
        }
        else //Do not fix by repeat
        {
            stNewGap.stBamGap.iStartPos = strRefinedSeq.length();
            for(int i=0; i<itrGap->iLen; i++)
            {
                strRefinedSeq += "N";
            }
            stNewGap.stBamGap.iEndPos = strRefinedSeq.length() - 1;
            stNewGap.stBamGap.iLen = itrGap->iLen;
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
        InitKmerInfoForOneGap(itr->stBamGap, strRefinedSeq, 100); // in face we just conside the 30 characters: 1 kmer will be very good
    }
    //<---
}

//用于更新该scaffold的质量，通过考虑能够成功match到其上面的reads，以及成功产生的softcip的reads的数目，如果过少，很显然质量就不应该很高
//通过strRefinedSeq进行quality的判定
//这里一定要注意：实际上两边的softclip占的比例还挺大的
void St_ScaffoldUnit::UpdateQuality(int iMatechedReads, int iSoftClipReads)
{
    if(strRefinedSeq == "")
        return;
    int iN = 0;
    int iLastValidPos = 0;
    int iNPos = -1;
    while((iNPos = strRefinedSeq.find('N', iLastValidPos)) != (int)string::npos)
    {
        iN++;
        iLastValidPos = iNPos + 1;
    }
    int iValidChar = strRefinedSeq.length() - iN;
    float fBound1AnyMatch = iValidChar * .5;
    float fBound2AnyMatch = iValidChar * .25;
    // 50 意味着首位总公产生50个soft clip,然后我们考虑到整体的有效的字符数，然后取较小的那个值
    float fBoundPartMatch = (iSoftClipReads > 200 ? 200 : iSoftClipReads) * .5 +
                            50 > (iValidChar/4) ? (iValidChar/4) : 50;
    if(iMatechedReads >= fBound1AnyMatch)
    {
        if(iSoftClipReads > fBoundPartMatch)
            enQuality = sqHigh;
        else
            enQuality = sqNorm;
    }
    else if(iMatechedReads >= fBound2AnyMatch)
    {
        if(iSoftClipReads > fBoundPartMatch)
            enQuality = sqNorm;
        else
            enQuality = sqLow;
    }
    else
        enQuality = sqLow;
}

//////////////////////////St_ScaffoldFile/////////////////////////////////////////////////
void St_ScaffoldFile::Init(string strFilePath, FastaParser* pFastaParse, int iMaxRepLen)
{
    vector<St_Fasta> vFastaData;
    pFastaParse->ReadFasta(strFilePath, vFastaData);
    //Transfer vData to vector<St_ScaffoldUnit>
    //首先需要预估计相应的长度
    unsigned int iDataNum = 0;
    for(vector<St_Fasta>::iterator itrFasta = vFastaData.begin();
        itrFasta != vFastaData.end(); itrFasta++)
    {
        if(itrFasta->strName.find("scaffold") == string::npos) // 在名字里面没有scaffold,证明不是scaffold
            break;
        else
            iDataNum++;
    }
    cout << IntToStr(iDataNum) << endl;
    vData.resize(iDataNum);
    vector<St_Fasta>::iterator itrFasta = vFastaData.begin();
    int iIndex = 0;
    //Set value to each scaffold unit
    for(vector<St_ScaffoldUnit>::iterator itr = vData.begin();
        itr != vData.end() && itrFasta != vFastaData.end(); itr++, itrFasta++)
    {
        if(itrFasta->strName.find("scaffold") == string::npos) // 在名字里面没有scaffold,证明不是scaffold
            break;
        if(itrFasta->strName.find(' ') == string::npos) // discard the name after ' '
            itr->strName = itrFasta->strName;
        else
        {
            int iValidLen = itrFasta->strName.find(' ');
            itr->strName = itrFasta->strName.substr(0, iValidLen);
        }

        itr->strSeq = itrFasta->strSeq;
        itr->iID = iIndex;
        itr->FindGaps();
        //-->InitKmerInfoForEachGap
        for(vector<St_Gap>::iterator itrGap = itr->vGap.begin();
            itrGap != itr->vGap.end(); itrGap++)
        {
            itr->InitKmerInfoForOneGap(*itrGap, itr->strSeq, iMaxRepLen);
        }
        //<--
        cout << IntToStr(iIndex) << "/" << IntToStr(iDataNum) << endl;
        iIndex++;
        //if(iIndex > 5000) // 为了测试青蛙的真实数据，防治被kill做的临时性代码
        //   break;
    }
}
//<--------
//////////////////////////////////////////////////////////////////////////////////

/*the strucute relevant with fill the gap by pair-end reads
 *St_SoftClipReads: since currently, we just consider the softclip reads. that's why we buid this structure for calculation
 *St_FillResultBySCReads: Record the result after gap filling by softclip reads
 */
string St_SoftClipReads::GetExtdSeq()
{
    if(!bRevsStrand)
        return strExtendSeq;
    return ClsAlgorithm::GetInstance().GetReverseCompelement(strExtendSeq);
}
//////////////////////////////////////////////////////////
