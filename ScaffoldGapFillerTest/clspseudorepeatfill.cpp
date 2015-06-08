#include "clspseudorepeatfill.h"
#include "local_alignment.h"

ClsPseudoRepeatFill::ClsPseudoRepeatFill()
{
    m_strRepeatPath = "";
    m_pRptBuild = NULL;
}

ClsPseudoRepeatFill::~ClsPseudoRepeatFill()
{
    if(m_pRptBuild)
    {
        delete m_pRptBuild;
        m_pRptBuild = NULL;
    }
}

void ClsPseudoRepeatFill::BaseInit()
{}

void ClsPseudoRepeatFill::Init(string strRepeatPath, string strScafPath, string strReads1Path, string strReads2Path)
{
    //Init its parent class
    BaseInit();
    //Init the ariables of itself
    m_strRepeatPath = strRepeatPath;
    m_stRepeat.Init(m_strRepeatPath, m_pFastaParse); //Init Repeat File
    if(m_pRptBuild == NULL)
    {
        m_pRptBuild = new ClsRepeatBuild();
        m_pRptBuild->Init(strScafPath, strReads1Path, strReads2Path);
    }    
}

int ClsPseudoRepeatFill::GetMaxRepeatLen()
{
    return m_stRepeat.GetMaxLen();
}

void ClsPseudoRepeatFill::FillScafUnitWithSet(St_ScaffoldUnit& stScafUnit, vector<St_Fasta>& vScafContigSet,
                                              vector<St_ScaffoldUnit>& vScafSet)
{
    //--->
    m_pRptBuild->BuildPseudoRepeat(stScafUnit, vScafContigSet, vScafSet);
    //m_pRptBuild->GapFillByAbnPE(stScafUnit, vScafContigSet);
    //<---
    //--->Fill the Gap if there is the repeat reference file
    //ActFillGapByRepeat(stScafUnit);

    //We need to update the refined sequence
    stScafUnit.UpdateRefinedSeqByRepeats();
}

void ClsPseudoRepeatFill::ActFillGapByRepeat(St_ScaffoldUnit& stScafUnit)
{
    //Go!!!!-->
    for(vector<St_Gap>::iterator itrGap = stScafUnit.vGap.begin();
        itrGap != stScafUnit.vGap.end(); itrGap++)
    {
        if(!itrGap->bAbnormalLarge) //do not fill by repeats if the length is not abnormal large
            continue;

        if(!itrGap->CouldFill())
            continue;

        //Possible fill then ---->
        //Use Normal Repeat and complement reverse repeat to fill the Gap
        //(both those two type has been contained by St_RepeatUnit)
        for(vector<St_RepeatUnit>::iterator itr = m_stRepeat.vData.begin();
            itr != m_stRepeat.vData.end(); itr++)
        {
            //For each repeat --> We need to check if it could be used to fix the gap
            if(itrGap->bFilled)
                continue;
            if(m_pQualityCtr->EstimateFillGapByRepeat(*itr, *itrGap))
            {
                FillGapByRepeat(*itr, *itrGap);
                if(itrGap->bFilled)
                {
                    cout<<"=================Gap Filler===============" << endl;
                    cout<<"Scaffold Name: " << stScafUnit.strName << endl;
                    cout<<"Filled String: " << itrGap->strFillSeq << endl;
                    cout<<"From the \'Start\' Pos of Repeat: " << IntToStr(itrGap->iReptStart) << endl;
                    cout<<"From the \'End\' Pos of Repeat: " << IntToStr(itrGap->iReptEnd) << endl;
                    cout<<"Filled \'Length\': " << IntToStr(itrGap->strFillSeq.length()) << endl;
                }
            }
        }
    }
}

void ClsPseudoRepeatFill::FillGapByRepeat(St_RepeatUnit& stRepatUnit, St_Gap& stGap)
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
            iRightFillPos = GetPosFromLocalAlignment(strOrg, stGap.strRefinedRightFlank, gtLeft);
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

int ClsPseudoRepeatFill::GetPosFromLocalAlignment(string& OrgStr, string& strCmp, En_GapType enType, int iTag)
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


