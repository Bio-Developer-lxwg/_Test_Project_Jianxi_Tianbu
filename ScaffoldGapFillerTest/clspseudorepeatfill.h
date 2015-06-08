#ifndef CLSPSEUDOREPEATFILL_H
#define CLSPSEUDOREPEATFILL_H
#include "clsbasescaffoldfill.h"
#include "clscorestructure.h"
#include "string.h"
#include "clsrepeatbuild.h"
using namespace std;

class ClsPseudoRepeatFill:public ClsBaseScaffoldFill
{
public:
    ClsPseudoRepeatFill();
    virtual ~ClsPseudoRepeatFill();
    virtual void FillScafUnitWithSet(St_ScaffoldUnit& stScafUnit, vector<St_Fasta>& vScafContigSet,
                                     vector<St_ScaffoldUnit>& vScafSet);
    virtual void BaseInit();

public:
    void Init(string strRepeatPath, string strScafPath, string strReads1Path, string strReads2Path);
    int GetMaxRepeatLen();    

private:
    void FillGapByRepeat(St_RepeatUnit& stRepatUnit, St_Gap& stGap);
    int GetPosFromLocalAlignment(string& OrgStr, string& strCmp, En_GapType enType, int iTag = 0);
    void ActFillGapByRepeat(St_ScaffoldUnit& stScafUnit);

private:
     St_RepeatFile m_stRepeat;
     string m_strRepeatPath;
     ClsRepeatBuild*m_pRptBuild;
};

#endif // CLSPSEUDOREPEATFILL_H
