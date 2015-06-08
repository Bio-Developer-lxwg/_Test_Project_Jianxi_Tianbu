#ifndef CLSDRAFTGENOFILL_H
#define CLSDRAFTGENOFILL_H
#include "clsbasescaffoldfill.h"

struct St_FlankPair
{
    string strLeft;
    string strRight;
    string strName; //Use GapIndex to Express its name
    St_FlankPair():strLeft(""), strRight(""), strName("")
    {}
};

class ClsDraftGenoFill: public ClsBaseScaffoldFill
{
public:    
    ClsDraftGenoFill();
    virtual ~ClsDraftGenoFill();
    virtual void FillScafUnit(St_ScaffoldUnit& stScafUnit);
    virtual void BaseInit();
    void Init(vector<St_Fasta>& vFastaData);

private:
    void CalcByBWA(St_ScaffoldUnit& stScafUnit);

private:
    St_Contigs m_stContigs;
};

#endif // CLSDRAFTGENOFILL_H
