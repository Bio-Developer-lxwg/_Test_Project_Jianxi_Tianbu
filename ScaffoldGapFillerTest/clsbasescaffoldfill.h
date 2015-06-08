#ifndef CLSBASESCAFFOLDFILL_H
#define CLSBASESCAFFOLDFILL_H
#include "clscorestructure.h"
#include "clsqualitycontrol.h"

class ClsBaseScaffoldFill
{
public:
    ClsBaseScaffoldFill();
    virtual ~ClsBaseScaffoldFill();

public:
    virtual void FillScafUnit(St_ScaffoldUnit& stScafUnit);
    virtual void FillScafUnitWithSet(St_ScaffoldUnit& stScafUnit, vector<St_Fasta>& vScafContigSet,
                                     vector<St_ScaffoldUnit>& vScafSet);
    virtual void FillScafSet(vector<St_ScaffoldUnit>& vScafSet);
    virtual void BaseInit() = 0;

protected:
    ClsQualityControl* m_pQualityCtr; //Strategy: Quality Control
    FastaParser* m_pFastaParse;
};

#endif // CLSBASESCAFFOLDFILL_H
