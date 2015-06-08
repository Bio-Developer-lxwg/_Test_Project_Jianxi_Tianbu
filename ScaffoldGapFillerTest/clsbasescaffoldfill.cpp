#include "clsbasescaffoldfill.h"

ClsBaseScaffoldFill::ClsBaseScaffoldFill(): m_pQualityCtr(NULL), m_pFastaParse(NULL)
{
}

ClsBaseScaffoldFill::~ClsBaseScaffoldFill()
{
    if(m_pQualityCtr)
    {
        delete m_pQualityCtr;
        m_pQualityCtr = NULL;
    }

    if(m_pFastaParse)
    {
        delete m_pFastaParse;
        m_pFastaParse = NULL;
    }
}

void ClsBaseScaffoldFill::FillScafUnit(St_ScaffoldUnit& stScafUnit)
{}

void ClsBaseScaffoldFill::FillScafUnitWithSet(St_ScaffoldUnit& stScafUnit, vector<St_Fasta>& vScafContigSet,
                                              vector<St_ScaffoldUnit>& vScafSet)
{}

void ClsBaseScaffoldFill::FillScafSet(vector<St_ScaffoldUnit>& vScafSet)
{}

void ClsBaseScaffoldFill::BaseInit()
{
    if(m_pQualityCtr == NULL)
    {
        m_pQualityCtr = new ClsQualityControl();
    }
    if(m_pFastaParse == NULL)
    {
        m_pFastaParse = new FastaParser();
    }
}
