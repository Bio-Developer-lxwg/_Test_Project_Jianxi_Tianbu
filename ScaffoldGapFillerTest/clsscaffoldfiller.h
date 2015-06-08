#ifndef CLSSCAFFOLDFILLER_H
#define CLSSCAFFOLDFILLER_H
#include "clscorestructure.h"
#include"fasta_parser.h"
#include "KmerUtils.h"
#include "api/BamReader.h"

#include "clspseudorepeatfill.h"
#include "clsdraftgenofill.h"
#include "clspairendfill.h"
#include "clsqualitycontrol.h"

using namespace BamTools;

class ClsScaffoldFiller
{
public:
    ClsScaffoldFiller();
    ~ClsScaffoldFiller();
public:
    void Init(char **argv);
    void FillScaffold();
    //void FillScaffoldUnit(St_ScaffoldUnit& stScaffoldUnit);

private:  
    //---->Reads Path
    string m_strReads1Path;
    string m_strReads2Path;
    //<----
    FastaParser* m_pFastaParse;
    //St_RepeatFile m_stRepeat;
    St_ScaffoldFile m_stScaffold;
    St_Contigs m_stContigs;
    //St_BamFile m_stBamFile; //Just collect the valid bam reads
    vector<St_FinalScaffoldUnit> m_vFinalResult; // the final result
    vector<St_Fasta> m_vDraftGenoFaData; // include both scaffold and contig

//private:
//    void CheckFillingResult();
private:
    ClsPseudoRepeatFill* m_pClPseudoRpt; //Strategy 1: Pseudo Repeat Filling
    ClsDraftGenoFill* m_pClDraftGeno; //Strategy 2: Draft Geno Filling
    ClsPairEndFill* m_pClPEReads; //Strategy 3: Pair_End Reads Filling
    ClsQualityControl* m_pQualityCtr; //Strategy: Quality Control    
};

#endif // CLSSCAFFOLDFILLER_H
