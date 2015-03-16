#ifndef CLSSCAFFOLDFILLER_H
#define CLSSCAFFOLDFILLER_H
#include "clscorestructure.h"
#include"fasta_parser.h"
#include "KmerUtils.h"
#include "api/BamReader.h"

using namespace BamTools;

struct St_CompEtdSeq //Compare the extend sequence --> it could be divided into the left extend and the right extend
{
    map<char, int> mpNdCount;
    char cValidNd;
    int iValidCount;
    St_CompEtdSeq():cValidNd('N'), iValidCount(-1)
    {}
};

class ClsScaffoldFiller
{
public:
    ClsScaffoldFiller();
    ~ClsScaffoldFiller();
public:
    void Init(char **argv);
    void FillScaffold();
    void FillScaffoldUnit(St_ScaffoldUnit& stScaffoldUnit);

    //Step 1:
    void FillScaffoldUnitByRepeats(St_ScaffoldUnit& stScaffoldUnit);
    bool EstimateFillGapByRepeat(St_RepeatUnit& stRepatUnit, St_Gap& stGap);
    void FillGapByRepeat(St_RepeatUnit& stRepatUnit, St_Gap& stGap);
    int GetPosFromLocalAlignment(string& OrgStr, string& strCmp, En_GapType enType, int iTag = 0);

    //Step 2:
    void FillScaffoldUnitByPairEndReads(St_ScaffoldUnit& stScaffoldUnit);
    string CreateBamFile(string& strFaName, string& strRefSeq);
    void ParseBamFileForGapFill(string strBamFilePath, St_ScaffoldUnit& stScaffoldUnit);
    void UpdateBamReads(BamAlignment& al, St_ScaffoldUnit& stScaffoldUnit,
                        int icurClipLen, int icurClipPos, int iRefPos);
    void FillGapBySoftClipReads();
    void AddValidBamReads(BamAlignment& al, St_GapRefinedByRept& stGapRefinedByRept,
                          int icurClipLen, int icurClipPos, int iRefPos);
    string CombineExt(string strExt1, string strExt2/*, int iStart1, int iStart2, int iEnd1, int iEnd2*/);   

    //Step 3: check the result to see if the filler is reasonable
    void CheckFillQuality(St_ScaffoldUnit& stScaffoldUnit);
    void QualityEstimate(string strBamFilePath, vector<St_FillUnit>& vFillUnit, string& strRef);

private:  
    //---->Reads Path
    string m_strReads1Path;
    string m_strReads2Path;
    //<----
    FastaParser* m_pFastaParse;
    St_RepeatFile m_stRepeat;
    St_ScaffoldFile m_stScaffold;
    St_BamReads m_stBamReads; //Just collect the valid bam reads
    vector<St_FinalScaffoldUnit> m_vFinalResult; // the final result
};

#endif // CLSSCAFFOLDFILLER_H
