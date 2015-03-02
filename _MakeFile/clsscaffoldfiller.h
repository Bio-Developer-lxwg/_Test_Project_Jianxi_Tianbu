#ifndef CLSSCAFFOLDFILLER_H
#define CLSSCAFFOLDFILLER_H
#include"fasta_parser.h"
#include "KmerUtils.h"
#include "api/BamReader.h"

using namespace BamTools;

struct St_RepeatUnit
{
    string strName;
    string strNormalSeq;
    string strCompleRevSeq;
    int iKmerLen;
    vector<KmerTypeShort> vNormalKmerValue;
    vector<KmerTypeShort> vCompleRevKmerValue;

    St_RepeatUnit():iKmerLen(10)
    {}
};

struct St_RepeatFile
{
    vector<St_RepeatUnit> vData;
    void Init(string strFilePath, FastaParser* pFastaParse);
    unsigned int GetMaxLen();
    St_RepeatFile()
    {}
};

enum En_MatchType{mtNone=0, mtNormal, mtCompleRev, mtMax};
enum En_GapType{gtNone=0, gtLeft, gtCenter, gtRight, gtMax};
struct St_Gap
{
    int iStartPos; //in scaffold
    int iEndPos;   //in scaffold
    int iLen;

    float fRatioTolerent;
    int iKmerLen;
    string strLeftFlank;
    int iLeftFlankStart;
    string strRightFlank;
    int iRightFlankEnd;

    string strRefinedLeftFlank;
    string strRefinedRightFlank;
    vector<KmerTypeShort> vLeftKmer;
    vector<KmerTypeShort> vRightKmer;

    string strFillSeq;
    bool bFilled;
    En_MatchType enMatchType;
    En_GapType enGapType;
    int iReptStart; //in repeat (if is could be fill in repeat)
    int iReptEnd;

    bool CouldFill();
    void RefinedByN(string& strOrg, string& strDest);
    St_Gap():iStartPos(-1), iEndPos(-1), iLen(-1), fRatioTolerent(0.3), iKmerLen(10),
             strFillSeq(""), bFilled(false), enMatchType(mtNone), enGapType(gtNone),
             iReptStart(-1), iReptEnd(-1)
    {}
};

struct St_GapRefinedByRept
{
    St_Gap* pOrgGap;
    St_Gap stGap;
};

struct St_ScaffoldUnit
{
    string strName;
    string strSeq;
    vector<St_Gap> vGap; //Record Gap Start and End

    string strRefinedSeq;
    vector<St_GapRefinedByRept> vGapRefinedByRept;

    void FindGaps();
    void InitKmerInfoForOneGap(St_Gap& stGap, string& strRefSeq, int iMaxReptLen);
    void RefineFlank(St_Gap& stCurGap);
    void UpdateRefinedSeqByRepeats();
    bool IsMissing(char nt)
    {
        return nt == 'N' || nt == 'n';
    }
};

struct St_ScaffoldFile
{
    vector<St_ScaffoldUnit> vData;
    FastaParser* pFastaParse;
    void Init(string strFilePath, FastaParser* pFastaParse, int iMaxRepLen);
    St_ScaffoldFile()
    {}
};

////////////////////////我们现在来进行使用pair-end reads中的softclip来进行相应的gap filling 的extend
enum En_ClipPart {cpLeft, cpRight, cpNone, cpMax};
struct St_SoftClipReads
{
    string strSeq;
    string strClipSeq;
    int iClipStart;    
    int iClipLen;
    int iRefPos;
    En_ClipPart enClipPart;    
    bool bRevsStrand;
    bool bSecondMate;

    int iMissingStart;
    int iMissingEnd;

    int iAlignLeft;
    int iAlignRight;
    string strExtendSeq;

    St_GapRefinedByRept* pGap;

    string GetExtdSeq();

    St_SoftClipReads():iMissingStart(-1), iMissingEnd(-1), enClipPart(cpNone), pGap(NULL),
                       iAlignLeft(-1), iAlignRight(-1), bRevsStrand(false), bSecondMate(false)
    {}
};

struct St_FillResultBySCReads
{
    St_GapRefinedByRept* pGap;
    string strFillSeq;
    string strRevCompFillSeq;
    St_FillResultBySCReads(St_GapRefinedByRept* pGapValue, string strFillSeqValue, string strRevCompFillSeqValue)
    {
        pGap = pGapValue;
        strFillSeq = strFillSeqValue;        
        strRevCompFillSeq = strRevCompFillSeqValue;
    }
};

struct St_BamReads
{
    vector<St_SoftClipReads> vSCReads; //soft clip reads --> they are very useful!!!
    vector<St_FillResultBySCReads> vSCResult;
};

struct St_CompEtdSeq //Compare the extend  sequence --> it could be divided into the left extend and the right extend
{
    map<char, int> mpNdCount;
    char cValidNd;
    int iValidCount;
    St_CompEtdSeq():cValidNd('N'), iValidCount(-1)
    {}
};

enum En_FillQuality{fqNone = 0, fqLow, fqNormal, fqHigh, fqMax};
struct St_FillUnit
{
    int iStart;
    int iEnd;
    int iCoverage;
    En_FillQuality enQuality;

    St_FillUnit():iStart(-1), iEnd(-1), iCoverage(0), enQuality(fqNone)
    {}
};

struct St_FinalScaffoldUnit
{
    string strFillByNorm;
    string strFillByRevComp;
    string strFill;
    string strName;
    vector<St_FillUnit> vNormGapRecord;
    vector<St_FillUnit> vRevCompGapRecord;
    vector<St_FillUnit> vFinalFillGapRecord;
};

class ClsScaffoldFiller
{
public:
    ClsScaffoldFiller();
    ~ClsScaffoldFiller();
public:
    void Init(string strRepeatFilePath, string strScaffoldPath);
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
    FastaParser* m_pFastaParse;
    St_RepeatFile m_stRepeat;
    St_ScaffoldFile m_stScaffold;
    St_BamReads m_stBamReads; //Just collect the valid bam reads
    vector<St_FinalScaffoldUnit> m_vFinalResult; // the final result
};

#endif // CLSSCAFFOLDFILLER_H
