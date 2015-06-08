#ifndef CLSCORESTRUCTURE_H
#define CLSCORESTRUCTURE_H
#include"fasta_parser.h"
#include "KmerUtils.h"
#include "api/BamReader.h"
using namespace BamTools;

/*
 * The very basical structure of gap filling: it defines what the gap it is.
 */
enum En_MatchType{mtNone=0, mtNormal, mtCompleRev, mtMax};
enum En_GapType{gtNone=0, gtLeft, gtCenter, gtRight, gtMax};
enum En_FillResult{frNone=0, frAllPass, frLeftFail, frRighFail, frAllFail, frMax}; //Fill result identification

//Basical Structure --> Gneral Math Algorithm
enum En_RangeRltsp{rrDepart=0, rrLeftOverlap, rrRightOverlap, rrContain, rrAbnormal, rrMax}; // The relationship of two ranges
struct St_RltBTRange //the struct of relationship between two ranges
{
    En_RangeRltsp enRangeType;
    int iOverlapStart;
    int iOverlapEnd;
    int iGapIndex;

    St_RltBTRange():enRangeType(rrAbnormal), iOverlapStart(-1), iOverlapEnd(-1), iGapIndex(-1)
    {}

    //-->Function
    int GetOverlapLen() // get the overlap length
    {
        return iOverlapEnd - iOverlapStart + 1;
    }
    //<--
};

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
    En_MatchType enLFMatchType; //Left Flank
    En_MatchType enRFMatchType; //Right Flank
    En_GapType enGapType;
    int iReptStart; //in repeat (if is could be fill in repeat)
    int iReptEnd;

    //评估之后的长度
    int iEstimateLen;
    bool bAbnormalLarge;

    //定义最后的结果是否已经被成功填补完成
    En_FillResult enFillResult;

    bool CouldFill();
    void RefinedByN(string& strOrg, string& strDest);
    St_Gap():iStartPos(-1), iEndPos(-1), iLen(-1), fRatioTolerent(0.3), iKmerLen(10),
             strFillSeq(""), bFilled(false), enLFMatchType(mtNone), enRFMatchType(mtNone),
             enGapType(gtNone), iReptStart(-1), iReptEnd(-1), iEstimateLen(-1), bAbnormalLarge(false),
             enFillResult(frNone)
    {}
};

struct St_BamGap: St_Gap  //Design for the gap of bam file (after refined by repeats )
{
    string strRevCompFillSeq;
    string strLeftExt;
    string strRightExt;
};
////////////////////////////////////////////////////////////////////////

/* This structure focus on record all of Repeats info from the "repeat file"
 * We use the "St_RepeatUnit" as each Unit, and "St_RepeatFile" to manage the whole repeats
 */
struct St_RepeatUnit
{
    string strName;
    string strNormalSeq; // fa文件中存储的序列
    string strCompleRevSeq; //fa文件中存储序列的反向逆序列 -->目的是为了解决可能反向逆才是真正的gap相应的填充序列，“但是目前没有使用”
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
////////////////////////////////////////////////////////////////////////
/*
 * we need to identify the type of third party tools to restrict parsing method
 */
enum En_AssemblyTool{atSoapDenovo=0, atMax};

/*The section of Draft Geno 1: Scaffold structure identification:
 * 1: St_GapRefinedByRept: record the updated gap after filled by repeat
 * 2: St_ScaffoldUnit: this is the scaffold uint, designing for the single scaffold
 * 3: St_ScaffoldFile: manage all of scaffold
 */
struct St_GapRefinedByRept
{
    St_Gap* pOrgGap;
    St_BamGap stBamGap;
    St_GapRefinedByRept(): pOrgGap(NULL)
    {}
    St_GapRefinedByRept(St_Gap* pGapValue, string strFillSeqValue, string strRevCompFillSeqValue)
    {
        pOrgGap = pGapValue;
        stBamGap.strFillSeq = strFillSeqValue;
        stBamGap.strRevCompFillSeq = strRevCompFillSeqValue;
    }
};

enum En_ScaffoldQuality {sqNone=0, sqLow, sqNorm, sqHigh, sqMax};

struct St_CPoint
{
    int iStart;
    int iEnd;
    St_CPoint(): iStart(-1), iEnd(-1)
    {}

    St_CPoint(int iStartValue, int iEndValue)
    {
        iStart = iStartValue;
        iEnd = iEndValue;
    }
};

struct St_OverLapState
{
    St_CPoint stBorder;
    St_RltBTRange stRltBTRange;
    //Default constructure
    St_OverLapState()
    {}
    //special construture
    St_OverLapState(St_CPoint stPointValue, St_RltBTRange stRltBTRangeValue)
    {
        stBorder = stPointValue;
        stRltBTRange = stRltBTRangeValue;
    }
};

struct St_MappingResult
{
    int iCount;
    vector<St_OverLapState> vMapReads; // Just Start and End

    St_MappingResult(): iCount(0)
    {}
};

struct St_FillingBorder
{
    int iStart;
    int iEnd;
    int iBKPEAndRepeats; // break point between Pair-End Reads and Repeats
    //Default Constucture
    St_FillingBorder(): iStart(-1), iEnd(-1), iBKPEAndRepeats(-1)
    {}
    //Customized Constructure
    St_FillingBorder(int iStartValue, int iEndValue, int iBKPEAndRepeatsValue)
    {
        iStart = iStartValue;
        iEnd = iEndValue;
        iBKPEAndRepeats = iBKPEAndRepeatsValue;
    }
};

struct St_FinalFilling // 以一个scaffold为单位进行考虑
{
    vector<St_FillingBorder> vBorder; //This is the gap border
    string strSeq; //The final filling sequence: notice we do not add any NNNNN into it.
    vector<int> vCoverage; // the coverage calculated by PE reads for each gap

    //Function
    void QualityCheckByPEReads(string& strScafName,
                               string& strReads1Path, string& strReads2Path); // 使用pair end reads，检查填补的质量
    St_RltBTRange GetOverlapStates(int iMapStartPos, int iMapEndPos, int iBaseLen, bool bRealCase = false);
};

struct St_ScaffoldUnit
{
    string strName;
    string strSeq;
    int iID; //Used for Bam File parsing
    vector<St_Gap> vGap; //Record Gap Start and End

    string strRefinedSeq;
    vector<St_GapRefinedByRept> vGapRefinedByRept;

    En_ScaffoldQuality enQuality;

    //---->Record the final filling result
    St_FinalFilling stFinalRst; // Final Filling Result
    //<----

    void FindGaps();
    void InitKmerInfoForOneGap(St_Gap& stGap, string& strRefSeq, int iMaxReptLen);
    void RefineFlank(St_Gap& stCurGap);
    void UpdateRefinedSeqByRepeats();
    //Update Final Filling Result-->
    void CheckFillingResult(string& strReads1Path, string& strReads2Path);
    void UpdateFinalSeq();
    //<--
    bool IsMissing(char nt)
    {
        return nt == 'N' || nt == 'n';
    }
    //calc quality-->
    void UpdateQuality(int iMatechedReads, int iSoftClipReads);
    //<--
    St_ScaffoldUnit():iID(-1000), enQuality(sqNone)
    {}
};

struct St_ScaffoldFile
{
    vector<St_ScaffoldUnit> vData;
    En_AssemblyTool enAsmbTool; //Assembly Tool
    //FastaParser* pFastaParse;
    void Init(vector<St_Fasta>& vFastaData, int iMaxRepLen);
    St_ScaffoldFile(): enAsmbTool(atSoapDenovo)
    {}
};

/*
 * Section of Draft Geno 2: Contigs
 */
struct St_ContigUnit
{
    string strName;
    string strSeq;
    int iID;
    St_ContigUnit():strName(""), strSeq(""), iID(-1)
    {}
};

struct St_Contigs
{
    En_AssemblyTool enAsmbTool; //Assembly Tool
    vector<St_ContigUnit> vData;

    void Init(vector<St_Fasta>& vFastaData);
    St_Contigs(): enAsmbTool(atSoapDenovo)
    {}
};
////////////////////////////////////////////////////////////////////////////////////////

/*the strucute relevant with fill the gap by pair-end reads
 * St_SoftClipReads:       since currently, we just consider the softclip reads. that's why we buid this structure for calculation
 * St_FillResultBySCReads: Record the result after gap filling by softclip reads
 * St_BamFile:            manage the bam file which contains all of target softclip reads
 */
////////////////////////我们现在来进行使用pair-end reads中的softclip来进行相应的gap filling 的extend
enum En_ClipPart {cpLeft, cpRight, cpNone, cpMax};
struct St_SoftClipReads
{
    string strSeq; // Aligned Part
    string strClipSeq; // cliped part

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

    St_SoftClipReads():enClipPart(cpNone),  bRevsStrand(false),  bSecondMate(false),
                       iMissingStart(-1), iMissingEnd(-1), iAlignLeft(-1), iAlignRight(-1), pGap(NULL)
    {}    
};

/*
struct St_BamFile
{
    vector<St_SoftClipReads> vSCReads; //soft clip reads --> they are very useful!!!
    //vector<St_FillResultBySCReads> vSCResult;
    //vector<St_GapRefinedByRept*> vpGapAfterRept; //The Gap after the phase of repeats filling
};*/

//////////////////////////////////////////////////////////////////////////

/* The structure which used to record and express the final result
 * St_FillUnit:          each result unit of final gap filling result
 * St_FinalScaffoldUnit: manage the final results
 * */
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
///////////////////////////////////////////////////////////////////////////
/// \brief The St_GapReads struct
struct St_NoneMap
{
    int iMatPos;
    St_NoneMap():iMatPos(-1)
    {}
    St_NoneMap(int iMatPosValue)
    {
        iMatPos = iMatPosValue;
    }
};

struct St_GapReads
{
    vector<St_SoftClipReads> vRightClip;
    vector<St_SoftClipReads> vLeftClip;
    //vector<St_NoneMap> vNoneMap;
    int iLength;
    bool bAbnormalLarge;
    St_GapReads():iLength(0),bAbnormalLarge(false)
    {}
};

struct St_GapAlnSet //This could be considered as one scaffold
{
    vector<St_GapReads> vGapReadsSet;
    //bool bAddNoneMap;
    St_GapAlnSet()//:bAddNoneMap(false)
    {}
};

#endif // CLSCORESTRUCTURE_H
