#ifndef CLSREPEATBUILD_H
#define CLSREPEATBUILD_H
#include "clscorestructure.h"
#include "clsalgorithm.h"

enum En_AnchorSide {asLeft=0, asRight, asMax};

class ClsRepeatBuild
{
public:
    ClsRepeatBuild();
    ~ClsRepeatBuild();
    void Init(string strScafPath, string strReads1Path, string strReads2Path);

public:
    //Gap fill by abnormal Pair-End Reads
    void GapFillByAbnPE(St_ScaffoldUnit& stScafUnit, vector<St_Fasta>& vScafContigSet);

    //The most important part in my program
public:
    void BuildPseudoRepeat(St_ScaffoldUnit& stScaf, vector<St_Fasta>& vDraftSet,
                           vector<St_ScaffoldUnit>& vScafSet); // we do not need the info of pair-end reads currently
private:
    string KmerCounting(string strDraftGeno);
    void FillCurGapByDraftGeno(St_Gap& stGap, string& strJF, vector<St_Fasta>& vDraftSet,
                               vector<St_ScaffoldUnit>& vScafSet);
    int GetKmerFreq(string& strKmer, string& strJF);
    void GetAnchorInfo(string strFlank, string& strJF, string& strAnchor, int& iFreq, En_AnchorSide enSide);

    void AnchorLeftExtend(string strAnchor, int iFreq,
                          vector<St_Fasta>& vDraftSet, string& strJF, string& strCmb);

    void AnchorRightExtend(string strAnchor, int iFreq,
                           vector<St_Fasta>& vDraftSet, string& strJF, string& strCmb);

    void AnchorExtdByBWA(string strAnchor, int iFreq, vector<St_Fasta>& vDraftSet);
    string FindRemainSeqByOrg(string& strSeq, string& strAnchor);    
    void GetCombinedSeq(string& strCmbSeq, vector<string>& vSeqSet, string strJF,
                        int iFreq, En_AlignType enAlignType);
    //Return iCurFreq ===>
    //int GetBreakPos(string strSubStr, string& strCmbSeq, string strJF,
    //                int iFreq, En_AlignType enAlignType); // Kmer & SubSeq

private:
    void FillCurGapByPEReads(St_ScaffoldUnit& stScafUnit, St_Gap& stGap, vector<St_Fasta>& vScafContigSet);    
    bool FillByCopyFrag(string strAnchor, vector<St_Fasta>& vDraftSet, vector<St_ScaffoldUnit>& vScafSet,
                        string& strCmb, En_AnchorSide enAnchorSide); //Fill the current gap by the copy of other filling result

private:
    string m_strBamPath;
    string m_strScafPath;
    string m_strReads1Path;
    string m_strReads2Path;    
};

#endif // CLSREPEATBUILD_H
