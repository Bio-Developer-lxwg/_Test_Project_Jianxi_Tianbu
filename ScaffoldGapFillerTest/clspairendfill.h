#ifndef CLSPAIRENDFILL_H
#define CLSPAIRENDFILL_H
#include "clsbasescaffoldfill.h"

struct St_CompEtdSeq //Compare the extend sequence --> it could be divided into the left extend and the right extend
{
    map<char, int> mpNdCount;
    char cValidNd;
    int iValidCount;
    St_CompEtdSeq():cValidNd('N'), iValidCount(-1)
    {}
};

class ClsPairEndFill: public ClsBaseScaffoldFill
{
public:    
    ClsPairEndFill();
    virtual ~ClsPairEndFill();
    virtual void FillScafUnit(St_ScaffoldUnit& stScafUnit);
    virtual void FillScafSet(vector<St_ScaffoldUnit>& vScafSet);
    virtual void BaseInit();

public:
    void Init(string strReads1Path, string strReads2Path);

private:
    void ParseBamFileForGapFill(string strBamFilePath, vector<St_ScaffoldUnit>& vScaffoldSet);
    void UpdateBamReads(BamAlignment& al, St_ScaffoldUnit&  stScaffoldUnit,
                        int icurClipLen, int icurClipPos, int iRefPos,
                        vector<St_SoftClipReads>& vSCReads);
    void FillGapBySoftClipReads(vector<St_SoftClipReads>& vSCReads);
    void AddValidSCReads(BamAlignment& al, St_GapRefinedByRept& stGapRefinedByRept,
                         int icurClipLen, int icurClipPos, int iRefPos,
                         vector<St_SoftClipReads>& vSCReads);
    string CombineExtByWholeExt(string strLeftExt, string strRightExt,
                                vector<string>& vWholeLeftExt, vector<string>& vWholeRightExt);
    string CombineExt(string strExt1, string strExt2/*, int iStart1, int iStart2, int iEnd1, int iEnd2*/);

private:
    string m_strReads1Path;
    string m_strReads2Path;
};

#endif // CLSPAIRENDFILL_H
