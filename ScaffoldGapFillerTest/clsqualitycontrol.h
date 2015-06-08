#ifndef CLSQUALITYCONTROL_H
#define CLSQUALITYCONTROL_H
#include "clscorestructure.h"

enum En_AlnState{apLeft=0, apMiddle, apRight, apWhole, apMax};

class ClsQualityControl
{
public:
    ClsQualityControl();
    ~ClsQualityControl();

public:
    void Init(string strReads1Path, string strReads2Path);
    //1:对于Gap Size的评估
    void EstimateGapLength(St_ScaffoldFile& stScaffold, string strReads1Path,string strReads2Path);
    int CalGapLenghByOverlap(string& strLeftSeq, string& strRightSeq, int iLStart,
                             int iLEnd, int iRStart, int iREnd);
    bool EstimateFillGapByRepeat(St_RepeatUnit& stRepatUnit, St_Gap& stGap);

private:
    //Step 3: check the result to see if the filler is reasonable
    void CheckFillQuality(St_ScaffoldUnit& stScaffoldUnit);
    void QualityEstimate(string strBamFilePath, vector<St_FillUnit>& vFillUnit, string& strRef);

private:
    string m_strReads1Path;
    string m_strReads2Path;
};

#endif // CLSQUALITYCONTROL_H
