#ifndef CLSALGORITHM_H
#define CLSALGORITHM_H
#include <string>
#include "clscorestructure.h"
using namespace std;

enum En_ScaffoldDataType{sdtNormal, sdtByRepeat, sdtByDraftGene, stdMax};

class ClsAlgorithm
{
public:
    static ClsAlgorithm& GetInstance()
    {
        static ClsAlgorithm instance;
        return instance;
    }
private:
    ClsAlgorithm(){}
    ClsAlgorithm(const ClsAlgorithm&);
    ClsAlgorithm& operator = (const ClsAlgorithm&);

public:
    //Folder Operation
    string GetCurExePath();
    string GetCurExeFolderPath();
    string GetHigherFolderPath(string strCurPath, int iLevel=1);
    //String Operation
    char GetComplement(char bp);
    bool IsMissing(char nt);
    string GetReverseCompelement(string strOrg, bool bRevsCompelement = true);
    //Bio Algorithm
    void GlobalAlignment(char* seq_a, char* seq_b);
    //Set
    int CheckInterSection(int iStart1, int iEnd1, int iStart2, int iEnd2);
    //Bam File
    string CreateBamFile(string& strFaName, string& strRefSeq,
                         string& strReads1Path, string& strReads2Path);
    string CreateBamFileByMultiScaff(vector<St_ScaffoldUnit>& vScaffoldUnit,
                                     string& strReads1Path, string& strReads2Path,                                     
                                     string strBamFileName = "", En_ScaffoldDataType enDataType = sdtNormal);
    //string calculation
    //This is just for DNA char, since we just considered 4 cases: A, T, G, C
    //MaxLen means the max length of string in the container of vOrgStr
    //The target of this function is try to find the max freqency char for each postion of a goup string and record its relevant frequency valye
    //在这里，我们默认是右对齐
    void GetMostFreqChar(vector<string>& vOrgStr, vector<char>& vMostFreqChar, vector<int>& vFreqValue, int iMaxLen, int iAlignType); // AlignType: 0 means 左对齐，1 means 右对齐
};

#endif // CLSALGORITHM_H
