#ifndef CLSALGORITHM_H
#define CLSALGORITHM_H
#include <string>
#include "clscorestructure.h"
using namespace std;

enum En_ScaffoldDataType{sdtNormal, sdtByRepeat, sdtByDraftGene, stdMax};
enum En_AlignType{atLeft=0, atRight, atAny};

bool sortfunction(St_SoftClipReads* pStI, St_SoftClipReads* pStJ);
//{ return pStI->iMissingStart < pStJ->iMissingStart; } // from small to large

bool sortfunctionLength(string str1, string str2);
//{ return str1.length() < str2.length(); } // from short to long

//Some Small Tool Functions
string::size_type FindStrPosByNS(string strBaseSeq, string strSeed, int iQueryPos = 0); // None Sensetive

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
    string CreateBamFile(string strFaName, string& strRefSeq,
                         string& strReads1Path, string& strReads2Path, bool bMapAll = true);
    string CreateBamFileByMultiScaff(vector<St_ScaffoldUnit>& vScaffoldUnit,
                                     string& strReads1Path, string& strReads2Path,                                     
                                     string strBamFileName = "", En_ScaffoldDataType enDataType = sdtNormal,
                                     bool bMapAll = true);
    //The relevant fasta file and fastq file has been created before
    string CreateBamFileForMultiRefPEReads(string& strRefPath, string& strRead1Path,
                                           string& strRead2Path, string strBamFileName = "", bool bMapAll = true);
    string CreateBamFileForMultiRefSingleReads(string& strRefPath, string& strReadsPath,
                                               string strBamFileName = "", bool bMapAll = true);
    //string calculation
    //This is just for DNA char, since we just considered 4 cases: A, T, G, C
    //MaxLen means the max length of string in the container of vOrgStr
    //The target of this function is try to find the max freqency char for each postion of a goup string and record its relevant frequency valye
    //在这里，我们默认是右对齐
    void GetMostFreqChar(vector<string>& vOrgStr, vector<char>& vMostFreqChar, vector<int>& vFreqValue,
                         int iMaxLen, En_AlignType enAlignType, bool bTight = true); // AlignType: 0 means 左对齐，1 means 右对齐
    St_RltBTRange GetTwoRangeRelationship(int iStart1, int iEnd1, int iStart2, int iEnd2); // Check the relationship between two ranges

    //statment checking
    bool IsWholeUpperCase(string& strValue);

    //record string in file ()
    void DisplayString(ofstream& ofs, string& strValue, int iLenPerLine=120);
};

#endif // CLSALGORITHM_H
