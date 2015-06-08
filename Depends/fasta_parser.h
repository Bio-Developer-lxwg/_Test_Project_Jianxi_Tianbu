#ifndef _H_RESRAP_ATSAF_
#define _H_RESRAP_ATSAF_

#include<string>
#include<vector>
#include<utility>
using namespace std;

string IntToStr(int iValue);
bool IsMissing(char nt);
char GetComplement(char bp);

struct St_Fasta
{
    string strName;
    string strSeq;

    St_Fasta():strName(""),strSeq("")
    {}

    ~St_Fasta(){}
};

class FastaParser
{
public:
	FastaParser();	

public:	
    //Step 1: Output the customized sequence ------------>lxwg
    void OutPutSubSeq(string strInputFastPath, string strOutputFastPath, int iStartPos, int iLength);
    int ReadFasta(string strPath, vector<St_Fasta>& vFasta);
    //<-----------------
    //Step 2:Get the random sequence
    void OutputRandomSeq(string strOutputPath, string strShortRepName, unsigned int iLength);
    void RevsereComplement(vector<St_Fasta>& vOrgFasta, vector<St_Fasta>& vDestFasta);
    string GetRevCompleString(string& strOrg);
};

#endif
