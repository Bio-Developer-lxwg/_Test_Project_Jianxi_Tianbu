#ifndef CLSALGORITHM_H
#define CLSALGORITHM_H
#include <string>
using namespace std;

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
};

#endif // CLSALGORITHM_H
