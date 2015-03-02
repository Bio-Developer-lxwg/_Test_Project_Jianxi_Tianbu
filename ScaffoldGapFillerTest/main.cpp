#include <iostream>
#include "clsscaffoldfiller.h"
#include "clsalgorithm.h"

using namespace std;
void GapFilling(string& strRepeatFilePath, string& strScaffoldFilePath);

bool CheckParameters(int argc, char **argv)
{
    if(argc != 3)
        return false;
    return true;
}

int main(int argc, char **argv)
{
    //Step 1: Check Parameters
    //Parameter 1: Repeat File
    //Parameter 2: scaffoldFile
    if(!CheckParameters(argc, argv))
        return -1;

    //Step 2:
    string strRepeatFilePath = argv[1];
    string strScaffoldFilePath =argv[2];
    GapFilling(strRepeatFilePath, strScaffoldFilePath);
    return 0;
}

void GapFilling(string& strRepeatFilePath, string& strScaffoldFilePath)
{
    ClsScaffoldFiller* pScaffoldFiller = new ClsScaffoldFiller();
    pScaffoldFiller->Init(strRepeatFilePath, strScaffoldFilePath);
    pScaffoldFiller->FillScaffold();

    delete pScaffoldFiller;
    pScaffoldFiller = NULL;
}
