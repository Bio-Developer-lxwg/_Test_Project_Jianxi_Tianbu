#include <iostream>
#include "clsscaffoldfiller.h"
#include "clsalgorithm.h"
#include "unistd.h"

using namespace std;
void GapFilling(char **argv);

const char* czFileType[] = {"Repeats File", "Scaffold File", "Reads_1", "Reads_2"};
bool CheckParameters(int argc, char **argv)
{
    if(argc != 5) //0 is itself, 1 --> 4 are the files need to be used for gap filling
        return false;
    //check if the file is exsited
    for(int i=1; i<5; i++) // 1-->4
    {
        if(::access(argv[i], 0) != 0)
        {
            cout << "Error: " << czFileType[i] << " is not exsited!" << endl;
            return false;
        }
    }
    return true;
}

int main(int argc, char **argv)
{
    //Step 1: Check Parameters
    //------------------->
    //Parameter 1: Repeat File
    //Parameter 2: scaffoldFile
    //Parameter 3: reads1
    //Parameter 4: reads2
    //<-------------------
    if(!CheckParameters(argc, argv))
        return -1;

    //Step 2:
    GapFilling(argv);
    return 0;
}

void GapFilling(char **argv)
{
    ClsScaffoldFiller* pScaffoldFiller = new ClsScaffoldFiller();
    cout << "Start Init" <<endl;
    pScaffoldFiller->Init(argv);
    cout << "Init Finished" << endl;
    pScaffoldFiller->FillScaffold();

    delete pScaffoldFiller;
    pScaffoldFiller = NULL;
}
