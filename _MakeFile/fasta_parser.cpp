#include"fasta_parser.h"
#include<fstream>
#include "unistd.h"
#include <iostream>
#include <cstdlib>
#include <string.h>

using namespace std;

string IntToStr(int iValue)
{
    char czStr[20];
    ::memset(czStr, 0, 20);
    ::sprintf(czStr, "%d", iValue);
    return czStr;
}

bool IsMissing(char nt)
{
    return nt == 'N' || nt == 'n';
}

char GetComplement(char bp)
{
    //
    if( IsMissing(bp) )
    {
        return bp;
    }
    char bpUse = toupper(bp);
    if( bpUse == 'A')
    {
        return 'T';
    }
    //else if(bp == 'a')
    //{
    //    return 't';
    //}
    else if( bpUse == 'T')
    {
        return 'A';
    }
    //else if( bpUse == 't')
    //{
    //    return 'a';
    //}
    else if( bpUse == 'G')
    {
        return 'C';
    }
    //else if(bp == 'g')
    //{
    //    ;
    //}
    else if(bpUse == 'C')
    {
        return 'G';
    }
    return 'N';
}

FastaParser::FastaParser()
{}

int FastaParser::ReadFasta(string& strPath, vector<St_Fasta>& vFasta)
{
    if(::access(strPath.c_str(), 0) != 0)
    {
        cout << "Contig Generated Failed" << endl;
        return 1;
    }
    vFasta.clear();
    //Parse the contig file
    ifstream infile;
    infile.open(strPath.c_str(), ios::in);
    string strLine = "";
    St_Fasta stFasta;
    while(!infile.eof()) // check if reached the end
    {
        getline(infile, strLine);
        string::size_type sztp = strLine.find('>', 0);
        if(sztp != string::npos) // this line is the name of contig
        {
            if(stFasta.strName != "")
            {
                vFasta.push_back(stFasta);
                stFasta.strName = "";
                stFasta.strSeq = "";
            }
            //Set Name
            stFasta.strName = strLine.substr(1, strLine.length()-1);
            strLine = "";
            continue;
        }
        //Set Sequence
        if(infile.eof()) // for the last item
        {
            vFasta.push_back(stFasta);
        }
        else
        {
            stFasta.strSeq += strLine;
        }        
    }
    infile.close();
    strLine = "";
    return -1;
}

//Output the customized sequence
void FastaParser::OutPutSubSeq(string strInputFastaPath, string strOutputFastPath, int iStartPos, int iLength)
{
    vector<St_Fasta> vFasta;
    this->ReadFasta(strInputFastaPath, vFasta);
    if(vFasta.empty())
    {
        cout << "Parse Failed" << endl;
        return;
    }
    if((int)vFasta[0].strSeq.length() < iStartPos + iLength)
    {
        cout << "Invalid customized value of start pos and length" << endl;
        return;
    }
    //Output file
    string strName = vFasta[0].strName;
    string strSeq = "";
    unsigned int iOffset = 0;
    while( vFasta[0].strSeq.length() > iOffset && vFasta[0].strSeq[iOffset] == 'N')
    {
        iOffset++;
    }
    if(iStartPos + iLength + iOffset >= vFasta[0].strSeq.length())
    {
        cout << "Error cased by N" << endl;
        return;
    }

    strSeq = vFasta[0].strSeq.substr(iStartPos+iOffset, iLength);
    ofstream ofs(strOutputFastPath.c_str());
    if(ofs.is_open())
    {
        ofs << ">" << strName << endl;
        unsigned int iRow = 1;
        while(iRow*50 < strSeq.length()) //output 50 charactors per line
        {
            ofs << strSeq.substr((iRow-1)*50, 50) << endl;
            iRow++;
        }
        if((iRow-1)*50 < strSeq.length())
            ofs << strSeq.substr((iRow-1)*50, strSeq.length() - (iRow-1)*50) << endl;
    }
    ofs.close();
    cout << "Generated succesfully" << endl;
    return;
}

void FastaParser::OutputRandomSeq(string strOutputPath, string strShortRepName, unsigned int iLength)
{
    srand(time(NULL));
    string strSeq = "";
    unsigned int iTime = 0;
    while(iTime < iLength)
    {
        switch(rand()%4)
        {
            case 0:
                strSeq += "A";
                break;
            case 1:
                strSeq += "T";
                break;
            case 2:
                strSeq += "G";
                break;
            case 3:
                strSeq += "C";
                break;
            default:
                break;
        }
        iTime++;
    }
    if(strSeq.length() != iLength)
    {
        cout << "Wrong Shot Repeats generat" <<endl;
        return;
    }
    ofstream ofs(strOutputPath.c_str());
    if(ofs.is_open())
    {
        ofs << strShortRepName << endl;
        unsigned int iRow = 1;
        while(iRow*50 < iLength) //output 50 charactors per line
        {
            ofs << strSeq.substr((iRow-1)*50, 50) << endl;
            iRow++;
        }
        if((iRow-1)*50 < iLength)
            ofs << strSeq.substr((iRow-1)*50, iLength - (iRow-1)*50) << endl;
    }
}

void FastaParser::RevsereComplement(vector<St_Fasta>& vOrgFasta, vector<St_Fasta>& vDestFasta)
{
    // do the reverse complement
    vDestFasta.clear();
    vDestFasta.resize(vOrgFasta.size());
    vector<St_Fasta>::iterator itrOrgFa = vOrgFasta.begin();
    for(vector<St_Fasta>::iterator itr = vDestFasta.begin();
        itr != vDestFasta.end() && itrOrgFa != vOrgFasta.end(); itr++, itrOrgFa++)
    {
        itr->strName = itrOrgFa->strName;
        //Reverse Complementary the org seq
        for(int i=itrOrgFa->strSeq.length()-1; i>=0; --i)
        {
            char bpComp = GetComplement(itrOrgFa->strSeq[i]);
            itr->strSeq.push_back(bpComp);
        }
    }
}

string FastaParser::GetRevCompleString(string& strOrg)
{
    // do the reverse complement
    string strNew;
    for(int i=(int)strOrg.length()-1; i>=0; --i)
    {
        char bpComp = GetComplement(strOrg[i]);
        strNew.push_back(bpComp);
    }
    // assign the new string
    return strNew;
}
