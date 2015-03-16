#include "clsalgorithm.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <vector>
#include <cctype>
#include <algorithm>
#include "needleman_wunsch.h"
#include <dirent.h>
#include <iostream>
#include <fstream>

#define MAXBUFSIZE 260
using namespace std;

string ClsAlgorithm::GetCurExePath()
{
    char buf[MAXBUFSIZE];
    memset(buf, 0, MAXBUFSIZE);
    int iCount;
    iCount = readlink("/proc/self/exe", buf, MAXBUFSIZE);
    if(iCount < 0 || iCount >= MAXBUFSIZE) {
        printf( "Failed\n" );
        return "";
    }
    return buf;
}

string ClsAlgorithm::GetCurExeFolderPath()
{
    //::get_current_dir_name();
    char buf[MAXBUFSIZE];
    memset(buf, 0, MAXBUFSIZE);
    ::getcwd(buf, MAXBUFSIZE);
    return buf;
}

string ClsAlgorithm::GetHigherFolderPath(string strCurPath, int iLevel)
{
    vector<int> vLevel;
    int iPos = 0;
    while((iPos = strCurPath.find('/', iPos)) != string::npos)
    {        
        vLevel.push_back(iPos);
        iPos += 1;
    }
    if((int)vLevel.size() < iLevel)
        return "";
    if(iLevel == 0)
        return  strCurPath.at(strCurPath.length()-1) != '/' ? strCurPath + "/" : strCurPath;
    return strCurPath.substr(0, vLevel[vLevel.size()-iLevel]+1);
}

bool ClsAlgorithm::IsMissing(char nt)
{
    return nt == 'N' || nt == 'n';
}

char ClsAlgorithm::GetComplement(char bp)
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

string ClsAlgorithm::GetReverseCompelement(string strOrg, bool bRevsCompelement)
{
    if(!bRevsCompelement) //Do not need the reverse complementary sequence
        return strOrg;
    //wanna the reverse complementary sequence
    //1: Get reverse
    reverse(strOrg.begin(), strOrg.end());
    //2: Get Compement
    for(unsigned int i=0; i<strOrg.length(); i++)
    {
        strOrg[i] = GetComplement(strOrg[i]);
    }
    return strOrg;
}

void ClsAlgorithm::GlobalAlignment(char* seq_a, char* seq_b)
{
    // Variables to store alignment result
    char *alignment_a, *alignment_b;

    // malloc the above variables
    // (seq1 and seq2 are used to figure out how much memory may be needed)
    nw_alloc_mem(seq_a, seq_b, &alignment_a, &alignment_b);

    // Decide on scoring
    int match = 1;
    int mismatch = -1;//-2;
    int gap_open = 0;//-4;
    int gap_extend = 0;//-1;

    // Don't penalise gaps at the start
    // ACGATTT
    // ----TTT would score +3 (when match=+1)
    char no_start_gap_penalty = 1;

    // ..or gaps at the end e.g.
    // ACGATTT
    // ACGA--- would score +4 (when match=+1)
    char no_end_gap_penalty = 1;

    // Compare character case-sensitively (usually set to 0 for DNA etc)
    char case_sensitive = 0;

    SCORING_SYSTEM* scoring = scoring_create(match, mismatch,
                                             gap_open, gap_extend,
                                             no_start_gap_penalty,
                                             no_end_gap_penalty,
                                             case_sensitive);

    // Add some special cases
    // x -> y means x in seq1 changing to y in seq2
    scoring_add_mutation(scoring, 'a', 'c', -2); // a -> c give substitution score -2
    scoring_add_mutation(scoring, 'c', 'a', -1); // c -> a give substitution score -1

    // We could also prohibit the aligning of characters not given as special cases
    // scoring->use_match_mismatch = 0;

    int score = needleman_wunsch(seq_a, seq_b, alignment_a, alignment_b, scoring);

    printf("seqA: %s\n", alignment_a);
    printf("seqB: %s\n", alignment_b);
    printf("alignment score: %i\n", score);

    // Free memory used to store scoring preferences
    scoring_free(scoring);

    free(alignment_a);
    free(alignment_b);
}

//Return intersection number
int ClsAlgorithm::CheckInterSection(int iStart1, int iEnd1, int iStart2, int iEnd2)
{
    if(iStart1 > iEnd1 ||
       iStart2 > iEnd2)
    {
        return 0;
    }
    if(iEnd2 < iStart1 || iStart2 > iEnd1)
        return 0;
    else
    {
        if(iStart2 <= iStart1)
        {
            if(iEnd2 <= iEnd1)
            {
                return iEnd2 - iStart1 + 1;
            }
            else
                return iEnd1 - iStart1 + 1;
        }
        else if(iStart2 >= iStart1 && iStart2 <= iEnd1)
        {
            if(iEnd2 <= iEnd1)
                return iEnd2 - iStart2 + 1;
            else
                return iEnd1 - iStart2 + 1;
        }
    }
    return 0;
}

string ClsAlgorithm::CreateBamFile(string& strFaName, string& strRefSeq,
                                   string& strReads1Path, string& strReads2Path)
{
    cout << "Create Bam File" << endl;
    //1: Save such scaffold as the reference file
    //Get the exe path and the the fill path  for testing   ==>Try to seach internet!!!!
    string strRootPath = ClsAlgorithm::GetInstance().GetHigherFolderPath(get_current_dir_name());
    strRootPath += "TempFile/";
    //Clear the files under this folder
    string strCmd = "";//"rm -rf " + strRootPath + "*";
    //system(strCmd.c_str());
    //Set the valuefor Reference Fasta file
    string strRefPath = strRootPath + strFaName + ".fa";

    ofstream ofs;
    ofs.open(strRefPath.c_str());
    ofs << ">" << strFaName << endl;
    ofs << strRefSeq;
    ofs.close();

    //2: Build Index for this Reference Fa File
    strCmd = "bwa index -a bwtsw " + strRefPath;
    system(strCmd.c_str());

    //3: Alignment-->Build SAI
    //string strRead1Path = m_strReads1Path//strRootPath + "read1.fq";
    //string strRead2Path = strRootPath + "read2.fq";
    string strSai1Path = strRootPath + "read1.sai";
    strCmd = "bwa aln -1 " + strRefPath + " " + strReads1Path + " > " + strSai1Path; //Read 1
    system(strCmd.c_str());
    string strSai2Path = strRootPath + "read2.sai";
    strCmd = "bwa aln -2 " + strRefPath + " " + strReads2Path + " > " + strSai2Path; //Read 2
    system(strCmd.c_str());

    //4: Generate SAM file
    string strSamPath = strRootPath + "Read.sam";
    strCmd = "bwa sample -f " + strSamPath + " " +
             strRefPath + " " + strSai1Path + " " + strSai2Path + " " +
             strReads1Path + " " + strReads2Path;
    system(strCmd.c_str());

    //5: transfer sam to bam
    string strBamPath = strRootPath + "Read.bam";
    strCmd = "samtools view -bS " + strSamPath + " > " + strBamPath;
    system(strCmd.c_str());

    //6: Sort bam file
    string strSortedBamPath = strRootPath + "Read.sorted.bam";
    strCmd = "bamtools sort -in " +
            strBamPath + " -out " + strSortedBamPath;
    system(strCmd.c_str());

    //7:build index file for bam file
    strCmd = "bamtools index -in " + strSortedBamPath;
    system(strCmd.c_str());

    return strSortedBamPath;
}

string ClsAlgorithm::CreateBamFileByMultiScaff(vector<St_ScaffoldUnit>& vScaffoldUnit,
                                               string& strReads1Path, string& strReads2Path,
                                               string strBamFileName, En_ScaffoldDataType enDataType)
{
    cout << "Create Bam File" << endl;
    //1: Save such scaffold as the reference file
    //Get the exe path and the the fill path  for testing   ==>Try to seach internet!!!!
    string strRootPath = ClsAlgorithm::GetInstance().GetHigherFolderPath(get_current_dir_name());
    strRootPath += "TempFile/";
    //Clear the files under this folder
    string strCmd = "";//"rm -rf " + strRootPath + "*";
    //system(strCmd.c_str());
    //Set the valuefor Reference Fasta file
    string strFaName = "ScaffoldSet";
    string strRefPath = strRootPath + strFaName + ".fa";

    ofstream ofs;
    ofs.open(strRefPath.c_str());
    for(vector<St_ScaffoldUnit>::iterator itr = vScaffoldUnit.begin(); itr != vScaffoldUnit.end(); itr++)
    {
        ofs << ">" << itr->strName << endl;
        switch((int)enDataType)
        {
            case sdtNormal:
                ofs << itr->strSeq << endl;
                break;
            case sdtByRepeat:
            case sdtByDraftGene:
                ofs << itr->strRefinedSeq << endl;
                break;
        }
    }
    ofs.close();

    //2: Build Index for this Reference Fa File
    strCmd = "bwa index -a bwtsw " + strRefPath;
    system(strCmd.c_str());

    //3: Alignment-->Build SAI
    //string strRead1Path = m_strReads1Path//strRootPath + "read1.fq";
    //string strRead2Path = strRootPath + "read2.fq";
    /*string strSai1Path = strRootPath + "read1.sai";
    strCmd = "bwa aln -1 " + strRefPath + " " + strReads1Path + " > " + strSai1Path; //Read 1
    system(strCmd.c_str());
    string strSai2Path = strRootPath + "read2.sai";
    strCmd = "bwa aln -2 " + strRefPath + " " + strReads2Path + " > " + strSai2Path; //Read 2
    system(strCmd.c_str());*/

    //4: Generate SAM file
    string strSamPath = strRootPath + (strBamFileName == "" ? "Read.sam" : (strBamFileName + ".sam"));
    strCmd = "bwa mem " + strRefPath + " " + strReads1Path + " " + strReads2Path + " > " + strSamPath;
    system(strCmd.c_str());

    //5: transfer sam to bam
    string strBamPath = strRootPath + (strBamFileName == "" ? "Read.bam" : (strBamFileName + ".bam"));
    strCmd = "samtools view -bS " + strSamPath + " > " + strBamPath;
    system(strCmd.c_str());

    //6: Sort bam file
    string strSortedBamPath = strRootPath +
                              (strBamFileName == "" ? "Read.sorted.bam" :
                                                      (strBamFileName + ".sorted.bam"));
    strCmd = "bamtools sort -in " +
            strBamPath + " -out " + strSortedBamPath;
    system(strCmd.c_str());

    //7:build index file for bam file
    strCmd = "bamtools index -in " + strSortedBamPath;
    system(strCmd.c_str());

    //8: delete the bamfile and samfile since they are too large to be stored
    //a) Delete Sam File
    strCmd = "rm -f " + strSamPath;
    system(strCmd.c_str());
    //b) Delete Bam File
    strCmd = "rm -f " + strBamPath;
    system(strCmd.c_str());

    return strSortedBamPath;
}


