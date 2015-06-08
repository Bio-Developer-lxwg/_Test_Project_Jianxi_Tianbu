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
#include <math.h>
#define MAXBUFSIZE 260
using namespace std;


bool sortfunction (St_SoftClipReads* pStI, St_SoftClipReads* pStJ)
{
    return pStI->iMissingStart < pStJ->iMissingStart;
} // from small to large

bool sortfunctionLength (string str1, string str2)
{
    return str1.length() < str2.length();
} // from short to long

//Notice: Do not use reference, since transform will change the original value
string::size_type FindStrPosByNS(string strBaseSeq, string strSeed, int iQueryPos) //None Sensetive
{
    std::transform(strBaseSeq.begin(), strBaseSeq.end(), strBaseSeq.begin(), ::toupper);
    std::transform(strSeed.begin(), strSeed.end(), strSeed.begin(), ::toupper);
    return strBaseSeq.find(strSeed, iQueryPos);
}

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
    size_t iPos = 0;
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
    char bpUse = bp;//toupper(bp);
    if( bpUse == 'A')
    {
        return 'T';
    }
    else if(bpUse == 'a')
    {
        return 't';
    }
    else if( bpUse == 'T')
    {
        return 'A';
    }
    else if( bpUse == 't')
    {
        return 'a';
    }
    else if( bpUse == 'G')
    {
        return 'C';
    }
    else if(bpUse == 'g')
    {
        return 'c';
    }
    else if(bpUse == 'C')
    {
        return 'G';
    }
    else if(bpUse == 'c')
    {
        return 'g';
    }
    return 'N';
}

string ClsAlgorithm::GetReverseCompelement(string strOrg, bool bRevsCompelement) //We need keep case sensitive
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

string ClsAlgorithm::CreateBamFile(string strFaName, string& strRefSeq,
                                   string& strReads1Path, string& strReads2Path, bool bMapAll)
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
    string strSpcialFileOfRefIndex = strRefPath + ".bwt";
    if(::access(strSpcialFileOfRefIndex.c_str(), 0) != 0) //This means such file DO NOT exsited
    {
        strCmd = "bwa index -a bwtsw " + strRefPath;
        system(strCmd.c_str());
    }

    /*
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
    strCmd = "bwa sampe -f " + strSamPath + " " +
             strRefPath + " " + strSai1Path + " " + strSai2Path + " " +
             strReads1Path + " " + strReads2Path;
    system(strCmd.c_str());*/
    string strSamPath = strRootPath + "Read.sam";
    if(bMapAll) //Use map all
        strCmd = "bwa mem -a " + strRefPath + " " + strReads1Path + " " + strReads2Path + " > " + strSamPath;
    else
        strCmd = "bwa mem " + strRefPath + " " + strReads1Path + " " + strReads2Path + " > " + strSamPath;
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

    //8:Remove the sam file and the bam file
    strCmd = "rm -f " + strSamPath;
    system(strCmd.c_str());
    //b) Delete Bam File
    strCmd = "rm -f " + strBamPath;
    system(strCmd.c_str());

    return strSortedBamPath;
}

/*
 *Notice:
 * 因为在这里我们有很多的scaffold candidate， 所以我们需要把map all 打开，这样能够保证有最多数的相应的read能够被很好的map上去
 *
 */
string ClsAlgorithm::CreateBamFileByMultiScaff(vector<St_ScaffoldUnit>& vScaffoldUnit,
                                               string& strReads1Path, string& strReads2Path,
                                               string strBamFileName, En_ScaffoldDataType enDataType,
                                               bool bMapAll)
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
    //Notice: we do not need to build the index everytime.
    //As a result, we could just create it at the first time by detect some special files
    //-->Detect if need to Build Index, build the index if it is needed
    string strSpcialFileOfRefIndex = strRefPath + ".bwt";
    if(::access(strSpcialFileOfRefIndex.c_str(), 0) != 0) //This means such file DO NOTexsited
    {
        strCmd = "bwa index -a bwtsw " + strRefPath;
        system(strCmd.c_str());
    }
    //<--

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
    if(bMapAll) //Use map all
        strCmd = "bwa mem -a " + strRefPath + " " + strReads1Path + " " + strReads2Path + " > " + strSamPath;
    else
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

string ClsAlgorithm::CreateBamFileForMultiRefPEReads(string& strRefPath, string& strReads1Path,
                                                     string& strReads2Path, string strBamFileName, bool bMapAll)
{
    string strRootPath = ClsAlgorithm::GetInstance().GetHigherFolderPath(get_current_dir_name());
    strRootPath += "TempFile/";
    //2: Build Index for this Reference Fa File
    string strCmd = "";
    string strSpcialFileOfRefIndex = strRefPath + ".bwt";
    if(::access(strSpcialFileOfRefIndex.c_str(), 0) != 0) //This means such file DO NOTexsited
    {
        strCmd = "bwa index -a bwtsw " + strRefPath;
        system(strCmd.c_str());
    }

    //4: Generate SAM file
    string strSamPath = strRootPath + (strBamFileName == "" ? "Read.sam" : (strBamFileName + ".sam"));
    if(bMapAll) //Use map all
        strCmd = "bwa mem -a " + strRefPath + " " + strReads1Path + " " + strReads2Path + " > " + strSamPath;
    else
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

string ClsAlgorithm::CreateBamFileForMultiRefSingleReads(string& strRefPath, string& strReadsPath,
                                                         string strBamFileName, bool bMapAll)
{
    string strRootPath = ClsAlgorithm::GetInstance().GetHigherFolderPath(get_current_dir_name());
    strRootPath += "TempFile/";
    //2: Build Index for this Reference Fa File
    string strCmd = "";
    string strSpcialFileOfRefIndex = strRefPath + ".bwt";
    if(::access(strSpcialFileOfRefIndex.c_str(), 0) != 0) //This means such file DO NOT exsited
    {
        strCmd = "bwa index -a bwtsw " + strRefPath;
        system(strCmd.c_str());
    }

    //4: Generate SAM file
    string strSamPath = strRootPath + (strBamFileName == "" ? "Read.sam" : (strBamFileName + ".sam"));
    if(bMapAll) //Use map all
        strCmd = "bwa mem -a " + strRefPath + " " + strReadsPath + " > " + strSamPath;
    else
        strCmd = "bwa mem " + strRefPath + " " + strReadsPath + " > " + strSamPath;
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

//This is just for DNA char, since we just considered 4 cases: A, T, G, C
// AlignType: 0 means 左对齐，1 means 右对齐
// Core Idea: Collect the most frequency bp from ACGTN -->So the charactor different from ACGTN will not be took into consideration
void ClsAlgorithm::GetMostFreqChar(vector<string>& vOrgStr, vector<char>& vMostFreqChar,
                                   vector<int>& vFreqValue, int iMaxLen,
                                   En_AlignType enAlignType, bool bTight)
{
    vMostFreqChar.clear();
    vFreqValue.clear();
    int iTotalNum = (int)vOrgStr.size();
    for(int i=0; i <iMaxLen; i++)
    {
        char cDNA[vOrgStr.size()];
        int iIndex = 0;
        for(vector<string>::iterator itr = vOrgStr.begin(); itr != vOrgStr.end(); itr++)
        {
            if(enAlignType == atRight) // 如果是右对齐
            {
                if((int)itr->length() > i)
                    cDNA[iIndex] = itr->at(itr->length()-i-1);
                else
                    cDNA[iIndex] = 'N';
            }
            else //否则则默认为0，属于左对齐的case
            {
                if(i < (int)itr->length())
                    cDNA[iIndex] = itr->at(i);
                else
                    cDNA[iIndex] = 'N';
            }
            iIndex++;
        }
        //找到出现频率最高的字符
        int aryCount[4] = {0, 0, 0, 0};
        char aryDNA[4] = {'A', 'T', 'G', 'C'};
        for(int iChar=0; iChar<(int)vOrgStr.size(); iChar++)
        {
            switch(cDNA[iChar])
            {
                case 'A':
                    aryCount[0]++;
                    break;
                case 'T':
                    aryCount[1]++;
                    break;
                case 'G':
                    aryCount[2]++;
                    break;
                case 'C':
                    aryCount[3]++;
                    break;
                case 'N':
                    break;
            }
        }
        int iMaxIndex = 0;
        int iMaxCount = aryCount[0];
        int iValidSum = aryCount[0]; //这个主要是用来过滤的，考虑到有的位置仅仅是某一个(或者少数几个)sequence中是有值的，那么我们统计发现概率很小的话，我们就不去考虑这个值
        for(int i=1; i<4; i++)
        {
            if(iMaxCount < aryCount[i])
            {
                iMaxIndex = i;
                iMaxCount = aryCount[i];
            }
            iValidSum += aryCount[i];
        }
        //这个选项用于控制psuedo repeat中的那种较loss得到combination
        if(bTight)
        {
            //在这里避免低概率的字符被选择-->这里暂时制定的是20%，如果低于20%，那么我们就丢弃,
            //这里是为了解决将少数派的碱基也作为有效碱基插入，导致错误的extend
            if((float)iValidSum / iTotalNum < .2)
                continue;
        }
        //不同的case，插入的方式不一样：右对齐是在头部插入，做对齐是在尾部插入
        if(enAlignType == atRight) //  如果是右对其
        {
            if(iMaxCount != 0) //防止盲目的加入，因为考虑到都为"空"的情况
            {
                vMostFreqChar.insert(vMostFreqChar.begin(), aryDNA[iMaxIndex]);
                vFreqValue.insert(vFreqValue.begin(), iMaxCount);
            }
            else
            {
                vMostFreqChar.insert(vMostFreqChar.begin(), NULL);
                vFreqValue.insert(vFreqValue.begin(), -1);
            }
        }
        else //否则则默认为0，属于左对齐的case
        {
            if(iMaxCount != 0) //防止盲目的加入，因为考虑到都为"空"的情况
            {
                vMostFreqChar.push_back(aryDNA[iMaxIndex]);
                vFreqValue.push_back(iMaxCount);
            }
            else
            {
                vMostFreqChar.push_back(NULL);
                vFreqValue.push_back(-1);
            }
        }
    }
}

/*Check the relationship between two ranges
 * get the overlap
 * case 0: depart
 * case 1: the mapping range be contained by the org gap
 * case 2: org gap range be contained by the mapping range
 * case 3: the later part of mapping range overlapped with org gap
 * case 4: the former part of mapping range overlapped with org gap
 */
St_RltBTRange ClsAlgorithm::GetTwoRangeRelationship(int iStart1, int iEnd1, int iStart2, int iEnd2)
{
    St_RltBTRange stRlt; //the struct of relationship between two ranges
    if(iEnd1 <= iStart1 || iEnd2 <= iStart2)
        return stRlt;
    //case 1: no over lap
    if(iEnd1 <= iStart2 || iStart1 >= iEnd2)
    {
        stRlt.enRangeType = rrDepart;
        return stRlt;
    }
    //Case 2: contained by each other
    if(iStart1 >= iStart2 && iEnd1 <= iEnd2)
    {
        stRlt.iOverlapStart = iStart1;
        stRlt.iOverlapEnd = iEnd1;
        stRlt.enRangeType = rrContain;
        return stRlt;

    }
    if(iStart2 >= iStart1 && iEnd2 <= iEnd1)
    {
        stRlt.iOverlapStart = iStart2;
        stRlt.iOverlapEnd = iEnd2;
        stRlt.enRangeType = rrContain;
        return stRlt;
    }
    //Case 3: overlapped by each other --> both former part overlap and later part overlap
    if(iStart1 < iStart2 && iEnd1 > iStart2 && iEnd1 < iEnd2)
    {
        //Update
        stRlt.iOverlapStart = iStart2;
        stRlt.iOverlapEnd = iEnd1;
        stRlt.enRangeType = rrLeftOverlap;
        return stRlt;
    }
    if(iStart2 < iStart1 && iEnd2 > iStart1 && iEnd2 < iEnd1)
    {
        //Update
        stRlt.iOverlapStart = iStart1;
        stRlt.iOverlapEnd = iEnd2;
        stRlt.enRangeType = rrRightOverlap;
        return stRlt;
    }
}

bool ClsAlgorithm::IsWholeUpperCase(string& strValue) // just check a,t,g,c
{
    if(strValue.find('a') != string::npos ||
       strValue.find('t') != string::npos ||
       strValue.find('g') != string::npos ||
       strValue.find('c') != string::npos)
        return false;
    else
        return true;
}

//Display the string in file --> multi-lines
void ClsAlgorithm::DisplayString(ofstream& ofs, string& strValue, int iLenPerLine)
{
    int iLineNum = ceil((float)strValue.length() / iLenPerLine);
    for(int i=0; i<iLineNum; i++)
    {
        if(i + 1 == iLineNum) // This is the last item
            ofs << strValue.substr(i*iLenPerLine, strValue.length() - i*iLenPerLine) << endl;
        else //This is not the last item
            ofs << strValue.substr(i*iLenPerLine, iLenPerLine) << endl;
    }
}
