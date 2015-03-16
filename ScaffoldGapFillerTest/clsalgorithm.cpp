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

#define MAXBUFSIZE 260

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
    while((iPos = (int)strCurPath.find('/', iPos)) != string::npos)
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


