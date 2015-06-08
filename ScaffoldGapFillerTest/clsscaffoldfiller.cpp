#include "clsscaffoldfiller.h"
#include "local_alignment.h"
#include "smith-waterman.h"
#include "stdlib.h"
#include <unistd.h>
#include "clsalgorithm.h"
#include "algorithm"
#include "math.h"

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#define __Testing

ClsScaffoldFiller::ClsScaffoldFiller(): m_pFastaParse(NULL)
{
    m_pClPseudoRpt = NULL;
    m_pClDraftGeno = NULL;
    m_pClPEReads = NULL;
    m_pQualityCtr = NULL;
}

ClsScaffoldFiller::~ClsScaffoldFiller()
{
    //Release m_pFastaParse
    if(m_pFastaParse)
    {
        delete m_pFastaParse;
        m_pFastaParse = NULL;
    }
    //Release m_pClPseudoRpt
    if(m_pClPseudoRpt)
    {
        delete m_pClPseudoRpt;
        m_pClPseudoRpt = NULL;
    }
    //Release m_pClDraftGeno
    if(m_pClDraftGeno)
    {
        delete m_pClDraftGeno;
        m_pClDraftGeno = NULL;
    }
    //Release m_pClPEReads
    if(m_pClPEReads)
    {
        delete m_pClPEReads;
        m_pClPEReads = NULL;
    }
    //Release m_pQualityCtr
    if(m_pQualityCtr)
    {
        delete m_pQualityCtr;
        m_pQualityCtr = NULL;
    }
}

//Parameter 1: Repeat File
//Parameter 2: scaffoldFile
//Parameter 3: reads1
//Parameter 4: reads2
void ClsScaffoldFiller::Init(char **argv)
{
    //m_pFastaParse
    if(m_pFastaParse == NULL)
    {
        m_pFastaParse = new FastaParser();
    }
    //m_pClPseudoRpt
    if(m_pClPseudoRpt == NULL) //Repeat File
    {
        m_pClPseudoRpt = new ClsPseudoRepeatFill();
        m_pClPseudoRpt->Init(argv[1], argv[2], argv[3], argv[4]);
    }
    //m_pClDraftGeno
    if(m_pClDraftGeno == NULL)
    {
        m_pClDraftGeno = new ClsDraftGenoFill();        
    }
    //m_pClPEReads
    if(m_pClPEReads == NULL)
    {
        m_pClPEReads = new ClsPairEndFill();
        m_pClPEReads->Init(argv[3], argv[4]);
    }
    //m_pQualityCtr
    if(m_pQualityCtr == NULL)
    {
        m_pQualityCtr = new ClsQualityControl();
        m_pQualityCtr->Init(argv[3], argv[4]);
    }
    int iMaxLen = m_pClPseudoRpt->GetMaxRepeatLen();
    //Init DraftGeno, including Scaffold and Contigs
    //这里我们还需要一个枚举型来定义如何解析相应的contig和scaffold ---->
    //vector<St_Fasta> vFastaData;
    m_vDraftGenoFaData.clear();
    m_pFastaParse->ReadFasta(argv[2], m_vDraftGenoFaData);
    m_stScaffold.Init(m_vDraftGenoFaData, iMaxLen); //Init Scaffold File
    m_pClDraftGeno->Init(m_vDraftGenoFaData); //Init Contig File
    //<----
    m_strReads1Path = argv[3]; //Set Reads1
    m_strReads2Path = argv[4]; //Set Reads2
    cout << "m_stScaffold Init finished" << endl;
    //手动清理容器，用于释放内存
    //vFastaData.clear();
}

void ClsScaffoldFiller::FillScaffold()
{
    if(m_stScaffold.vData.empty())
        return;
    m_vFinalResult.clear();
    //move the file -->Which record the result of gap filling
    string strResultPath = ClsAlgorithm::GetInstance().GetCurExeFolderPath() + "/GapFillerInfo";
    ::remove(strResultPath.c_str());
    //move all of files in temperary folder
    string strCmd = "exec rm -r " +
                    ClsAlgorithm::GetInstance().GetHigherFolderPath(ClsAlgorithm::GetInstance().GetCurExeFolderPath()) +
                    "TempFile/*";
    system(strCmd.c_str());    
    //第一步，首先进行GapLength的评估---->这个是受到了那个paper的启示
    //m_pQualityCtr->EstimateGapLength(m_stScaffold, m_strReads1Path, m_strReads2Path);
    //将结果输出来
#ifdef __Testing
    int iScaffoldIndex = 1;
    for( vector<St_ScaffoldUnit>::iterator itr = m_stScaffold.vData.begin();
         itr != m_stScaffold.vData.end(); itr++)
    {
        int iGapIndex = 0;
        for(vector<St_Gap>::iterator itrGap = itr->vGap.begin(); itrGap != itr->vGap.end(); itrGap++)
        {
            if(itrGap->bAbnormalLarge)
            {
                cout << "Scaffold[" << IntToStr(iScaffoldIndex) << "]: " << endl;
                cout <<"______Gap[" << IntToStr(iGapIndex) << "] Abmormal: " << "Yes" << endl;
            }
            iGapIndex++;
        }
        iScaffoldIndex++;
    }
#endif

    //FIlling Notice: -->
    //1: We do not try to fill every gap!!! The gap with a small number of N will not be considered to be repared!!!!!
    //<--
    //在此处我想在服务器上做一个实验，也就是

    //然后我们需要做的是对剩下的相应的scaffold进行repeat的填补
    //Step 2: Fill the gap for each scaffold unit
    //整体逻辑应该是
    //(1)一次性先把scaffold 用repeat填补完    
    //for(vector<St_ScaffoldUnit>::iterator itr = m_stScaffold.vData.begin();
    //    itr != m_stScaffold.vData.end(); itr++)
    for(vector<St_ScaffoldUnit>::iterator itr = m_stScaffold.vData.end() - 1;
        itr >= m_stScaffold.vData.begin(); itr--)
    {
        if(itr->vGap.empty()) //Do not need gap filling
            continue;
        //Step 1: try to use repeat to fill the gap
        //Here we use the pseudo repeat to fill the gap by using the Draft Geno!!!!!!!!  -->Works fine!!!!
        m_pClPseudoRpt->FillScafUnitWithSet(*itr, m_vDraftGenoFaData, m_stScaffold.vData);

        //Step 2: 用其他scaffold里面的去填补完成
        //在这里把使用draft geno填补的逻辑给补充完成 ----> Go!!!!!!!!!!
        //m_pClDraftGeno->FillScafUnit(*itr);

        //Step 3: Use Pair End Reads to fill the gap --> In this case we fill the gap one by one rather than build just one BAM File for all of those gaps
        m_pClPEReads->FillScafUnit(*itr);

        //Step 4: 用pair-end reads去进行相应的结果验证
        //(3)然后用pair-end reads第三次的填补: 分开的原因是因为，这一步需要作bwa，
        //   而input的数据是需要通过上面两个算法进行fixed后的数据
        //m_pClPEReads->FillScafSet(m_stScaffold.vData); //---> do bwa everytime for each scaffold
        //(4)最后再是用pair-end reads去进行相应的结果验证
        itr->CheckFillingResult(this->m_strReads1Path, this->m_strReads2Path);
    }

    //(5) 将填充结果展现出来
    ofstream ofs;
    string strRootPath = ClsAlgorithm::GetInstance().GetCurExeFolderPath();
    ofs.open((strRootPath + "/GapFillerInfo").c_str()); //Do not need to use the "a+", since we write such file just one time

    //---> We need to write file
    for(vector<St_ScaffoldUnit>::iterator itr = m_stScaffold.vData.begin();
        itr != m_stScaffold.vData.end(); itr++)
    {
        ofs << itr->strName << endl;
        ofs << "Filling Seq: " << endl;
        ClsAlgorithm::GetInstance().DisplayString(ofs, itr->stFinalRst.strSeq);

        int iIndex = 0;
        for(vector<int>::iterator itrCoverage = itr->stFinalRst.vCoverage.begin();
            itrCoverage != itr->stFinalRst.vCoverage.end(); itrCoverage++, iIndex++)
        {
            ofs << "Gap[" << IntToStr(iIndex) << "]: " << IntToStr(*itrCoverage) << endl;
            ofs << "Gap[" << IntToStr(iIndex) << "] " << "Fill By Repeats: " << endl;
            ClsAlgorithm::GetInstance().DisplayString(ofs, itr->vGapRefinedByRept[iIndex].pOrgGap->strFillSeq);
            ofs << "Gap[" << IntToStr(iIndex) << "] " << "Fill By PE Reads: " << endl;
            ClsAlgorithm::GetInstance().DisplayString(ofs, itr->vGapRefinedByRept[iIndex].stBamGap.strFillSeq);
        }
        ofs << "******************" << endl;
    }
    //<--

    /*
    for(vector<St_ScaffoldUnit>::iterator itr = m_stScaffold.vData.begin();
        itr != m_stScaffold.vData.end(); itr++)
    {
        cout << ">" << itr->strName << endl;
        int i = 0;
        for(vector<St_GapRefinedByRept>::iterator subItr = itr->vGapRefinedByRept.begin();
            subItr != itr->vGapRefinedByRept.end(); subItr++)
        {
            cout << "Gap Index [" << IntToStr(i) << "]:" << endl;
            cout << "Fill by Repeats ==>" << subItr->pOrgGap->strFillSeq << endl;
            cout << "Final_PE: Fill by PE Reads: Norm Ore ==> " << subItr->stBamGap.strFillSeq << endl;
            cout << "Final_PE: Fill by PE Reads: Reverse Complementary ==> " << subItr->stBamGap.strRevCompFillSeq << endl;
            cout << "Temp_PE: left Extend: " << subItr->stBamGap.strLeftExt << endl;
            cout << "Temp_PE: right Extend: " << subItr->stBamGap.strRightExt << endl;
            cout << "=========================================" << endl;
            i++;
        }
    }*/
    ofs.close();
}

/*
void ClsScaffoldFiller::CheckFillingResult() // --->Go!!!!!!!!!!!!!
{    
    for(vector<St_ScaffoldUnit>::iterator itr = m_stScaffold.vData.begin();
        itr != m_stScaffold.vData.end(); itr++)
    {        
        //Step1: Arrange the filling information and generate the real filling result
        itr->UpdateFinalFillResult(); //Finished
        //Step2: Create fasta file & Step3: Use pair-end reads map back to fasta file
        string strBamFilePath = ClsAlgorithm::GetInstance().CreateBamFile(itr->strName, itr->strRefinedSeq,
                                                                          m_strReads1Path, m_strReads2Path);
        //Step4: Analysis the result -->how to get it
        itr->stFinalRst.QualityCheckByPEReads(strBamFilePath);

        //Step5: If the result is ok, refine the result again by considering all the reads
        //which maps back to the rang of gap
    }
}*/

/*
void ClsScaffoldFiller::FillScaffoldUnit(St_ScaffoldUnit& stScaffoldUnit)
{
    //Step 1: try to use repeat to fill the gap
    //FillScaffoldUnitByRepeats(stScaffoldUnit);
    //Step 2: Use Pair end Reads to fix the gap:
    //FillScaffoldUnitByPairEndReads(stScaffoldUnit);

    //--------->Just for temperay output: we try to print out those result
    //In this phase: we just want to list all the gaps which could be filled by repeats and pair-end reads
    ofstream ofs;
    string strRootPath = ClsAlgorithm::GetInstance().GetCurExeFolderPath();
    ofs.open((strRootPath + "/GapFillerInfo").c_str(), _S_app); //Do not need to use the "a+", since we write such file just one time
    {
        ofs << ">" << stScaffoldUnit.strName << endl;
        int i = 0;
        for(vector<St_GapRefinedByRept>::iterator itr = stScaffoldUnit.vGapRefinedByRept.begin();
            itr != stScaffoldUnit.vGapRefinedByRept.end(); itr++)
        {
            ofs << "Gap Index [" << IntToStr(i) << "]:" << endl;
            ofs << "Fill by Repeats ==>" << itr->pOrgGap->strFillSeq << endl;
            ofs << "Final_PE: Fill by PE Reads: Norm Ore ==> " << itr->stBamGap.strFillSeq << endl;
            ofs << "Final_PE: Fill by PE Reads: Reverse Complementary ==> " << itr->stBamGap.strRevCompFillSeq << endl;
            ofs << "Temp_PE: left Extend: " << itr->stBamGap.strLeftExt << endl;
            ofs << "Temp_PE: right Extend: " << itr->stBamGap.strRightExt << endl;
            ofs << "=========================================" << endl;
            i++;
        }
    }
    //<----------------

    //Step 3: Check the Quality of fill result    
    //CheckFillQuality(stScaffoldUnit);
}*/

