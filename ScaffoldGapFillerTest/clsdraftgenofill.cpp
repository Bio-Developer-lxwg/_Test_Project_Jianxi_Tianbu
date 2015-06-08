#include "clsdraftgenofill.h"
#include "clsalgorithm.h"
#include "unistd.h"

ClsDraftGenoFill::ClsDraftGenoFill()
{}

ClsDraftGenoFill::~ClsDraftGenoFill()
{}

void ClsDraftGenoFill::BaseInit()
{}

void ClsDraftGenoFill::Init(vector<St_Fasta>& vFastaData)
{
    BaseInit();
    m_stContigs.Init(vFastaData); //Init Contig File    
}

/*
 * 我们的策略如下:
 * 1:取得Gap两边的部分，然后进行马匹萍到相应的draft geno 中
 * 2:将的到mapping 排序过后的结果，从原始的里面进行查找，找到属于一个pair的
 * 3:如果这样的pair被找到了，那么我们中间的部分就可以直接填补相应的gap
 * 4: 实际上用draft geno 填补最理想的情况就是一锤子买卖
 * *5：实际上我们可能通过mapping 上的边的邻接部分，去进行相应的pair end reads之类的矫正
 */
//这里是存在限制的，因为这种bwa找不到少部分的match，对于reads的最小单元，它的要求需要时30bps以上才能进行相应的识别，这个跟我们的初衷是向违背的，
//我们并不需要那么多的去匹配(overlap)！！！！！！
void ClsDraftGenoFill::FillScafUnit(St_ScaffoldUnit& stScafUnit)
{
    CalcByBWA(stScafUnit); //this is not work, since the overlap range could not reach to 30 bps
}

void ClsDraftGenoFill::CalcByBWA(St_ScaffoldUnit& stScafUnit)
{
    unsigned int iCutLen = 30; //基本上30是最小能够使得bwa work的单元了
    vector<St_FlankPair> vFlankReads;
    vFlankReads.resize(stScafUnit.vGap.size());
    //Step 1: Collect all of pair flank
    unsigned int iIndex = 0;
    for(vector<St_Gap>::iterator itr = stScafUnit.vGap.begin(); itr != stScafUnit.vGap.end(); itr++)
    {
        //For left flank 这个时候应该是取后面的X个
        if(itr->strRefinedLeftFlank.length() > iCutLen)
            vFlankReads[iIndex].strLeft = itr->strRefinedLeftFlank.substr(itr->strRefinedLeftFlank.length() - iCutLen,
                                                                          iCutLen);
        else
            vFlankReads[iIndex].strLeft = itr->strRefinedLeftFlank;
        //For right flank 直接取前面的X个
        if(itr->strRefinedRightFlank.length() > iCutLen)
            //vFlankReads[iIndex].strRight = ClsAlgorithm::GetInstance().GetReverseCompelement(itr->strRefinedRightFlank.substr(0, iCutLen));
            vFlankReads[iIndex].strRight = itr->strRefinedRightFlank.substr(0, iCutLen);
        else
            //vFlankReads[iIndex].strRight = ClsAlgorithm::GetInstance().GetReverseCompelement(itr->strRefinedRightFlank);
            vFlankReads[iIndex].strRight = itr->strRefinedRightFlank;
        //For its Name
        vFlankReads[iIndex].strName = stScafUnit.strName + "&" +IntToStr(iIndex) + "#";
        iIndex++;
    }
    //相当于生成一个pair end reads
    //1: 首先生成相应的fa文件
    string strRootPath = ClsAlgorithm::GetInstance().GetHigherFolderPath(get_current_dir_name());
    string strTempFileFolderPath = strRootPath + "TempFile/";
    //Clear the files under this folder
    string strCmd = "";
    //Set the valuefor Reference Fasta file
    string strFa1Name = "FlankReads1";
    string strFa2Name = "FlankReads2";
    string strReads1Path = strTempFileFolderPath + strFa1Name + ".fa";
    string strFastq1Path = strTempFileFolderPath + strFa1Name + ".fq";
    string strReads2Path = strTempFileFolderPath + strFa2Name + ".fa";
    string strFastq2Path = strTempFileFolderPath + strFa2Name + ".fq";

    ofstream ofs1;
    ofs1.open(strReads1Path.c_str());
    ofstream ofs2;
    ofs2.open(strReads2Path.c_str());
    //Name
    for(vector<St_FlankPair>::iterator itr = vFlankReads.begin(); itr != vFlankReads.end(); itr++)
    {
        //Name
        ofs1 << ">" << itr->strName << "/1" << endl;
        ofs2 << ">" << itr->strName << "/2" << endl;
        //Value
        ofs1 << itr->strLeft << endl; //Reads 1
        ofs2 << itr->strRight << endl; //Reads 2
    }
    ofs1.close();
    ofs2.close();

    //2: 然后将fa文件转换成fq文件
    strCmd = "perl " + strRootPath + "/ThirdPartyTools/fasta_to_fastq.pl " +
             strReads1Path + " > " + strFastq1Path;
    system(strCmd.c_str());
    strCmd = "perl " + strRootPath + "/ThirdPartyTools/fasta_to_fastq.pl " +
             strReads2Path + " > " + strFastq2Path;
    system(strCmd.c_str());

    //3: 将所有的draft geno 合并成一个 fasta文件用于之后的bwa
    string strDraftFa = "DraftContigs";
    string strDraftPath = strTempFileFolderPath + strDraftFa + ".fa"; //This could be viewed as reference Geno
    ofstream ofs;
    ofs.open(strDraftPath.c_str());
    for(vector<St_ContigUnit>::iterator itr = m_stContigs.vData.begin();
        itr != m_stContigs.vData.end(); itr++)
    {
        //Name
        ofs << ">" << itr->strName << endl;
        //Value
        ofs << itr->strSeq << endl;
    }
    ofs.close();

    //string strBamFilePath = ClsAlgorithm::GetInstance().CreateBamFileForMultiRefSingleReads(strDraftPath,
    //                                                                                    strFastq1Path,
    //                                                                                    "DraftGeno");
    //
    //4: Align the PE Pseudo Reads Back to Pseudo Reference
    string strBamFilePath = ClsAlgorithm::GetInstance().CreateBamFileForMultiRefPEReads(strDraftPath,
                                                                                        strFastq1Path,
                                                                                        strFastq2Path,
                                                                                        "DraftGeno");
    //5: Parse the bam file to find our if there are something could be used for filling directely (the seq between overlap)
    BamReader* pBamReader = new BamReader();
    pBamReader->Open(strBamFilePath);
    pBamReader->OpenIndex(strBamFilePath + ".bai");
    BamAlignment al;
    string strFillSeq = "";
    while(pBamReader->GetNextAlignment(al))  //现在最简单的方式：找到一个就可以了--->
    {
        if(al.IsMapped())
            cout << "GoodA" << endl;
        if(al.IsMateMapped())
            cout << "GoodB" << endl;
        // "al.AlignedBases != " proves that current reads mapped
        // "al.IsMateMapped()" means its mate matches fine
        // "al.RefID == al.MateRefID" means they mapped into the same draft geno
        if(al.IsMapped() && al.IsMateMapped() &&
           al.RefID == al.MateRefID)
        {
            //找到他是属于哪一个Gap
            int iStartPos = al.Filename.find('&') + 1;
            int iLen = al.Filename.find('#') - al.Filename.find('&') - 1;
            int iGapIndex = atoi(al.Filename.substr(iStartPos, iLen).c_str());

            //然后check 这个gap是不是已经填充完毕了
            if(stScafUnit.vGap[iGapIndex].bFilled)
                continue;

            //若没有填充那么进行相应的填充
            int alPos = al.GetEndPosition();
            int alMatPos = al.MatePosition;
            iLen = abs(alPos - alMatPos);
            iStartPos = alPos > alMatPos ? alMatPos :alPos;
            strFillSeq = m_stContigs.vData[al.RefID].strSeq.substr(iStartPos, iLen);
            stScafUnit.vGap[iGapIndex].strFillSeq = strFillSeq;
            stScafUnit.vGap[iGapIndex].bFilled = true;
        }
    }
    delete pBamReader;
    pBamReader = NULL;
}




