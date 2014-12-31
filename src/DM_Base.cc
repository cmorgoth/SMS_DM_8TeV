#include "DM_Base.hh"
#include <iostream>
#include <iomanip>

BaseDM::BaseDM(){
  MRMin = 200.0;//nominal 150.
  RSQMin = 0.5;//nomimal 0.5
  //Does nothing
};


BaseDM::BaseDM( const char* FileName, TString pName, float sig, int MetIndex ){
  
  F = TFile::Open(FileName);
  T = (TTree*)F->Get("outTree");

  this->sigma = sig;//Sets cross sections from derived class namely (DYJets, WJets, TTJets, Data)
  this->processName = pName;
  this->metIndex = MetIndex;
  
  TFile* file = new TFile("/afs/cern.ch/work/c/cpena/scratch_DM/CMSSW_5_2_3/src/VecbosApp/V07_5_3_X/MuHad-Run2012AB/out/hlt_eff_mr200_MoreBin_ABCD_PT80v2Muon.root");

  eff = (TEfficiency*)file->Get("Eff2d");
  TFile* file1 = new TFile("/afs/cern.ch/work/c/cpena/scratch_DM/CMSSW_5_2_3/src/VecbosApp/V07_5_3_X/EleHad-Run2012AB/out/hlt_eff_mr200_MoreBin_PT80.root");

  eff_ele = (TEfficiency*)file->Get("Eff2d");
  
  if( btagIndex == 0 || btagIndex == 1 ){
    this->BtagBranch = "nBtag";
  }else{
    this->BtagBranch = "nBtagTight";
  }
  
  std::cout << "----------------Branch Name:  " << BtagBranch << std::endl;
  
  CalcWeight();
  MRMin = 200.0;//nominal 150.
  RSQMin = 0.5;//nomimal 0.5
};

BaseDM::BaseDM( const char* FileName, TString pName , int MetIndex ){
  
  F = TFile::Open(FileName);
  T = (TTree*)F->Get("outTree");
  
  this->processName = pName;
  this->weight = 1.0;
  this->sigma = -99.;
  this->metIndex = MetIndex;
  
  TFile* file = new TFile("/afs/cern.ch/work/c/cpena/scratch_DM/CMSSW_5_2_3/src/VecbosApp/V07_5_3_X/MuHad-Run2012AB/out/hlt_eff_mr200_MoreBin_ABCD_PT80v2Muon.root");
  eff = (TEfficiency*)file->Get("Eff2d");

  TFile* file1 = new TFile("/afs/cern.ch/work/c/cpena/scratch_DM/CMSSW_5_2_3/src/VecbosApp/V07_5_3_X/EleHad-Run2012AB/out/hlt_eff_mr200_MoreBin_PT80.root");
  eff_ele = (TEfficiency*)file->Get("Eff2d");
   
  if( btagIndex == 0 || btagIndex == 1 ){
    this->BtagBranch = "nBtag";
  }else{
    this->BtagBranch = "nBtagTight";
  }
  
  std::cout << "----------------Branch Name:  " << BtagBranch << std::endl;
    
  MRMin = 200.0;//nominal 150.
  RSQMin = 0.5;//nomimal 0.5
  
};


BaseDM::~BaseDM(){
  //delete T;
};


/*bool BaseDM::FillVecBadHcalEvents(){
  
  ifstream file ;
  file.open("AllBadHCALLaser/AllBadHCALLaser.txt" );
  std::string buffer;
  
  if( file.is_open() ){
    while( file.good() ){
      file >> buffer;
      if( file.eof() )break;
      HcalBadRuns.insert( std::pair<std::string, bool>(buffer, true) );
      //std::cout << "=== buffer: " << buffer << std::endl;
    }
    
  }else{
    std::cout << "==========unable to open bad Hcal runs=======" << std::endl;
  }
  
  file.close();
  
};
*/

bool BaseDM::PrintEvents(){
  std::setprecision(15);
  std::cout << std::fixed;
  
  double NtotGen = 0, Nt_PV = 0, Nt_2J = 0, Nt_0b = 0, Nt_LepVeto = 0, N_In, N_PV, N_2J, N_0b, N_LepVeto;
  TTree* effT = (TTree*)F->Get("effTree");
  effT->SetBranchAddress("Npassed_In", &N_In);
  effT->SetBranchAddress("Npassed_PV", &N_PV);
  effT->SetBranchAddress("Npassed_2Jet", &N_2J);
  effT->SetBranchAddress("Npassed_0btag", &N_0b);
  effT->SetBranchAddress("Npassed_LepVeto", &N_LepVeto);
  
  for (int i = 0; i < effT->GetEntries(); i++){
    effT->GetEntry(i);
    NtotGen += N_In;
    Nt_PV += N_PV;
    Nt_2J += N_2J;
    Nt_0b += N_0b;
    Nt_LepVeto += N_LepVeto;
  }

  std::cout << "========================================" << std::endl;
  std::cout << "============== " << this->processName << " ==================" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "weighted Nt_In: " << NtotGen*weight << std::endl;
  std::cout << "weighted Nt_PV: " << Nt_PV*weight << std::endl;
  std::cout << "weighted Nt_2J: " << Nt_2J*weight << std::endl;
  std::cout << "weighted Nt_0b: " << Nt_0b*weight << std::endl;
  std::cout << "weighted Nt_LepVeto: " << Nt_LepVeto*weight << std::endl;

  std::cout << "NEntries weighted: " << T->GetEntries()*weight << "\n\n" << std::endl; 

  std::cout << "========================================" << std::endl;
  std::cout << "============== " << this->processName << " ==================" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "MC Nt_In: " << NtotGen << std::endl;
  std::cout << "MC Nt_PV: " << Nt_PV << std::endl;
  std::cout << "MC Nt_2J: " << Nt_2J << std::endl;
  std::cout << "MC Nt_0b: " << Nt_0b << std::endl;
  std::cout << "MC Nt_LepVeto: " << Nt_LepVeto << std::endl;

  std::cout << "NEntries weighted: " << T->GetEntries() << "\n\n" << std::endl;

  double RSQ[4], MR[4]/*, run, ls, evNum*/;
  int BOX, nBtag;
  
  SetBrachStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag);
  //T->SetBranchAddress("run", &run);
  //T->SetBranchAddress("ls", &ls);
  //T->SetBranchAddress("evNum", &evNum);
  
  long Nt_MR_RSQ_cut = .0,  Nt_2muBox = .0, Nt_1muBox = .0, Nt_0muBox = .0,	\
    Nt_MR_RSQ_cut0BTag = .0,  Nt_2muBox0BTag = .0, Nt_1muBox0BTag = .0, Nt_0muBox0BTag = .0;
   
  
     
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    /*if( isData ){
      ss << run << ":" << ls << ":" << evNum;
      
      if( HcalBadRuns.find(ss.str()) != HcalBadRuns.end() )std::cout << "======RJE=====" << std::endl;
      
      }*/
    
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin ){
      Nt_MR_RSQ_cut++;
      if( BOX == 0 )Nt_0muBox++;
      if( BOX == 1 )Nt_1muBox++;
      if( BOX == 2 )Nt_2muBox++;
      
      if( nBtag == 0 ){
	Nt_MR_RSQ_cut0BTag++;
	if( BOX == 0 )Nt_0muBox0BTag++;
	if( BOX == 1 )Nt_1muBox0BTag++;
	if( BOX == 2 )Nt_2muBox0BTag++;
	
      }
      
      
    }
  }

  std::cout << "weighted Nt_MR_RSQ_cut: " << Nt_MR_RSQ_cut*weight << std::endl;
  std::cout << "weighted Nt_0muBox: " << Nt_0muBox*weight << std::endl;
  std::cout << "weighted Nt_1muBox: " << Nt_1muBox*weight << std::endl;
  std::cout << "weighted Nt_2muBox: " << Nt_2muBox*weight << std::endl;
  
  std::cout << "weighted Nt_MR_RSQ_cut0BTag: " << Nt_MR_RSQ_cut0BTag*weight << std::endl;
  std::cout << "weighted Nt_0muBox0BTag: " << Nt_0muBox0BTag*weight << std::endl;
  std::cout << "weighted Nt_1muBox0BTag: " << Nt_1muBox0BTag*weight << std::endl;
  std::cout << "weighted Nt_2muBox0BTag: " << Nt_2muBox0BTag*weight << std::endl;
  
  std::cout << this->processName << "Btag EVENTS: " << (Nt_MR_RSQ_cut - Nt_MR_RSQ_cut0BTag)*weight << std::endl;

  std::cout << "========================================" << std::endl;
  std::cout << "========================================" << "\n\n" << std::endl;
  
  return true;
  
};


bool BaseDM::CalcWeight(){
  
  double Ntot = 0, N_In;
  TTree* effT = (TTree*)F->Get("effTree");
  effT->SetBranchAddress("Npassed_In", &N_In);
  for (int i = 0; i < effT->GetEntries(); i++){
    effT->GetEntry(i);
    Ntot += N_In;
  }
  
  weight = Lumi*sigma/Ntot;
  std::cout << "BaseDM weight Process " << this->processName << " is ----> " << weight << std::endl;
  return true;
  
};

TH1F BaseDM::PlotMR_1Box(){
  
  TString BoxCut = "BOX == 1 && ";
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  
  TH1F* MR1 = new TH1F("MR1", "MR1BOX", MR_Bins, MR_BinArr);
  
  SetBrachStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    
    if( BOX==1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex] ){
      if(this->processName == "Data"){
	MR1->Fill(MR[metIndex], weight);
      }else{
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
	if( hltWeight == 0.0 )hltWeight = 1.0;
	MR1->Fill(MR[metIndex], weight*hltWeight);
      }
    }
  }
  
  return *MR1;
};

TH1F BaseDM::PlotMR_0Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  
  TH1F* MR0 = new TH1F("MR0", "MR0BOX", MR_Bins, MR_BinArr);
  
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );

    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if(this->processName == "Data"){
        MR0->Fill(MR[metIndex], weight);
      }else{
	hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	MR0->Fill(MR[metIndex], weight*hltWeight);
      }
    }
  }
  return *MR0;
};


TH1F  BaseDM::PlotRSQ_1Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  
  TH1F* RSQ1 = new TH1F("RSQ1", "RSQ1BOX", RSQ_Bins, RSQ_BinArr);
  
  SetBrachStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    
    if( BOX == 1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if(RSQ[metIndex] > 1.5)std::cout << "out of limits" << std::endl; 
      if(this->processName == "Data"){
        RSQ1->Fill(RSQ[metIndex], weight);
      }else{
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	RSQ1->Fill(RSQ[metIndex], weight*hltWeight);
      }
    }
  }
  
  return *RSQ1;
  
};


TH1F  BaseDM::PlotRSQ_0Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  
  TH1F* RSQ0 = new TH1F("RSQ0", "RSQ0BOX", RSQ_Bins, RSQ_BinArr);
 
  SetBrachStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    
    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if(this->processName == "Data"){
        RSQ0->Fill(RSQ[metIndex], weight);
      }else{
	hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	RSQ0->Fill(RSQ[metIndex], weight*hltWeight);
      }
    }
  }
  
  return *RSQ0;
  
};


TH2F BaseDM::PlotRSQ_vs_MR_0Box(){
  
  double metX[4], metcorrX[4], metY[4], metcorrY[4], ht, RSQ[4], MR[4];
  double mht[3], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  
  TH2F* h = new TH2F("RSQ_MR_0BOX_data", "RSQ_VS_MR_0BOX_data", MR_Bins, MR_BinArr, RSQ_Bins, RSQ_BinArr);
  h->Sumw2();
  
  SetBrachStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("ht", &ht);
  T->SetBranchAddress("mht", mht);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metCorrX", metcorrX);
  T->SetBranchAddress("metY", metY);
  T->SetBranchAddress("metCorrY", metcorrY);
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    
    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if(this->processName == "Data"){
        h->Fill(MR[metIndex], RSQ[metIndex], weight);
      }else{
	hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	h->Fill(MR[metIndex], RSQ[metIndex], weight*hltWeight);
      }
    }
  }
  
  return *h;
};

TH2F BaseDM::PlotRSQ_vs_MR_1Box(){
  
  double metX[4], metcorrX[4], metY[4], metcorrY[4], ht, RSQ[4], MR[4];
  double mht[3], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  
  TH2F* h = new TH2F("RSQ_MR_1BOX_data", "RSQ_VS_MR_1BOX_data", MR_Bins, MR_BinArr, RSQ_Bins, RSQ_BinArr);
  h->Sumw2();
  
  SetBrachStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("ht", &ht);
  T->SetBranchAddress("mht", mht);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metCorrX", metcorrX);
  T->SetBranchAddress("metY", metY);
  T->SetBranchAddress("metCorrY", metcorrY);
  std::cout << "Nentries:  " << T->GetEntries() << std::endl;
  for(long j = 0; j < T->GetEntries(); j++){
    T->GetEntry(j);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    
    if( BOX == 1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if(this->processName == "Data"){
        h->Fill(MR[metIndex], RSQ[metIndex], weight);
      }else{
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	h->Fill(MR[metIndex], RSQ[metIndex], weight*hltWeight);
      }
    }
  }
  
  return *h;
};

TH2F BaseDM::PlotRSQ_vs_MR_2Box(){
  
  double metX[4], metcorrX[4], metY[4], metcorrY[4], ht, RSQ[4], MR[4];
  double mht[3], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  
  TH2F* h = new TH2F("RSQ_MR_2BOX_data", "RSQ_VS_MR_2BOX_data", MR_Bins, MR_BinArr, RSQ_Bins, RSQ_BinArr);
  h->Sumw2();
  
  SetBrachStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("ht", &ht);
  T->SetBranchAddress("mht", mht);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metCorrX", metcorrX);
  T->SetBranchAddress("metY", metY);
  T->SetBranchAddress("metCorrY", metcorrY);
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    
    if( BOX == 2 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex] ){
      if(this->processName == "Data"){
        h->Fill(MR[metIndex], RSQ[metIndex], weight);
      }else{
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	h->Fill(MR[metIndex], RSQ[metIndex], weight*hltWeight);
      }
    }
  }
  
  return *h;
  
};


TH1F  BaseDM::PlotRSQ_2Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  TH1F* RSQ2 = new TH1F("RSQ2_DY", "RSQ2BOX_DY", RSQ_Bins, RSQ_BinArr);
  SetBrachStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    
    if( BOX == 2 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex] ){
      if(this->processName == "Data"){
        RSQ2->Fill(RSQ[metIndex], weight);
      }else{
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	RSQ2->Fill(RSQ[metIndex], weight*hltWeight);
      }
    }
  }
  
  return *RSQ2;
  
};


TH1F BaseDM::PlotMR_2Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  TH1F* MR2 = new TH1F("MR2_DY", "MR2BOX_DY", MR_Bins, MR_BinArr);
  SetBrachStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );

    if( BOX == 2 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex] ){
      if(this->processName == "Data"){
        MR2->Fill(MR[metIndex], weight);
      }else{
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	MR2->Fill(MR[metIndex], weight*hltWeight);
      }
    }
  }
  
  return *MR2;
};

std::vector<TH1F*> BaseDM::PlotMETx(){
  double metX[4], metcorrX[4], metY[4], metcorrY[4], ht, RSQ[4], MR[4]/*, run, ls, evNum*/;
  double mht[3], CSV[30];
  int BOX, N_Jets, nBtag[2];
  
  std::vector< TH1F* > metvec;
  
  TH1F* METx[9];
  TString name;
  for(int l = 0; l < 3; l++ ){
    for( int m = 0; m < 3; m++ ){
      name = TString(Form("BaseDM_metXBox%d_plotType%d",l,m));
      METx[3*l + m] = new TH1F( name, name, 50, -300, 300 );
    }
  }
  
  SetBrachStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("ht", &ht);
  T->SetBranchAddress("mht", &mht[0]);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metCorrX", metcorrX);
  T->SetBranchAddress("metY", metY);
  T->SetBranchAddress("metCorrY", metcorrY);
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if( BOX == 0){
	METx[0]->Fill(metX[metIndex], weight);
	METx[1]->Fill(metcorrX[metIndex], weight);
	METx[2]->Fill(metcorrX[metIndex]-metX[metIndex], weight);
      }else if( BOX == 1 ){
	METx[3]->Fill(metX[metIndex], weight);
	METx[4]->Fill(metcorrX[metIndex], weight);
	METx[5]->Fill(metcorrX[metIndex]-metX[metIndex], weight);
      }else if( BOX == 2 ){
	METx[6]->Fill(metX[metIndex], weight);
	METx[7]->Fill(metcorrX[metIndex], weight);
	METx[8]->Fill(metcorrX[metIndex]-metX[metIndex], weight);
      }
    }
    
  }
  
  for(int j = 0; j < 9; j++){
    metvec.push_back(METx[j]);
  }
    
  return metvec;
};


std::vector<TH1F*> BaseDM::PlotMETy(){
  double metX[4], metcorrX[4], metY[4], metcorrY[4], ht, RSQ[4], MR[4]/*, run, ls, evNum*/;
  double mht[3], CSV[30];
  int BOX, N_Jets, nBtag[2];
  
  std::vector< TH1F* > metvec;
  TH1F* MET[9];
  TString name;
  for(int l = 0; l < 3; l++ ){
    for( int m = 0; m < 3; m++ ){
      name = TString(Form("BaseDM_metY_Box%d_plotType%d",l,m));
      MET[3*l + m] = new TH1F( name, name, 50, -300, 300 );
    }
  }
  
  SetBrachStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("ht", &ht);
  T->SetBranchAddress("mht", &mht[0]);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metCorrX", metcorrX);
  T->SetBranchAddress("metY", metY);
  T->SetBranchAddress("metCorrY", metcorrY);
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if( BOX == 0){
	MET[0]->Fill(metY[metIndex], weight);
	MET[1]->Fill(metcorrY[metIndex], weight);
	MET[2]->Fill(metcorrY[metIndex]-metY[metIndex], weight);
      }else if( BOX == 1 ){
	MET[3]->Fill(metY[metIndex], weight);
	MET[4]->Fill(metcorrY[metIndex], weight);
	MET[5]->Fill(metcorrY[metIndex]-metY[metIndex], weight);
      }else if( BOX == 2 ){
	MET[6]->Fill(metY[metIndex], weight);
	MET[7]->Fill(metcorrY[metIndex], weight);
	MET[8]->Fill(metcorrY[metIndex]-metY[metIndex], weight);
	//std::cout << "METY debug: " << "entry: "  <<  i << " RSQ: " << RSQ << " MR: " << MR \
	//	  << " nBtag:" << nBtag << std::endl;
      }
    }
    
  }
  
  for(int j = 0; j < 9; j++){
    metvec.push_back(MET[j]);
  }
  
  return metvec; 
};


std::vector<TH1F*> BaseDM::PlotMETmag(){
  double metX[4], metcorrX[4], metY[4], metcorrY[4], ht, RSQ[4], MR[4], run/*, ls, evNum*/;
  double mht[3], CSV[30];
  int BOX, nBtag[2], N_Jets;
  
  std::vector< TH1F* > metvec;
  TH1F* MET[12];
  TString name;
  for(int l = 0; l < 3; l++ ){
    for( int m = 0; m < 3; m++ ){
      name = TString(Form("BaseDM_METmag_Box%d_plotType%d",l,m));
      MET[3*l + m] = new TH1F( name, name, 50, 0, 1000);
    }
  }
  
  MET[9] = new TH1F( this->processName + "_NJETS0_Z", "NJETS 0 BOX", 9, 1, 10);
  MET[10] = new TH1F( this->processName + "_NJETS1_Z", "NJETS 1 BOX", 9, 1, 10);
  MET[11] = new TH1F( this->processName + "_NJETS2_Z", "NJETS 2 BOX", 9, 1, 10);
  
  SetBrachStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("ht", &ht);
  T->SetBranchAddress("mht", &mht[0]);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metCorrX", metcorrX);
  T->SetBranchAddress("metY", metY);
  T->SetBranchAddress("metCorrY", metcorrY);
  T->SetBranchAddress("N_Jets", &N_Jets);

  double metmag =0;
  double metmagcorr = 0;
  double wt = 1.;
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );

    metmag = sqrt(metX[metIndex]*metX[metIndex]+metY[metIndex]*metY[metIndex]);
    metmagcorr = sqrt(metcorrX[metIndex]*metcorrX[metIndex]+metcorrY[metIndex]*metcorrY[metIndex]);
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if( BOX == 0){
	if(this->processName == "Data"){
	  wt = 1.0;
	}else{
	  wt = HLTscaleEle(MR[metIndex], RSQ[metIndex]);
	  if( wt == 0.0)wt = 1.0;
	}
	MET[0]->Fill(metmag, weight*wt);
	MET[1]->Fill(metmagcorr, weight*wt);
	MET[2]->Fill(metmagcorr-metmag, weight*wt);
	MET[9]->Fill(N_Jets, weight*wt);
      }else if( BOX == 1 ){
	if(this->processName == "Data"){
          wt = 1.0;
        }else{
          wt = HLTscale(MR[metIndex], RSQ[metIndex]);
          if( wt == 0.0)wt = 1.0;
        }
	MET[3]->Fill(metmag, weight*wt);
	MET[4]->Fill(metmagcorr, weight*wt);
	MET[5]->Fill(metmagcorr-metmag, weight*wt);
	MET[10]->Fill(N_Jets,weight*wt);
      }else if( BOX == 2 ){
	if(this->processName == "Data"){
          wt = 1.0;
        }else{
          wt = HLTscale(MR[metIndex], RSQ[metIndex]);
          if( wt == 0.0)wt = 1.0;
        }
	MET[6]->Fill(metmag, weight*wt);
	MET[7]->Fill(metmagcorr, weight*wt);
	MET[8]->Fill(metmagcorr-metmag, weight*wt);
	MET[11]->Fill(N_Jets,weight*wt);
      }
    }
    
  }
  
  for(int j = 0; j < 12; j++){
    metvec.push_back(MET[j]);
  }
  
  return metvec; 
};

std::vector<TH1F*> BaseDM::PlotHT(){
  double metX[4], metcorrX[4], metY[4], metcorrY[4], ht, RSQ[4], MR[4]/*, run, ls, evNum*/;
  double mht[3], CSV[30];
  int BOX, N_Jets, nBtag[2];
  
  std::vector< TH1F* > metvec;
  
  TH1F* HT[3];
  TString name;
  for(int l = 0; l < 3; l++ ){
    name = TString(Form("BaseDM_HT_Box%d",l));
    HT[l] = new TH1F( name, name, 50, 135, 600 );
  }
  
  SetBrachStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("ht", &ht);
  T->SetBranchAddress("mht", &mht[0]);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metCorrX", metcorrX);
  T->SetBranchAddress("metY", metY);
  T->SetBranchAddress("metCorrY", metcorrY);
  //T->SetBranchAddress("run", &run);
  //T->SetBranchAddress("ls", &ls);
  //T->SetBranchAddress("evNum", &evNum);
  
  //std::cout << "Tree entries: " << T->GetEntries() << std::endl;
    for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    
    if(  RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if( BOX == 0 ){
	HT[0]->Fill(ht, weight);
      }else if( BOX == 1 ){
	HT[1]->Fill(ht, weight);
      }else if( BOX == 2 ){
	HT[2]->Fill(ht, weight);
      }
      
    }
    
  }
  
  for(int j = 0; j < 3; j++){
    metvec.push_back(HT[j]);
  }
  
  return metvec;
};

std::vector<TH2F*> BaseDM::Plot_2DRazor(){
  double RSQ[4], MR[4], CSV[30];
  int BOX, N_Jets, nBtag[2];

  std::vector< TH2F* > Razor2DVec;
  TH2F* Razor2D[3];
  TString name;
  double hltWeight;
  std::string type = (std::string)this->processName;
  std::cout << "=== TYPE: " << type << " ===" << std::endl;
  for(int l = 0; l < 3; l++ ){
    name = TString( Form("Razor2D_%dmu_Box", l) );
    name.Insert(9, this->processName);
    Razor2D[l] = new TH2F( name, name, 50, 200., 1500., 50, 0.5, 1.5);
  }

  SetBrachStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);

  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if( BOX == 0){
	if(this->processName == "Data"){
	  Razor2D[0]->Fill(MR[metIndex], RSQ[metIndex], weight);
	}else{
	  hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
	  if( hltWeight == 0.0 )hltWeight = 1.0;
	  Razor2D[0]->Fill(MR[metIndex], RSQ[metIndex], weight*hltWeight);
	}
      }else if( BOX == 1 ){
        if(this->processName == "Data"){
          Razor2D[1]->Fill(MR[metIndex], RSQ[metIndex], weight);
        }else{
          hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
          if( hltWeight == 0.0 )hltWeight = 1.0;
          Razor2D[1]->Fill(MR[metIndex], RSQ[metIndex], weight*hltWeight);
	}
      }else if( BOX == 2 ){
        if(this->processName == "Data"){
          Razor2D[2]->Fill(MR[metIndex], RSQ[metIndex], weight);
	}else{
          hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
          if( hltWeight == 0.0 )hltWeight = 1.0;
          Razor2D[2]->Fill(MR[metIndex], RSQ[metIndex], weight*hltWeight);
        }
      }
    }
    
  }
  
  for(int j = 0; j < 3; j++){
    Razor2DVec.push_back(Razor2D[j]);
  }
  
  return Razor2DVec;
  
};

bool BaseDM::pfJetPassCSVM(double btagOutput){
  if(btagOutput < 0.679)   return false;
  return true;
};

int BaseDM::pfJetPassCSVM(double* CSVM, int N_Jets){
  int nMBtag = 0;
  for(int i = 0; i < N_Jets; i++)if(CSVM[i] >= 0.679)nMBtag++;
  return nMBtag;
};

bool BaseDM::SetBrachStatus(){
  T->SetBranchStatus("*",0); //disable all branches
  T->SetBranchStatus("RSQ",1);
  T->SetBranchStatus("MR",1);
  T->SetBranchStatus("BOX_NUM",1);
  T->SetBranchStatus("nBtag",1);
  T->SetBranchStatus("nBtagTight",1);
  T->SetBranchStatus("N_Jets",1);
  T->SetBranchStatus("CSV",1);
  T->SetBranchStatus("ht",1);
  T->SetBranchStatus("mht",1);
  T->SetBranchStatus("metX",1);
  T->SetBranchStatus("metY",1);
  T->SetBranchStatus("metCorrX",1);
  T->SetBranchStatus("metCorrY",1);
  T->SetBranchStatus("N_Jets",1);
  //T->SetBranchStatus("run",1);
  //T->SetBranchStatus("ls",1);
  //T->SetBranchStatus("evNum",1);
  
};

double BaseDM::HLTscale(double MR, double R2){
  
  int MRbin = -1;
  int R2bin = -1;
  
  const double R2A[] = {0.35, 0.5, 0.6, 0.8, 1.5};
  const double MRA[] = {200., 300., 400. ,500., 2000.};
  
  int Nbins = 4;
  
  for(int j = 0; j <= Nbins; j++){
    if( R2 > R2A[j]){
      if(R2 < R2A[j + 1]){
        R2bin = j+1;
	break;
      }
    }
  }
  
  for(int j = 0; j <= Nbins; j++){
    if( MR > MRA[j]){
      if(MR < MRA[j + 1]){
        MRbin = j+1;
        break;
      }
    }
  }

  //return hlt->GetBinContent( MRbin, R2bin );
  return eff->GetEfficiency(eff->GetGlobalBin(MRbin, R2bin , 0));
  
};

double BaseDM::HLTscaleEle(double MR, double R2){ 

  int MRbin = -1;
  int R2bin = -1;
  
  const double R2A[] = {0.35, 0.5, 0.6, 0.8, 1.5};
  const double MRA[] = {200., 300., 400. ,500., 2000.};
  
  int Nbins = 4;

  for(int j = 0; j <= Nbins; j++){
    if( R2 > R2A[j]){
      if(R2 < R2A[j + 1]){
        R2bin = j+1;
	break;
      }
    }
  }

  for(int j = 0; j <= Nbins; j++){
    if( MR > MRA[j]){
      if(MR < MRA[j + 1]){
        MRbin = j+1;
        break;
      }
    }
  }
  

  //return hlt_ele->GetBinContent( MRbin, R2bin );
  return eff_ele->GetEfficiency(eff_ele->GetGlobalBin(MRbin, R2bin , 0));

};

/*
double BaseDM::HLTscale(double MR, double R2){
  double MRhigh = 2000.;
  double MRlow = 200.;
  int Nbins = 50;
  std::vector<double> MRlimit;
  double stepMR = (MRhigh - MRlow)/Nbins;

  double R2high = 1.5;
  double R2low = 0.5;
  std::vector<double> R2limit;
  double stepR2 = (R2high - R2low)/Nbins;

  for(int i = 0; i <= Nbins; i++){
    MRlimit.push_back(MRlow+i*stepMR);
    R2limit.push_back(R2low+i*stepR2);
  }

  int binH = Nbins;
  int binL = 0;
  int bin = (int)binH/2;
  int MRbin = -1;

  for(int j = 0; j < Nbins; j++){//Binary Search!!                                                              
    if( MR > MRlimit[bin]){
      if(MR < MRlimit[bin + 1]){
        MRbin = bin+1;
        break;
      }else{
        binL = bin;
        bin = binL + (int)(binH-binL)/2;
        continue;
      }
    }else if(  MR < MRlimit[bin] ){
      if(MR > MRlimit[bin - 1]){
        MRbin = bin;
        break;
      }else{
        binH = bin;
        bin = binL + (int)(binH-binL)/2;
        continue;
      }
    }else{
      std::cout << "number is:  " << MRlimit[bin] << std::endl;
      MRbin = bin+1;
      break;
    }

  }

  binH = Nbins;
  binL = 0;
  bin = (int)binH/2;
  int R2bin = -1;
  for(int j = 0; j < Nbins; j++){//Binary Search!!
    if( R2 > R2limit[bin]){
      if(R2 < R2limit[bin + 1]){
        R2bin = bin+1;
        break;
      }else{
        binL = bin;
        bin = binL + (int)(binH-binL)/2;
        continue;
      }
    }else if(  R2 < R2limit[bin] ){
      if(R2 > R2limit[bin - 1]){
        R2bin = bin;
        break;
      }else{
        binH = bin;
        bin = binL + (int)(binH-binL)/2;
        continue;
      }
    }else{
      std::cout << "number is:  " << MRlimit[bin] << std::endl;
      R2bin = bin+1;
      break;
    }

  }

  //std::cout << MRbin << " " << R2bin << std::endl;
  return hlt->GetBinContent( MRbin, R2bin );

};

double BaseDM::HLTscaleEle(double MR, double R2){
  double MRhigh = 2000.;
  double MRlow = 200.;
  int Nbins = 50;
  std::vector<double> MRlimit;
  double stepMR = (MRhigh - MRlow)/Nbins;

  double R2high = 1.5;
  double R2low = 0.5;
  std::vector<double> R2limit;
  double stepR2 = (R2high - R2low)/Nbins;

  for(int i = 0; i <= Nbins; i++){
    MRlimit.push_back(MRlow+i*stepMR);
    R2limit.push_back(R2low+i*stepR2);
  }

  int binH = Nbins;
  int binL = 0;
  int bin = (int)binH/2;
  int MRbin = -1;

  for(int j = 0; j < Nbins; j++){//Binary Search!!                                                          
    if( MR > MRlimit[bin]){
      if(MR < MRlimit[bin + 1]){
	MRbin = bin+1;
        break;
      }else{
        binL = bin;
        bin = binL + (int)(binH-binL)/2;
        continue;
      }
    }else if(  MR < MRlimit[bin] ){
      if(MR > MRlimit[bin - 1]){
        MRbin = bin;
        break;
      }else{
        binH = bin;
        bin = binL + (int)(binH-binL)/2;
        continue;
      }
    }else{
      std::cout << "number is:  " << MRlimit[bin] << std::endl;
      MRbin = bin+1;
      break;
    }

  }

  binH = Nbins;
  binL = 0;
  bin = (int)binH/2;
  int R2bin = -1;
  for(int j = 0; j < Nbins; j++){//Binary Search!!
    if( R2 > R2limit[bin]){
      if(R2 < R2limit[bin + 1]){
	R2bin = bin+1;
	break;
      }else{
        binL = bin;
        bin = binL + (int)(binH-binL)/2;
        continue;
      }
    }else if(  R2 < R2limit[bin] ){
      if(R2 > R2limit[bin - 1]){
        R2bin = bin;
	break;
      }else{
        binH = bin;
	bin = binL + (int)(binH-binL)/2;
        continue;
      }
    }else{

      std::cout << "number is:  " << MRlimit[bin] << std::endl;
      R2bin = bin+1;
      break;
    }

  }
  
  //std::cout << MRbin << " " << R2bin << std::endl;
  return hlt_ele->GetBinContent( MRbin, R2bin );

};
*/
