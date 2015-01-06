/*
Generates the input for combine-tool limit setting code
including data card and the root file with the histograms
for the systematics, signal is scale to the correct LUMI
and includes the ISR, PDF, JES systematic, as well as the
PU re-weighting, and the HLT weight for the turn-on curve
*/
//C++ INCLUDES
#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
//ROOT INCLUDES
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
//LOCAL INCLUDES
#include "hlt.hh"
#include "HelperFunctions.hh"
#include "StructDefinition.hh"

int main(){

  //gROOT->Reset();
  
  //Defining Binning for RUN1/8TeV
  const int r2B[4] = {11, 6, 6, 4};
  float c1B[] = {0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.9, 0.95, 1.0, 1.2};
  float c2B[] = {0.50, 0.575, 0.65, 0.75, 0.85, .950, 1.2};
  float c3B[] = {0.50, 0.575, 0.65, 0.75, 0.85, .950, 1.2};
  float c4B[] = {0.50, 0.60, 0.70, .950, 1.2};
  std::vector<float*> v;
  v.push_back(c1B);
  v.push_back(c2B);
  v.push_back(c3B);
  v.push_back(c4B);
  
  TCanvas* ca = new TCanvas("c","c", 640, 640);
  
  //Obtaining Trigger Efficiency File
  TFile* f2 = new TFile("~/Software/git/DM_Signal/trigger/hlt_eff_HTMHT_0mu_BOX_0bTag_Final.root");
  TEfficiency* hlt = (TEfficiency*)f2->Get("Eff2d");

  TH1F* dy[4];
  TH1F* z[4];
  TH1F* w[4];
  TH1F* tt[4];
  TH1F* bkg[4];
  TH1F* data[4];
  //Getting MR categories plots from Prediction
  TString dys, zs, ws, tts, bkgs, datas;
  double tt_N[4];//Total contribution #
  double dy_N[4];//Total contribution #
  double z_N[4];//Total contribution #
  double w_N[4];//Total contribution #
  double data_N[4];//Total contribution #
  //TFile* in = new TFile("/Users/cmorgoth/Software/git/BkgPredictionDM/PredFilesAN/MR_Cat_PredV2_NEW_kF.root");
  TFile* in = new TFile("/Users/cmorgoth/Software/git/BkgPredictionDM/PredFilesAN/MR_Cat_PredV2_NEW_kF_OriginalANA.root");
  for(int i = 0; i < 4; i++){
    dys = TString(Form("cat%d_dy_Pred",i+1));
    zs = TString(Form("cat%d_z_Pred",i+1));
    ws = TString(Form("cat%d_w_Pred",i+1));
    tts = TString(Form("cat%d_tt_Pred",i+1));
    bkgs = TString(Form("cat%d_1D_0mu_Box_Pred_sys",i+1));
    datas = TString(Form("Data_cat%d_1D_0mu_Box",i+1));
    dy[i] = (TH1F*)in->Get(dys);
    z[i] = (TH1F*)in->Get(zs);
    w[i] = (TH1F*)in->Get(ws);
    tt[i] = (TH1F*)in->Get(tts);
    bkg[i] = (TH1F*)in->Get(bkgs);
    data[i] = (TH1F*)in->Get(datas);
    dy_N[i] = dy[i]->Integral();
    z_N[i] = z[i]->Integral();
    w_N[i] = w[i]->Integral();
    tt_N[i] = tt[i]->Integral();
    data_N[i] = data[i]->Integral();
  }
  
  //Creating Histos for Shape Analysis
  TH1F* dy_up[4];
  TH1F* z_up[4];
  TH1F* w_up[4];
  TH1F* tt_up[4];
  
  TH1F* dy_down[4];
  TH1F* z_down[4];
  TH1F* w_down[4];
  TH1F* tt_down[4];
  
  TString bkgn;
  for(int i = 0; i < 4; i++){
    bkgn = TString(Form("dy_cat%d",i));
    dy_up[i] = new TH1F(bkgn, bkgn, r2B[i], v.at(i));
    bkgn = TString(Form("dy_down_cat%d",i));
    dy_down[i] = new TH1F(bkgn, bkgn, r2B[i], v.at(i));
    
    bkgn = TString(Form("z_up_cat%d",i));
    z_up[i] = new TH1F(bkgn, bkgn, r2B[i], v.at(i));
    bkgn = TString(Form("z_down_cat%d",i));
    z_down[i] = new TH1F(bkgn, bkgn, r2B[i], v.at(i));

    bkgn = TString(Form("w_up_cat%d",i));
    w_up[i] = new TH1F(bkgn, bkgn, r2B[i], v.at(i));
    bkgn = TString(Form("w_down_cat%d",i));
    w_down[i] = new TH1F(bkgn, bkgn, r2B[i], v.at(i));

    bkgn = TString(Form("tt_up_cat%d",i));
    tt_up[i] = new TH1F(bkgn, bkgn, r2B[i], v.at(i));
    bkgn = TString(Form("tt_down_cat%d",i));
    tt_down[i] = new TH1F(bkgn, bkgn, r2B[i], v.at(i));

    double min_norm = 0.0001;
    for(int j = 1; j <= dy_up[i]->GetNbinsX(); j++){
      dy_up[i]->SetBinContent(j,dy[i]->GetBinContent(j)+dy[i]->GetBinError(j));
      if(dy[i]->GetBinContent(j)-dy[i]->GetBinError(j) > 0.0){
	dy_down[i]->SetBinContent(j,dy[i]->GetBinContent(j)-dy[i]->GetBinError(j));
      }else{
	dy_down[i]->SetBinContent(j,min_norm);
      }
      
      z_up[i]->SetBinContent(j,z[i]->GetBinContent(j)+z[i]->GetBinError(j));
      if(z[i]->GetBinContent(j)-z[i]->GetBinError(j) > 0.0){
	z_down[i]->SetBinContent(j,z[i]->GetBinContent(j)-z[i]->GetBinError(j));
      }else{
	z_down[i]->SetBinContent(j,min_norm);
      }

      w_up[i]->SetBinContent(j,w[i]->GetBinContent(j)+w[i]->GetBinError(j));
      if(w[i]->GetBinContent(j)-w[i]->GetBinError(j) > 0.0){
	w_down[i]->SetBinContent(j,w[i]->GetBinContent(j)-w[i]->GetBinError(j));
      }else{
	w_down[i]->SetBinContent(j,min_norm);
      }
      
      tt_up[i]->SetBinContent(j,tt[i]->GetBinContent(j)+tt[i]->GetBinError(j));
      if(tt[i]->GetBinContent(j)-tt[i]->GetBinError(j) > 0.0){
	tt_down[i]->SetBinContent(j,tt[i]->GetBinContent(j)-tt[i]->GetBinError(j));
      }else{
	tt_down[i]->SetBinContent(j,min_norm);
      }
      
    }
    
  }

  //Obtaining Input File
  TFile* fin = new TFile("SMS-MadGraph_2J_T2cc_NoFilter_mStop-100to250_mLSP-20to230_8TeV.root");
  TTree* outT = (TTree*)fin->Get("outTree");
  float mssm[3];
  outT->SetBranchStatus("*", 0);
  outT->SetBranchStatus("mssm", 1);
  outT->SetBranchAddress("mssm", mssm);
  
  std::map<std::string, int> sms_aux;
  long nInitEntries = outT->GetEntries();
  for(long i = 0; i < nInitEntries; i++){
    outT->GetEntry(i);
    std::stringstream ss_aux;
    ss_aux << mssm[0] << " " << mssm[1];
    std::string s_aux = ss_aux.str();
    if(sms_aux.find(s_aux) == sms_aux.end()){
      sms_aux[s_aux] = 1;
      std::cout << "MASS PAIR: " << s_aux << std::endl;
    }
  }
 
  std::cout << "MAP SIZE: " << sms_aux.size() << std::endl;
  std::map<std::string, SignalStruct> sms_mass_pair;
  
  //Getting Number of Gen Events Histos
  TH2F* h_Nisr = (TH2F*)fin->Get("h_Nisr");
  TH2F* h_Nisr_up = (TH2F*)fin->Get("h_Nisr_up");
  TH2F* h_Nisr_down = (TH2F*)fin->Get("h_Nisr_down");
  
  //Define Liminosity and xsec
  double Lumi = 18.836;//(fb-1)Jan22Rereco
  double xsec = 30.0*1000.0;//xsec = 20 pb
  
  TFile* facc;
  TFile* f_xsec = new TFile("SMS_LimitsFile/stop.root");
  TH1F* stop = (TH1F*)f_xsec->Get("stop");
  for( auto tmp : sms_aux ){
    std::stringstream tmp_ss;
    tmp_ss.str(tmp.first);
    int mass1, mass2;
    tmp_ss >> mass1 >> mass2;
    std::stringstream s_fname;
    s_fname << mass1 << "_" << mass2;
    SignalStruct aux_struct;
    //xsec = (stop->GetBinContent(stop->FindBin((double)mass1))*1000.0)/50.0;//xsec fb
    std::cout << "Mass" << mass1 << " == xsec: " << xsec/1000.0 << " pb"<< std::endl;
    //if(mass1 >= 200.0 )xsec = 33.0*1000.0;//xsec = 33 pb
    aux_struct.Npassed_ISR = GetGenEvents(mass1, mass2, h_Nisr);
    aux_struct.scaleF = Lumi*xsec/aux_struct.Npassed_ISR;
    std::cout << "scaleF: " << aux_struct.scaleF << std::endl;  
    
    aux_struct.Npassed_ISR_up = GetGenEvents(mass1, mass2, h_Nisr_up);
    aux_struct.scaleF_up = Lumi*xsec/aux_struct.Npassed_ISR_up;
    
    aux_struct.Npassed_ISR_down = GetGenEvents(mass1, mass2, h_Nisr_down);
    aux_struct.scaleF_down = Lumi*xsec/aux_struct.Npassed_ISR_down;
    
    //Create Struct Histos
    CreateStructHistos(aux_struct, tmp.first);//Not using string argument at the moment
    //Getting PDF acceptance Histo
    TString sfile = "AccBIS/SMS_";
    sfile = sfile + s_fname.str().c_str() + ".root";
    std::cout << sfile << std::endl;
    facc = new TFile(sfile);
    for(int i = 0; i < 4; i++){
      TString shisto = TString(Form("PDF_SYS_cat%d",i+1));
      shisto = s_fname.str() + "_" + shisto;
      aux_struct.pdf_acc[i] = new TH1F(*((TH1F*)facc->Get(shisto)));
      //std::cout << aux_struct.pdf_acc[i]->GetBinContent(1) << std::endl;
    }
    std::cout << aux_struct.pdf_acc[0]->GetBinContent(1) << std::endl;
    //Inserting New Struct
    sms_mass_pair[tmp.first] = aux_struct;
    //delete facc;
  }
  
  //Here the program starts
  
  std::cout.precision(16);  
  TFile* f_acc;
  
  double mr[4], rsq[4], Jet_PT[20], Jet_Eta[20], Jet_Phi[20], metCorrX[4], metCorrY[4];
  double mr_up[4], rsq_up[4], Jet_PT_up[20], Jet_Eta_up[20], Jet_Phi_up[20], metCorrX_up[4], metCorrY_up[4];
  double mr_down[4], rsq_down[4], Jet_PT_down[20], Jet_Eta_down[20], Jet_Phi_down[20], metCorrX_down[4], metCorrY_down[4];
  double metX[4], metY[4];
      
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2;
  double pTHem1_up, pTHem2_up, etaHem1_up, etaHem2_up, phiHem1_up, phiHem2_up;
  double pTHem1_down, pTHem2_down, etaHem1_down, etaHem2_down, phiHem1_down, phiHem2_down;
  
  int btag, box, N_Jets;
  double pu_w, ISR, ISR_up, ISR_down;
  
  outT->SetBranchStatus("*", 0);
      
  outT->SetBranchStatus("pu_w", 1);
  outT->SetBranchStatus("ISR", 1);
  outT->SetBranchStatus("ISR_up", 1);
  outT->SetBranchStatus("ISR_down", 1);
      
  outT->SetBranchStatus("MR", 1);
  outT->SetBranchStatus("MR_up", 1);
  outT->SetBranchStatus("MR_down", 1);

  outT->SetBranchStatus("RSQ",1);
  outT->SetBranchStatus("RSQ_up",1);
  outT->SetBranchStatus("RSQ_down",1);
      
  outT->SetBranchStatus("nBtag", 1);
  outT->SetBranchStatus("BOX_NUM",1);
  outT->SetBranchStatus("N_Jets",1);
  
  outT->SetBranchStatus("Jet_PT",1);
  outT->SetBranchStatus("Jet_PT_up",1);
  outT->SetBranchStatus("Jet_PT_down",1);
      
  outT->SetBranchStatus("Jet_Phi",1);
  outT->SetBranchStatus("Jet_Phi_up",1);
  outT->SetBranchStatus("Jet_Phi_down",1);
      
  outT->SetBranchStatus("Jet_Eta",1);
  outT->SetBranchStatus("Jet_Eta_up",1);
  outT->SetBranchStatus("Jet_Eta_down",1);
      
  outT->SetBranchStatus("pTHem1",1);
  outT->SetBranchStatus("pTHem1_up",1);
  outT->SetBranchStatus("pTHem1_down",1);
      
  outT->SetBranchStatus("pTHem2",1);
  outT->SetBranchStatus("pTHem2_up",1);
  outT->SetBranchStatus("pTHem2_down",1);
  
  outT->SetBranchStatus("etaHem1",1);
  outT->SetBranchStatus("etaHem1_up",1);
  outT->SetBranchStatus("etaHem1_down",1);
  
  outT->SetBranchStatus("etaHem2",1);
  outT->SetBranchStatus("etaHem2_up",1);
  outT->SetBranchStatus("etaHem2_down",1);
  
  outT->SetBranchStatus("phiHem1",1);
  outT->SetBranchStatus("phiHem1_up",1);
  outT->SetBranchStatus("phiHem1_down",1);
  
  outT->SetBranchStatus("phiHem2",1);
  outT->SetBranchStatus("phiHem2_up",1);
  outT->SetBranchStatus("phiHem2_down",1);
      
  outT->SetBranchStatus("metCorrX",1);
  outT->SetBranchStatus("metX",1);
  outT->SetBranchStatus("metCorrX_up",1);
  outT->SetBranchStatus("metCorrX_down",1);
  
  outT->SetBranchStatus("metCorrY",1);
  outT->SetBranchStatus("metY",1);
  outT->SetBranchStatus("metCorrY_up",1);
  outT->SetBranchStatus("metCorrY_down",1);
  //SMS INFO
  outT->SetBranchStatus("mssm", 1);
  
  ///////////////////////////////
  ///////////Addresses///////////
  ///////////////////////////////
      
  outT->SetBranchAddress("pu_w", &pu_w);
  outT->SetBranchAddress("ISR", &ISR);
  outT->SetBranchAddress("ISR_up", &ISR_up);
  outT->SetBranchAddress("ISR_down", &ISR_down);
  
  outT->SetBranchAddress("MR", mr);
  outT->SetBranchAddress("MR_up", mr_up);
  outT->SetBranchAddress("MR_down", mr_down);
  
  outT->SetBranchAddress("RSQ", rsq);
  outT->SetBranchAddress("RSQ_up", rsq_up);
  outT->SetBranchAddress("RSQ_down", rsq_down);
  
  outT->SetBranchAddress("nBtag", &btag);
  outT->SetBranchAddress("BOX_NUM", &box);
  outT->SetBranchAddress("N_Jets", &N_Jets);
      
  outT->SetBranchAddress("Jet_PT", Jet_PT);
  outT->SetBranchAddress("Jet_PT_up", Jet_PT_up);
  outT->SetBranchAddress("Jet_PT_down", Jet_PT_down);
      
  outT->SetBranchAddress("Jet_Phi", Jet_Phi);
  outT->SetBranchAddress("Jet_Phi_up", Jet_Phi_up);
  outT->SetBranchAddress("Jet_Phi_down", Jet_Phi_down);
      
  outT->SetBranchAddress("Jet_Eta", Jet_Eta);
  outT->SetBranchAddress("Jet_Eta_up", Jet_Eta_up);
  outT->SetBranchAddress("Jet_Eta_down", Jet_Eta_down);
      
  outT->SetBranchAddress("pTHem1", &pTHem1);
  outT->SetBranchAddress("pTHem1_up", &pTHem1_up);
  outT->SetBranchAddress("pTHem1_down", &pTHem1_down);
      
  outT->SetBranchAddress("pTHem2", &pTHem2);
  outT->SetBranchAddress("pTHem2_up", &pTHem2_up);
  outT->SetBranchAddress("pTHem2_down", &pTHem2_down);
      
  outT->SetBranchAddress("etaHem1", &etaHem1);
  outT->SetBranchAddress("etaHem1_up", &etaHem1_up);
  outT->SetBranchAddress("etaHem1_down", &etaHem1_down);
      
  outT->SetBranchAddress("etaHem2", &etaHem2);
  outT->SetBranchAddress("etaHem2_up", &etaHem2_up);
  outT->SetBranchAddress("etaHem2_down", &etaHem2_down);
      
  outT->SetBranchAddress("phiHem1", &phiHem1);
  outT->SetBranchAddress("phiHem1_up", &phiHem1_up);
  outT->SetBranchAddress("phiHem1_down", &phiHem1_down);
      
  outT->SetBranchAddress("phiHem2", &phiHem2);
  outT->SetBranchAddress("phiHem2_up", &phiHem2_up);
  outT->SetBranchAddress("phiHem2_down", &phiHem2_down);
      
  outT->SetBranchAddress("metCorrX", metCorrX);
  outT->SetBranchAddress("metX", metX);
  outT->SetBranchAddress("metCorrX_up", metCorrX_up);
  outT->SetBranchAddress("metCorrX_down", metCorrX_down);
  
  outT->SetBranchAddress("metCorrY", metCorrY);
  outT->SetBranchAddress("metY", metY);
  outT->SetBranchAddress("metCorrY_up", metCorrY_up);
  outT->SetBranchAddress("metCorrY_down", metCorrY_down);
  
  //SMS info
  outT->SetBranchAddress("mssm", mssm);
  int N_out = outT->GetEntries();
  
  for(int j = 0; j < N_out; j++){
    outT->GetEntry(j);
    //GETTING SMS INFO
    std::stringstream ss_aux;
    ss_aux << mssm[0] << " " << mssm[1];
    std::string s_aux = ss_aux.str();
    double hlt_w = HLTscale(mr[2], rsq[2], hlt);
    //Nominal
    TLorentzVector j1;
    TLorentzVector j2;
    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
    double Dphi = j1.DeltaPhi(j2);
  
    double MET = sqrt(metX[2]*metX[2] + metY[2]*metY[2]);
    
    if(mr[2] >= 200.0 && rsq[2] >= 0.5 && btag == 0 && box == 0 && fabs(Dphi) < 2.5){
      if(mr[2] > 200.0 && mr[2] <= 300.0 ){
	if(rsq[2] > 1.2)rsq[2] = 1.1;
	sms_mass_pair[s_aux].h_rsq[0]->Fill(rsq[2], hlt_w*sms_mass_pair[s_aux].scaleF*pu_w*ISR);
	sms_mass_pair[s_aux].h_met[0]->Fill(MET, hlt_w*sms_mass_pair[s_aux].scaleF*pu_w*ISR);
	sms_mass_pair[s_aux].h_njets[0]->Fill(N_Jets, hlt_w*sms_mass_pair[s_aux].scaleF*pu_w*ISR);
	/*
	N_passed_ISR[0] += hlt_w*pu_w*ISR;
	if(Jet_PT[0] > 110.0 && N_Jets == 2 && MET > 500.0){
	  N_passed_ISR_MonoJet[0] += hlt_w*pu_w*ISR;
	}
	*/
      }else if(mr[2] > 300.0 && mr[2] <= 400.0){
	if(rsq[2] > 1.2)rsq[2] = 1.1;
	sms_mass_pair[s_aux].h_rsq[1]->Fill(rsq[2], hlt_w*sms_mass_pair[s_aux].scaleF*pu_w*ISR);
	sms_mass_pair[s_aux].h_met[1]->Fill(MET, hlt_w*sms_mass_pair[s_aux].scaleF*pu_w*ISR);
	sms_mass_pair[s_aux].h_njets[1]->Fill(N_Jets, hlt_w*sms_mass_pair[s_aux].scaleF*pu_w*ISR);
	/*
	N_passed_ISR[1] += hlt_w*pu_w*ISR;
	if(Jet_PT[0] > 110.0 && N_Jets == 2 && MET > 500.0){
	  N_passed_ISR_MonoJet[1] += hlt_w*pu_w*ISR;
	}
	*/
      }else if(mr[2] > 400.0 && mr[2] <= 600.0){
	if(rsq[2] > 1.2)rsq[2] = 1.1;
	sms_mass_pair[s_aux].h_rsq[2]->Fill(rsq[2], hlt_w*sms_mass_pair[s_aux].scaleF*pu_w*ISR);
	sms_mass_pair[s_aux].h_met[2]->Fill(MET, hlt_w*sms_mass_pair[s_aux].scaleF*pu_w*ISR);
	sms_mass_pair[s_aux].h_njets[2]->Fill(N_Jets, hlt_w*sms_mass_pair[s_aux].scaleF*pu_w*ISR);
	/*
	N_passed_ISR[2] += hlt_w*pu_w*ISR;
	if(Jet_PT[0] > 110.0 && N_Jets == 2 && MET > 500.0){
	  N_passed_ISR_MonoJet[2] += hlt_w*pu_w*ISR;
	  }
	*/
      }else if(mr[2] > 600.0 && mr[2] <= 3500.0){
	if(rsq[2] > 1.2)rsq[2] = 1.1;
	sms_mass_pair[s_aux].h_rsq[3]->Fill(rsq[2], hlt_w*sms_mass_pair[s_aux].scaleF*pu_w*ISR);
	sms_mass_pair[s_aux].h_met[3]->Fill(MET, hlt_w*sms_mass_pair[s_aux].scaleF*pu_w*ISR);
	sms_mass_pair[s_aux].h_njets[3]->Fill(N_Jets, hlt_w*sms_mass_pair[s_aux].scaleF*pu_w*ISR);
	/*
	N_passed_ISR[3] += hlt_w*pu_w*ISR;
	if(Jet_PT[0] > 110.0 && N_Jets == 2 && MET > 500.0){
	  N_passed_ISR_MonoJet[3] += hlt_w*pu_w*ISR;
	}
	*/
      }
      
      /*
      N_passed += hlt_w*pu_w;
      
      N_passed_ISR[4] += hlt_w*pu_w*ISR;
      if(Jet_PT[0] > 110.0 && N_Jets == 2 && MET > 500.0){
	N_passed_ISR_MonoJet[4] += hlt_w*pu_w*ISR;
      }
      */
    }
      
    //ISR UP
    if(mr[2] >= 200.0 && rsq[2] >= 0.5 && btag == 0 && box == 0 && fabs(Dphi) < 2.5){
      if(mr[2] > 200.0 && mr[2] <= 300.0 ){
	if(rsq[2] > 1.2)rsq[2] = 1.1;
	sms_mass_pair[s_aux].h_rsq_ISR_up[0]->Fill(rsq[2], hlt_w*sms_mass_pair[s_aux].scaleF_up*pu_w*ISR_up);
      }else if(mr[2] > 300.0 && mr[2] <= 400.0){
	if(rsq[2] > 1.2)rsq[2] = 1.1;
	sms_mass_pair[s_aux].h_rsq_ISR_up[1]->Fill(rsq[2], hlt_w*sms_mass_pair[s_aux].scaleF_up*pu_w*ISR_up);
      }else if(mr[2] > 400.0 && mr[2] <= 600.0){
	if(rsq[2] > 1.2)rsq[2] = 1.1;
	sms_mass_pair[s_aux].h_rsq_ISR_up[2]->Fill(rsq[2], hlt_w*sms_mass_pair[s_aux].scaleF_up*pu_w*ISR_up);
      }else if(mr[2] > 600.0 && mr[2] <= 3500.0){
	if(rsq[2] > 1.2)rsq[2] = 1.1;
	sms_mass_pair[s_aux].h_rsq_ISR_up[3]->Fill(rsq[2], hlt_w*sms_mass_pair[s_aux].scaleF_up*pu_w*ISR_up);
      }
    }
    
    //ISR DOWN
    if(mr[2] >= 200.0 && rsq[2] >= 0.5 && btag == 0 && box == 0 && fabs(Dphi) < 2.5){
      if(mr[2] > 200.0 && mr[2] <= 300.0 ){
	if(rsq[2] > 1.2)rsq[2] = 1.1;
	sms_mass_pair[s_aux].h_rsq_ISR_down[0]->Fill(rsq[2], hlt_w*sms_mass_pair[s_aux].scaleF_down*pu_w*ISR_down);
      }else if(mr[2] > 300.0 && mr[2] <= 400.0){
	if(rsq[2] > 1.2)rsq[2] = 1.1;
	sms_mass_pair[s_aux].h_rsq_ISR_down[1]->Fill(rsq[2], hlt_w*sms_mass_pair[s_aux].scaleF_down*pu_w*ISR_down);
      }else if(mr[2] > 400.0 && mr[2] <= 600.0){
	if(rsq[2] > 1.2)rsq[2] = 1.1;
	sms_mass_pair[s_aux].h_rsq_ISR_down[2]->Fill(rsq[2], hlt_w*sms_mass_pair[s_aux].scaleF_down*pu_w*ISR_down);
      }else if(mr[2] > 600.0 && mr[2] <= 3500.0){
	if(rsq[2] > 1.2)rsq[2] = 1.1;
	sms_mass_pair[s_aux].h_rsq_ISR_down[3]->Fill(rsq[2], hlt_w*sms_mass_pair[s_aux].scaleF_down*pu_w*ISR_down);
      }
    }
    
    //JES UP
    j1.SetPtEtaPhiE(pTHem1_up, etaHem1_up, phiHem1_up, pTHem1_up*cosh(etaHem1_up));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2_up, etaHem2_up, phiHem2_up, pTHem2_up*cosh(etaHem2_up));//Hemisphere
    Dphi = j1.DeltaPhi(j2);
    if(mr_up[2] >= 200.0 && rsq_up[2] >= 0.5 && btag == 0 && box == 0 && fabs(Dphi) < 2.5){
      if(mr_up[2] > 200.0 && mr_up[2] <= 300.0 ){
	if(rsq_up[2] > 1.2)rsq_up[2] = 1.1;
	sms_mass_pair[s_aux].h_rsq_JES_up[0]->Fill(rsq_up[2], hlt_w*sms_mass_pair[s_aux].scaleF*pu_w*ISR);
      }else if(mr_up[2] > 300.0 && mr_up[2] <= 400.0){
	if(rsq_up[2] > 1.2)rsq_up[2] = 1.1;
	sms_mass_pair[s_aux].h_rsq_JES_up[1]->Fill(rsq_up[2], hlt_w*sms_mass_pair[s_aux].scaleF*pu_w*ISR);
      }else if(mr_up[2] > 400.0 && mr_up[2] <= 600.0){
	if(rsq_up[2] > 1.2)rsq_up[2] = 1.1;
	sms_mass_pair[s_aux].h_rsq_JES_up[2]->Fill(rsq_up[2], hlt_w*sms_mass_pair[s_aux].scaleF*pu_w*ISR);
      }else if(mr_up[2] > 600.0 && mr_up[2] <= 3500.0){
	if(rsq_up[2] > 1.2)rsq_up[2] = 1.1;
	sms_mass_pair[s_aux].h_rsq_JES_up[3]->Fill(rsq_up[2], hlt_w*sms_mass_pair[s_aux].scaleF*pu_w*ISR);
      }
    }
  
    //JES DOWN
    j1.SetPtEtaPhiE(pTHem1_down, etaHem1_down, phiHem1_down, pTHem1_down*cosh(etaHem1_down));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2_down, etaHem2_down, phiHem2_down, pTHem2_down*cosh(etaHem2_down));//Hemisphere
    Dphi = j1.DeltaPhi(j2);
    if(mr_down[2] >= 200.0 && rsq_down[2] >= 0.5 && btag == 0 && box == 0 && fabs(Dphi) < 2.5){
      if(mr_down[2] > 200.0 && mr_down[2] <= 300.0 ){
	if(rsq_down[2] > 1.2)rsq_down[2] = 1.1;
	sms_mass_pair[s_aux].h_rsq_JES_down[0]->Fill(rsq_down[2], hlt_w*sms_mass_pair[s_aux].scaleF*pu_w*ISR);
      }else if(mr_down[2] > 300.0 && mr_down[2] <= 400.0){
	if(rsq_down[2] > 1.2)rsq_down[2] = 1.1;
	sms_mass_pair[s_aux].h_rsq_JES_down[1]->Fill(rsq_down[2], hlt_w*sms_mass_pair[s_aux].scaleF*pu_w*ISR);
      }else if(mr_down[2] > 400.0 && mr_down[2] <= 600.0){
	if(rsq_down[2] > 1.2)rsq_down[2] = 1.1;
	sms_mass_pair[s_aux].h_rsq_JES_down[2]->Fill(rsq_down[2], hlt_w*sms_mass_pair[s_aux].scaleF*pu_w*ISR);
      }else if(mr_down[2] > 600.0 && mr_down[2] <= 3500.0){
	if(rsq_down[2] > 1.2)rsq_down[2] = 1.1;
	sms_mass_pair[s_aux].h_rsq_JES_down[3]->Fill(rsq_down[2], hlt_w*sms_mass_pair[s_aux].scaleF*pu_w*ISR);
      }
    }
    
  }
      
  /*
	
      //std::cout << "Razor/MONOJET: " << (N_passed_ISR_MonoJet/N_passed_ISR)*100 << "%" << std::endl;
      //std::cout << "MONOJET ONLY: " <<  (NP_ISR_MonoJetONLY/Gen_Evts_isr)*100 << "%" << std::endl;
      //outfile << dm_sample << " & " << sample_eff*100 << "\\%" << "\\" << "\\" << "\n";
      double sample_eff[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
      for(int l = 0; l < 5; l++)sample_eff[l] = (N_passed_ISR_MonoJet[l]/N_passed_ISR[l])*100.0;
      outfile << dm_sample << " & " << sample_eff[4] << "\\%" << " & " << sample_eff[0] << "\\%" 
	      << " & " << sample_eff[1] << "\\%"
	      << " & " << sample_eff[2] << "\\%"
	      << " & " << sample_eff[3] << "\\%" << "\\\\"
	      << "\n";
      outfile << "\\hline" << std::endl;
      
      double sample_Eff_sample[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
      for(int l = 0; l < 5; l++)sample_Eff_sample[l] = (N_passed_ISR[l]/Gen_Evts_isr)*100.0;
      outfile2 << dm_sample << " & " << sample_Eff_sample[4] << "\\%" << " & " << sample_Eff_sample[0] << "\\%"
              << " & " << sample_Eff_sample[1] << "\\%"
              << " & " << sample_Eff_sample[2] << "\\%"
              << " & " << sample_Eff_sample[3] << "\\%" << "\\\\"
              << "\n";
      outfile2 << "\\hline" << std::endl;
      
      outfile1 << dm_sample << " & " << N_passed_ISR[4] << " & " << N_passed_ISR[0]
	       << " & " << N_passed_ISR[1]
	       << " & " << N_passed_ISR[2]
	       << " & " << N_passed_ISR[3]  << "\\\\"
	       << "\n";
      outfile1 << "\\hline" << std::endl;
      
      for(int i = 0; i < 4; i++){
	for(int j = 1; j <= s_up[xs_counter][i]->GetNbinsX(); j++){
	  s_up[xs_counter][i]->SetBinContent(j, h_rsq[xs_counter][i]->GetBinContent(j)+h_rsq[xs_counter][i]->GetBinError(j));
	  s_down[xs_counter][i]->SetBinContent(j, h_rsq[xs_counter][i]->GetBinContent(j)-h_rsq[xs_counter][i]->GetBinError(j));
	}
      }
  */
  
  
  for( auto tmp : sms_mass_pair){
    std::stringstream ss_aux;
    std::stringstream ss_aux2;
    ss_aux.str(tmp.first);
    int mass1, mass2;
    ss_aux >> mass1 >> mass2;
    ss_aux.str("");
    ss_aux2 << mass1 << "_" << mass2;
    TString sms_s = ss_aux2.str().c_str();
    /////1D files
    TFile* fo;
    for(int i = 0; i < 4; i++){
      TString s1 = sms_s;
      s1 = s1 + Form("_combine_rsq_cat%d.root",i+1);
      TString data_card_name1;
      data_card_name1 = "CombineFilesSMS/" + sms_s + Form("_rsq_cat%d.txt",i+1);
      std::ofstream data_card_f1(data_card_name1);
      data_card_f1 << "imax 1\njmax 4\nkmax 6\n";
      data_card_f1 << "------------------------------------------------------------------------------------------\n";
      data_card_f1 << "shapes * *\t" << s1 << "\t\t$PROCESS\t$PROCESS_$SYSTEMATIC\n";
      data_card_f1 << "------------------------------------------------------------------------------------------\n";
      data_card_f1 << "Observation\t" << data_N[i] << "\n";
      data_card_f1 << "------------------------------------------------------------------------------------------\n";
      data_card_f1 << "bin\t\tb1\t\tb1\t\tb1\t\tb1\t\tb1\n";
      data_card_f1 << "process\t\tsignal_rsq\ttt_rsq\t\tdy_rsq\t\tz_rsq\t\tw_rsq\n";
      data_card_f1 << "process\t\t0\t\t1\t\t2\t\t3\t\t4\n";
      data_card_f1 << "rate\t\t"<< tmp.second.h_rsq[i]->Integral() <<"\t\t"<< tt_N[i] <<"\t\t" << dy_N[i] << "\t" << z_N[i] << "\t\t" << w_N[i] << "\n";
      data_card_f1 << "------------------------------------------------------------------------------------------\n";
      data_card_f1 << "lumi\tlnN\t1.026\t\t1.0\t\t1.0\t\t1.0\t\t1.0\n";
      //data_card_f1 << "alpha\tshape\t1\t\t-\t\t-\t\t-\t\t-\n";
      data_card_f1 << "Isr\tshape\t1\t\t-\t\t-\t\t-\t\t-\n";
      data_card_f1 << "Jes\tshape\t1\t\t-\t\t-\t\t-\t\t-\n";
      data_card_f1 << "Pdf\tshape\t1\t\t-\t\t-\t\t-\t\t-\n";
      data_card_f1 << "beta"<<i+1<<"\tshape\t-\t\t1\t\t-\t\t-\t\t-\n";
      data_card_f1 << "gamma"<<i+1<<"\tshape\t-\t\t-\t\t1\t\t1\t\t1\n";
    
      //data_card_f1 << "gamma"<<i+1<<"\tshape\t-\t\t1\t\t1\t\t1\t\t1\n";
      
      //data_card_f1 << "gamma"<<i+1<<"\tshape\t-\t\t-\t\t1\t\t-\t\t-\n";
      //data_card_f1 << "epsilon"<<i+1<<"\tshape\t-\t\t-\t\t-\t\t1\t\t-\n";
      //data_card_f1 << "delta"<<i+1<<"\tshape\t-\t\t-\t\t-\t\t-\t\t1\n";
      data_card_f1.close();

      fo = new TFile("CombineFilesSMS/"+s1, "RECREATE");
      
      data[i]->Write("data_obs");
	
      tmp.second.h_rsq_PDF_up[i] = new TH1F( *tmp.second.h_rsq[i] );
      tmp.second.h_rsq_PDF_down[i] = new TH1F( *tmp.second.h_rsq[i] );
      
      for(int bin = 1; bin <= tmp.second.h_rsq_PDF_up[i]->GetNbinsX(); bin++){
	double pdf_up = (tmp.second.h_rsq_PDF_up[i]->GetBinContent(bin))*(1.0 + tmp.second.pdf_acc[i]->GetBinContent(bin));
	double pdf_down = (tmp.second.h_rsq_PDF_up[i]->GetBinContent(bin))*(1.0 - tmp.second.pdf_acc[i]->GetBinContent(bin));
	tmp.second.h_rsq_PDF_up[i]->SetBinContent(bin, pdf_up);
	tmp.second.h_rsq_PDF_down[i]->SetBinContent(bin, pdf_down);
      }
      
      tmp.second.h_rsq[i]->Write("signal_rsq");
      //s_up[i]->Write("signal_rsq_alphaUp");
      //s_down[i]->Write("signal_rsq_alphaDown");
      
      tmp.second.h_rsq_ISR_up[i]->Write("signal_rsq_IsrUp");
      tmp.second.h_rsq_ISR_down[i]->Write("signal_rsq_IsrDown");
      
      tmp.second.h_rsq_JES_up[i]->Write("signal_rsq_JesUp");
      tmp.second.h_rsq_JES_down[i]->Write("signal_rsq_JesDown");
      
      
      tmp.second.h_rsq_PDF_up[i]->Write("signal_rsq_PdfUp");
      tmp.second.h_rsq_PDF_down[i]->Write("signal_rsq_PdfDown");
      
      TString SysN = TString(Form("dy_rsq_gamma%dUp",i+1));
      dy[i]->Write("dy_rsq");
      dy_up[i]->Write(SysN);
      SysN = TString(Form("dy_rsq_gamma%dDown",i+1));
      dy_down[i]->Write(SysN);
      
      z[i]->Write("z_rsq");
      SysN = TString(Form("z_rsq_gamma%dUp",i+1));
      z_up[i]->Write(SysN);
      SysN = TString(Form("z_rsq_gamma%dDown",i+1));
      z_down[i]->Write(SysN);
      
      
      w[i]->Write("w_rsq");
      SysN = TString(Form("w_rsq_gamma%dUp",i+1));
      w_up[i]->Write(SysN);
      SysN = TString(Form("w_rsq_gamma%dDown",i+1));
      w_down[i]->Write(SysN);
      
      tt[i]->Write("tt_rsq");
      SysN = TString(Form("tt_rsq_beta%dUp",i+1));
      //SysN = TString(Form("tt_rsq_gamma%dUp",i+1));
      tt_up[i]->Write(SysN);
      SysN = TString(Form("tt_rsq_beta%dDown",i+1));
      //SysN = TString(Form("tt_rsq_gamma%dDown",i+1));
      tt_down[i]->Write(SysN);
      
      bkg[i]->Write("Pred");
      fo->Close();
      
    } 
  }
  /*
  mfile0.close();
  outfile << "\\end{tabular}\n\\end{center}\n\\label{default}\n\\end{table}\n";
  outfile.close();
  outfile1.close();
  outfile2.close();

  TFile* f4 = new TFile("SIGNAL_MET_NJETS.root", "RECREATE");
  for(int j = 0; j < 24; j++){
    for(int i = 0; i < 4; i++){
      tmp.second.h_met[j][i]->Write();
      tmp.second.h_njets[j][i]->Write();
    }
  }
  f4->Close();
  */
  return 0;
}
