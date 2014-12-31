/*
Generates the input for combine-tool limit setting code
including data card and the root file with the histograms
for the systematics, signal is scale to the correct LUMI
and includes the ISR, PDF, JES systematic, as well as the
PU re-weighting, and the HLT weight for the turn-on curve
*/
#include <iostream>
#include "math.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "hlt.hh"
#include <fstream>
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include <vector>
#include "DM_1DRatio.hh"

const int r2Bins = 5;
const float BaseDM::RSQ_BinArr[] = {0.5, 0.6, 0.725, 0.85,  2.50};
const float BaseDM::MR_BinArr[] = {200., 300., 400., 600.,  3500.};

int main(){

  gROOT->Reset();


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
  TFile* f2 = new TFile("trigger/hlt_eff_HTMHT_0mu_BOX_0bTag_Final.root");
  TEfficiency* hlt = (TEfficiency*)f2->Get("Eff2d");

  //Creating Histos for Shape Analysis
  
  TLegend* leg;

  TH1F* h_rsq[24][4];
  TH1F* N_j[24];
  TH1F* j1_Pt[24];
  TH1F* j2_Pt[24];
  TH1F* mR[24];//wacky name to avoid overlap
  TH1F* rSQ[24];//wacky name to avoid overlap
  TH1F* dPHI[24];
  TH2F* mR_rSQ[24];
  TH1F* mET[24];
  
  TString sn;
  for(int j = 0; j < 24; j++){
    sn = TString(Form("signal%d_NJ",j));
    N_j[j] = new TH1F(sn, sn, 6, 2, 8);
    sn = TString(Form("signal%d_j1PT",j));
    j1_Pt[j] = new TH1F(sn, sn, 50, 50, 1500);
    sn = TString(Form("signal%d_j2PT",j));
    j2_Pt[j] = new TH1F(sn, sn, 50, 50, 1500);
    
    sn = TString(Form("signal%d_mR",j));
    mR[j] = new TH1F(sn, sn, 50, 50, 2500);
    sn = TString(Form("signal%d_rSQ",j));
    rSQ[j] = new TH1F(sn, sn, 50, 0.0, 1.5);

    sn = TString(Form("signal%d_dPHI",j));
    dPHI[j] = new TH1F(sn, sn, 50, -3.3, 3.3);
    
    sn = TString(Form("signal%d_mR_rSQ",j));
    mR_rSQ[j] = new TH2F(sn, sn, 50, 50, 2500, 50, 0.0, 1.5);

    sn = TString(Form("signal%d_MET",j));
    mET[j] = new TH1F(sn, sn, 100, 0, 1500);

    for(int i = 0; i < 4; i++){
      sn = TString(Form("signal%d_cat%d",j,i));
      h_rsq[j][i] = new TH1F(sn, sn, r2B[i], v.at(i));
    }
  }

  //Here the program starts
  
  //std::ifstream mfile0("list_DM_BIS.list");
  std::ifstream mfile0("list_DM_BugFixed.list");
  std::ofstream outfile("eff_table_normal_R2_0p5_MR_200_Dphi_B_2p5_New.tex");
  
  outfile << "\\begin{table}[htdp]\n\\caption{default}\n\\begin{center}\n\\begin{tabular}{|c|c|}\n\\hline\n";
  
  std::string fname0;
  std::cout.precision(16);
  int xs_counter = 0;
  
  TFile* f_acc;
  
  if (mfile0.is_open()){
    while ( mfile0.good() ){
      mfile0 >> fname0;
      if(mfile0.eof())break;
      std::cout << fname0 << std::endl;
      int low_ = fname0.find("DMm");
      int high_ = fname0.find("_testMC_0.root") - low_;
      
      std::string dm_sample = fname0.substr(low_,high_);
      std::cout << "============ " << dm_sample << " ==================" << std::endl;
      
      TFile* f = new TFile(fname0.c_str());
      TTree* eff = (TTree*)f->Get("effTree");
      TTree* out = (TTree*)f->Get("outTree");
      
      double mr[4], rsq[4], Jet_PT[20], Jet_Eta[20], Jet_Phi[20], metCorrX[4], metCorrY[4];
      double mr_up[4], rsq_up[4], Jet_PT_up[20], Jet_Eta_up[20], Jet_Phi_up[20], metCorrX_up[4], metCorrY_up[4];
      double mr_down[4], rsq_down[4], Jet_PT_down[20], Jet_Eta_down[20], Jet_Phi_down[20], metCorrX_down[4], metCorrY_down[4];
      double metX[4], metY[4];
      
      double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2;
      double pTHem1_up, pTHem2_up, etaHem1_up, etaHem2_up, phiHem1_up, phiHem2_up;
      double pTHem1_down, pTHem2_down, etaHem1_down, etaHem2_down, phiHem1_down, phiHem2_down;
      
      int btag, box, N_Jets;
      double pu_w, ISR, ISR_up, ISR_down;
      
      double Npassed_In, Npassed_ISR, Npassed_ISR_up, Npassed_ISR_down;
      
      eff->SetBranchStatus("*", 0);
      eff->SetBranchStatus("Npassed_In", 1);
      eff->SetBranchStatus("Npassed_ISR", 1);
      eff->SetBranchStatus("Npassed_ISR_up", 1);
      eff->SetBranchStatus("Npassed_ISR_down", 1);
      
      eff->SetBranchAddress("Npassed_In", &Npassed_In);
      eff->SetBranchAddress("Npassed_ISR", &Npassed_ISR);
      eff->SetBranchAddress("Npassed_ISR_up", &Npassed_ISR_up);
      eff->SetBranchAddress("Npassed_ISR_down", &Npassed_ISR_down);
      
      
      int N_eff = eff->GetEntries();
      int Gen_Evts = 0;
      int Gen_Evts_isr = 0;
      int Gen_Evts_isrUp = 0;
      int Gen_Evts_isrDown = 0;
      
      for(int i = 0; i < N_eff; i++){
	eff->GetEntry(i);
	Gen_Evts += Npassed_In;
	Gen_Evts_isr += Npassed_ISR;
	Gen_Evts_isrUp += Npassed_ISR_up;
	Gen_Evts_isrDown += Npassed_ISR_down;
      }
      
      std::cout << "Gen_Events: " << Gen_Evts << std::endl;
      std::cout << "Gen_ISR: " << Gen_Evts_isr << std::endl;
      std::cout << "Gen_ISR_Up: " << Gen_Evts_isrUp << std::endl;
      std::cout << "Gen_ISR_Down: " << Gen_Evts_isrDown << std::endl;
      
      out->SetBranchStatus("*", 0);
      
      out->SetBranchStatus("pu_w", 1);
      out->SetBranchStatus("ISR", 1);
      out->SetBranchStatus("ISR_up", 1);
      out->SetBranchStatus("ISR_down", 1);
      
      out->SetBranchStatus("MR", 1);
      out->SetBranchStatus("MR_up", 1);
      out->SetBranchStatus("MR_down", 1);

      out->SetBranchStatus("RSQ",1);
      out->SetBranchStatus("RSQ_up",1);
      out->SetBranchStatus("RSQ_down",1);
      
      out->SetBranchStatus("nBtag", 1);
      out->SetBranchStatus("BOX_NUM",1);
      out->SetBranchStatus("N_Jets",1);
      
      out->SetBranchStatus("Jet_PT",1);
      out->SetBranchStatus("Jet_PT_up",1);
      out->SetBranchStatus("Jet_PT_down",1);
      
      out->SetBranchStatus("Jet_Phi",1);
      out->SetBranchStatus("Jet_Phi_up",1);
      out->SetBranchStatus("Jet_Phi_down",1);
      
      out->SetBranchStatus("Jet_Eta",1);
      out->SetBranchStatus("Jet_Eta_up",1);
      out->SetBranchStatus("Jet_Eta_down",1);
      
      out->SetBranchStatus("pTHem1",1);
      out->SetBranchStatus("pTHem1_up",1);
      out->SetBranchStatus("pTHem1_down",1);
      
      out->SetBranchStatus("pTHem2",1);
      out->SetBranchStatus("pTHem2_up",1);
      out->SetBranchStatus("pTHem2_down",1);
       
      out->SetBranchStatus("etaHem1",1);
      out->SetBranchStatus("etaHem1_up",1);
      out->SetBranchStatus("etaHem1_down",1);
      
      out->SetBranchStatus("etaHem2",1);
      out->SetBranchStatus("etaHem2_up",1);
      out->SetBranchStatus("etaHem2_down",1);
      
      out->SetBranchStatus("phiHem1",1);
      out->SetBranchStatus("phiHem1_up",1);
      out->SetBranchStatus("phiHem1_down",1);
      
      out->SetBranchStatus("phiHem2",1);
      out->SetBranchStatus("phiHem2_up",1);
      out->SetBranchStatus("phiHem2_down",1);
      
      out->SetBranchStatus("metCorrX",1);
      out->SetBranchStatus("metX",1);
      out->SetBranchStatus("metCorrX_up",1);
      out->SetBranchStatus("metCorrX_down",1);

      out->SetBranchStatus("metCorrY",1);
      out->SetBranchStatus("metY",1);
      out->SetBranchStatus("metCorrY_up",1);
      out->SetBranchStatus("metCorrY_down",1);

      ///////////////////////////////
      ///////////Addresses///////////
      ///////////////////////////////
      
      out->SetBranchAddress("pu_w", &pu_w);
      out->SetBranchAddress("ISR", &ISR);
      out->SetBranchAddress("ISR_up", &ISR_up);
      out->SetBranchAddress("ISR_down", &ISR_down);

      out->SetBranchAddress("MR", mr);
      out->SetBranchAddress("MR_up", mr_up);
      out->SetBranchAddress("MR_down", mr_down);
      
      out->SetBranchAddress("RSQ", rsq);
      out->SetBranchAddress("RSQ_up", rsq_up);
      out->SetBranchAddress("RSQ_down", rsq_down);
      
      out->SetBranchAddress("nBtag", &btag);
      out->SetBranchAddress("BOX_NUM", &box);
      out->SetBranchAddress("N_Jets", &N_Jets);
      
      out->SetBranchAddress("Jet_PT", Jet_PT);
      out->SetBranchAddress("Jet_PT_up", Jet_PT_up);
      out->SetBranchAddress("Jet_PT_down", Jet_PT_down);
      
      out->SetBranchAddress("Jet_Phi", Jet_Phi);
      out->SetBranchAddress("Jet_Phi_up", Jet_Phi_up);
      out->SetBranchAddress("Jet_Phi_down", Jet_Phi_down);
      
      out->SetBranchAddress("Jet_Eta", Jet_Eta);
      out->SetBranchAddress("Jet_Eta_up", Jet_Eta_up);
      out->SetBranchAddress("Jet_Eta_down", Jet_Eta_down);
      
      out->SetBranchAddress("pTHem1", &pTHem1);
      out->SetBranchAddress("pTHem1_up", &pTHem1_up);
      out->SetBranchAddress("pTHem1_down", &pTHem1_down);
      
      out->SetBranchAddress("pTHem2", &pTHem2);
      out->SetBranchAddress("pTHem2_up", &pTHem2_up);
      out->SetBranchAddress("pTHem2_down", &pTHem2_down);
      
      out->SetBranchAddress("etaHem1", &etaHem1);
      out->SetBranchAddress("etaHem1_up", &etaHem1_up);
      out->SetBranchAddress("etaHem1_down", &etaHem1_down);
      
      out->SetBranchAddress("etaHem2", &etaHem2);
      out->SetBranchAddress("etaHem2_up", &etaHem2_up);
      out->SetBranchAddress("etaHem2_down", &etaHem2_down);
      
      out->SetBranchAddress("phiHem1", &phiHem1);
      out->SetBranchAddress("phiHem1_up", &phiHem1_up);
      out->SetBranchAddress("phiHem1_down", &phiHem1_down);
      
      out->SetBranchAddress("phiHem2", &phiHem2);
      out->SetBranchAddress("phiHem2_up", &phiHem2_up);
      out->SetBranchAddress("phiHem2_down", &phiHem2_down);
      
      out->SetBranchAddress("metCorrX", metCorrX);
      out->SetBranchAddress("metX", metX);
      out->SetBranchAddress("metCorrX_up", metCorrX_up);
      out->SetBranchAddress("metCorrX_down", metCorrX_down);
      
      out->SetBranchAddress("metCorrY", metCorrY);
      out->SetBranchAddress("metY", metY);
      out->SetBranchAddress("metCorrY_up", metCorrY_up);
      out->SetBranchAddress("metCorrY_down", metCorrY_down);
      
      int N_out = out->GetEntries();
      double Lumi = 18.836;//Jan22Rereco
      double scaleF = Lumi*1000./Gen_Evts_isr;
      double scaleF_up = Lumi*1000./Gen_Evts_isrUp;
      double scaleF_down = Lumi*1000./Gen_Evts_isrDown;
      

      double N_passed = 0.0;
      double N_passed_ISR = 0.0;
      double N_passed_ISR_up = 0.0;
      double N_passed_ISR_down = 0.0;
      double N_passed_JES_up = 0.0;
      double N_passed_JES_down = 0.0;

      double N_passed_ISR_MonoJet = 0.0;
      double NP_ISR_MonoJetONLY = 0.0;
      
      for(int j = 0; j < N_out; j++){
	out->GetEntry(j);
	double hlt_w = HLTscale(mr[2], rsq[2], hlt);
	//hlt_w = 1.0;
	//ISR = 1.0;
	//Nominal
	TLorentzVector j1;
	TLorentzVector j2;
	j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
	j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
	double Dphi = j1.DeltaPhi(j2);

	TLorentzVector j1m;//monojet version
	TLorentzVector j2m;//monojet version
	j1.SetPtEtaPhiE(Jet_PT[0], Jet_Eta[0], Jet_Phi[0], Jet_PT[0]*cosh(Jet_Eta[0]));//Hemisphere
	j2.SetPtEtaPhiE(Jet_PT[1], Jet_Eta[1], Jet_Phi[1], Jet_PT[1]*cosh(Jet_Eta[1]));//Hemisphere
	double DphiM = j1m.DeltaPhi(j2m);//monojet dPHI

	double MET = sqrt(metX[2]*metX[2] + metY[2]*metY[2]);
	if(Jet_PT[0] > 110.0 && N_Jets == 2 && MET > 550.0 && btag == 0 && box == 0 && fabs(DphiM) < 2.5 && fabs(Jet_Eta[0]) < 2.4){
	  NP_ISR_MonoJetONLY += hlt_w*pu_w*ISR;
	}
	
	/*
	dPHI[xs_counter]->Fill(Dphi,hlt_w*scaleF*pu_w*ISR);
	mR_rSQ[xs_counter]->Fill(mr[2], rsq[2], hlt_w*scaleF*pu_w*ISR);
	N_j[xs_counter]->Fill(N_Jets,hlt_w*scaleF*pu_w*ISR);
	j1_Pt[xs_counter]->Fill(Jet_PT[0],hlt_w*scaleF*pu_w*ISR);
	j2_Pt[xs_counter]->Fill(Jet_PT[1],hlt_w*scaleF*pu_w*ISR);
	mR[xs_counter]->Fill(mr[2],hlt_w*scaleF*pu_w*ISR);
	rSQ[xs_counter]->Fill(rsq[2],hlt_w*scaleF*pu_w*ISR);
	mET[xs_counter]->Fill(MET,hlt_w*scaleF*pu_w*ISR);
	*/
	
	if(mr[2] >= 200.0 && rsq[2] >= 0.5 && btag == 0 && box == 0 && fabs(Dphi) < 2.5){
	  if(!(Jet_PT[0] > 110.0 && N_Jets == 2 && MET > 550.0 && fabs(Jet_Eta[0]) < 2.4 && fabs(DphiM) < 2.5))continue;
	  dPHI[xs_counter]->Fill(Dphi,hlt_w*scaleF*pu_w*ISR);
	  mR_rSQ[xs_counter]->Fill(mr[2], rsq[2], hlt_w*scaleF*pu_w*ISR);
	  N_j[xs_counter]->Fill(N_Jets,hlt_w*scaleF*pu_w*ISR);
	  j1_Pt[xs_counter]->Fill(Jet_PT[0],hlt_w*scaleF*pu_w*ISR);
	  j2_Pt[xs_counter]->Fill(Jet_PT[1],hlt_w*scaleF*pu_w*ISR);
	  mR[xs_counter]->Fill(mr[2],hlt_w*scaleF*pu_w*ISR);
	  rSQ[xs_counter]->Fill(rsq[2],hlt_w*scaleF*pu_w*ISR);
	  mET[xs_counter]->Fill(MET,hlt_w*scaleF*pu_w*ISR);
	  
	  if(mr[2] > 200.0 && mr[2] <= 300.0 ){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    h_rsq[xs_counter][0]->Fill(rsq[2], hlt_w*scaleF*pu_w*ISR);
	  }else if(mr[2] > 300.0 && mr[2] <= 400.0){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    h_rsq[xs_counter][1]->Fill(rsq[2], hlt_w*scaleF*pu_w*ISR);
	  }else if(mr[2] > 400.0 && mr[2] <= 600.0){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    h_rsq[xs_counter][2]->Fill(rsq[2], hlt_w*scaleF*pu_w*ISR);
	  }else if(mr[2] > 600.0 && mr[2] <= 3500.0){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    h_rsq[xs_counter][3]->Fill(rsq[2], hlt_w*scaleF*pu_w*ISR);
	  }
	  N_passed += hlt_w*pu_w;
	  N_passed_ISR += hlt_w*pu_w*ISR;

	  if(Jet_PT[0] > 110.0 && N_Jets == 2 && MET > 550.0 && fabs(Jet_Eta[0]) < 2.4 && fabs(DphiM) < 2.5){
	    N_passed_ISR_MonoJet += hlt_w*pu_w*ISR;
	  }
	}
	
      }
      
      double sample_eff_old = N_passed/Gen_Evts;
      double sample_eff = N_passed_ISR/Gen_Evts_isr;
      double s_eff_Mono = N_passed_ISR_MonoJet/Gen_Evts_isr;
            
      std::cout << "Sample Eff OLD: " << sample_eff_old*100 << "%" << std::endl;
      std::cout << "Sample Eff: " << sample_eff*100 << "%" << std::endl;
      std::cout << "Sample Eff MONOJET: " << s_eff_Mono*100 << "%" << std::endl;
      
      std::cout << "Razor/MONOJET: " << (N_passed_ISR_MonoJet/N_passed_ISR)*100 << "%" << std::endl;
      std::cout << "MONOJET ONLY: " <<  (NP_ISR_MonoJetONLY/Gen_Evts_isr)*100 << "%" << std::endl;
      outfile << dm_sample << " & " << sample_eff*100 << "\\%" << "\\" << "\\" << "\n";
      TString current;
      TString mASS;
      std::cout << "deb -1" << std::endl;
      if(dm_sample.find("AV") != std::string::npos){
	mASS = dm_sample.substr(3, dm_sample.find("AV")-3).c_str();
	current = dm_sample.substr(dm_sample.find("AV"),3).c_str();
      }else{
	mASS = dm_sample.substr(3, dm_sample.find("V")-3).c_str();
	current = dm_sample.substr(dm_sample.find("V"),3).c_str();
      }
      //outfile << dm_sample << " & " << (N_passed_ISR_MonoJet/N_passed_ISR)*100 << "\\%" << "\\" << "\\" << "\n";
      outfile << "\\hline" << std::endl;
      std::cout << "deb" << std::endl;
      TString DIR = "SignalPlots/RazorAndMonojet/"; 
      PlotCosmetics(mR[xs_counter],DIR+dm_sample+"_MR", "M_{R} GeV", "Events", "", mASS, current);
      PlotCosmetics(rSQ[xs_counter],DIR+dm_sample+"_RSQ", "R^{2}", "Events", "", mASS, current);
      PlotCosmetics(mET[xs_counter],DIR+dm_sample+"_MET", "#slash{E}_{T}  GeV", "Events", "", mASS, current);
      PlotCosmetics(j1_Pt[xs_counter],DIR+dm_sample+"_J1PT", "P^{j1}_{T} GeV", "Events", "", mASS, current);
      PlotCosmetics(j1_Pt[xs_counter],DIR+dm_sample+"_J2PT", "P^{j2}_{T} GeV", "Events", "", mASS, current);
      PlotCosmetics(N_j[xs_counter],DIR+dm_sample+"_NJets", "N_{jets}", "Events", "", mASS, current);
      PlotCosmetics(dPHI[xs_counter],DIR+dm_sample+"_dPHI", "#Delta#phi(J_{1}, J_{2})", "Events", "", mASS, current);
      PlotCosmetics2D(mR_rSQ[xs_counter],DIR+dm_sample+"_MR_RSQ", "M_{R} GeV", "R^{2}", "colz", mASS, current);
      xs_counter++;
      delete out;
      delete eff;
    }
    
  }else{
    std::cout << "Unable to open the file" << std::endl;
  }

  TFile* fo = new TFile("TEST_OUTPUT.root", "recreate");
  for(int j = 0; j < 24; j++){
    N_j[j]->Write();
    j1_Pt[j]->Write();
    j2_Pt[j]->Write();
    mR[j]->Write();
    rSQ[j]->Write();
    dPHI[j]->Write();
    mR_rSQ[j]->Write();
    mET[j]->Write();
  }
  
  mfile0.close();
  outfile << "\\end{tabular}\n\\end{center}\n\\label{default}\n\\end{table}\n";
  outfile.close();
  
  return 0;
}
