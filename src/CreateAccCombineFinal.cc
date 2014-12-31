/*
This code obatains the acceptance 
including the PDF systematic,
Calculates the systematic error band
for the pdf acceptance
it outpus a ROOT file with a histogram containing
the error band
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
#include <cmath>

int main(){

  //gROOT->Reset();


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
  TFile* f2 = new 
    TFile("~/Software/git/DM_Signal/trigger/hlt_eff_HTMHT_0mu_BOX_0bTag_Final.root");
  TEfficiency* hlt = (TEfficiency*)f2->Get("Eff2d");

  //Histograms for calculating acceptances
  TH1F* h_acc_pdf_CTEQ[24][45][4];
  TH1F* h_acc_pdf_MRST[24][31][4];
  TH1F* h_acc_pdf_NNPDF[24][101][4];
  
  TH1F* h_acc_pdf_CTEQ_up[24][4];
  TH1F* h_acc_pdf_MRST_up[24][4];
  TH1F* h_acc_pdf_NNPDF_up[24][4];
  
  TH1F* h_acc_pdf_CTEQ_down[24][4];
  TH1F* h_acc_pdf_MRST_down[24][4];
  TH1F* h_acc_pdf_NNPDF_down[24][4];

  TH1F* h_acc_CENTRAL[24][4];
  TH1F* h_acc_PDF_up[24][4];
  TH1F* h_acc_PDF_down[24][4];
  
  std::ifstream mfile0("list_DM_BIS_PDF.list");
  std::ofstream outfile("eff_table_normal_R2_0p5_MR_200_Dphi_B_2p5_New.tex");

  std::string fname0;
  std::cout.precision(16);
  int xs_counter = 0;
   
  if (mfile0.is_open()){
    while ( mfile0.good() ){
      mfile0 >> fname0;
      if(mfile0.eof())break;
      std::cout << fname0 << std::endl;
      int low_ = fname0.find("DMm");
      int high_ = fname0.find("_testMC_0.root") - low_;
      
      std::string dm_sample = fname0.substr(low_,high_);
   
      TFile* f = new TFile(fname0.c_str());
      TTree* eff = (TTree*)f->Get("effTree");
      TTree* out = (TTree*)f->Get("outTree");
      
      double mr[4], rsq[4], Jet_PT[20], Jet_Eta[20],
	Jet_Phi[20];
      double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2;
      int btag, box, N_Jets;
      double pu_w, ISR;
      
      double Npassed_ISR;
      int nCTEQ66, nMRST2006NNLO, nNNPDF10100;
      double wCTEQ66[100], wMRST2006NNLO[100], wNNPDF10100[150];
      double N_CTEQ66_isr[100], N_MRST2006NNLO_isr[100], N_NNPDF10100_isr[150];
      
      //Defining Branches to be used.
      eff->SetBranchStatus("*", 0);
      eff->SetBranchStatus("Npassed_ISR", 1);
      
      eff->SetBranchStatus("nCTEQ66", 1);
      eff->SetBranchStatus("N_pdf_CTEQ66_isr", 1);
      
      eff->SetBranchStatus("nMRST2006NNLO", 1);
      eff->SetBranchStatus("N_pdf_MRST2006NNLO_isr", 1);
      
      eff->SetBranchStatus("nNNPDF10100", 1);
      eff->SetBranchStatus("N_pdf_NNPDF10100_isr", 1);
      
      //Addresses 
      eff->SetBranchAddress("Npassed_ISR", &Npassed_ISR);
      
      eff->SetBranchAddress("nCTEQ66", &nCTEQ66);
      eff->SetBranchAddress("N_pdf_CTEQ66_isr", N_CTEQ66_isr);
      
      eff->SetBranchAddress("nMRST2006NNLO", &nMRST2006NNLO);
      eff->SetBranchAddress("N_pdf_MRST2006NNLO_isr", N_MRST2006NNLO_isr); 
      
      eff->SetBranchAddress("nNNPDF10100", &nNNPDF10100);
      eff->SetBranchAddress("N_pdf_NNPDF10100_isr", N_NNPDF10100_isr);
      
      int N_eff = eff->GetEntries();
            
      double N_pdf_CTEQ66_isr_Tot[100];
      double N_pdf_MRST2006NNLO_isr_Tot[100];
      double N_pdf_NNPDF10100_isr_Tot[150];
      
      for(int k = 0; k < 100; k++){
	N_pdf_CTEQ66_isr_Tot[k] = 0.0;
	N_pdf_MRST2006NNLO_isr_Tot[k] = 0.0;
      }
      for(int k = 0; k < 150; k++){
	N_pdf_NNPDF10100_isr_Tot[k] = 0.0;
      }
      
      for(int i = 0; i < N_eff; i++){
	eff->GetEntry(i);
	std::cout << "nCTEQ66: " << nCTEQ66 << std::endl; 
	for(int l = 0; l < 45; l++){
	  N_pdf_CTEQ66_isr_Tot[l] += N_CTEQ66_isr[l];
	}
	for(int l = 0; l < 31; l++){
	  N_pdf_MRST2006NNLO_isr_Tot[l] += N_MRST2006NNLO_isr[l];
	}
	for(int l = 0; l < 101; l++){
	  N_pdf_NNPDF10100_isr_Tot[l] += N_NNPDF10100_isr[l];
	}
      }
      
      TString DM_Name = dm_sample.c_str();
      TString sn;      
      for(int i = 0; i < 4; i++){
	for(int l = 0; l < 45; l++){
	  sn = TString(Form("_acc_pdfCTEQ_cat%d_index%d",i+1,l));
	  h_acc_pdf_CTEQ[xs_counter][l][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	}
	sn = TString(Form("_acc_pdfCTEQ_up_cat%d",i+1));
	h_acc_pdf_CTEQ_up[xs_counter][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	sn = TString(Form("_acc_pdfCTEQ_down_cat%d",i+1));
	h_acc_pdf_CTEQ_down[xs_counter][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	
	for(int l = 0; l < 31; l++){
	  sn = TString(Form("_acc_pdfMRST_cat%d_index%d",i+1, l));
	  h_acc_pdf_MRST[xs_counter][l][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	}
	sn = TString(Form("_acc_pdfMRST_up_cat%d",i+1));
	h_acc_pdf_MRST_up[xs_counter][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	sn = TString(Form("_acc_pdfMRST_down_cat%d",i+1));
	h_acc_pdf_MRST_down[xs_counter][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	
	for(int l = 0; l < 101; l++){
	  sn = TString(Form("_acc_pdfNNPDF_cat%d_index%d",i+1, l));
	  h_acc_pdf_NNPDF[xs_counter][l][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	}
	sn = TString(Form("_acc_pdfNNPDF_up_cat%d",i+1));
	h_acc_pdf_NNPDF_up[xs_counter][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	sn = TString(Form("_acc_pdfNNPDF_down_cat%d",i+1));
	h_acc_pdf_NNPDF_down[xs_counter][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));

	sn = TString(Form("_acc_CENTRAL_cat%d",i+1));
	h_acc_CENTRAL[xs_counter][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	sn = TString(Form("_acc_TOTAL_PDF_up_cat%d",i+1));
	h_acc_PDF_up[xs_counter][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	sn = TString(Form("_acc_TOTAL_PDF_down_cat%d",i+1));
	h_acc_PDF_down[xs_counter][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	
      }
      
      //Defining needed variables
      out->SetBranchStatus("*", 0);
      out->SetBranchStatus("pu_w", 1);
      out->SetBranchStatus("ISR", 1);
            
      out->SetBranchStatus("MR", 1);
      out->SetBranchStatus("RSQ",1);
      out->SetBranchStatus("nBtag", 1);
      out->SetBranchStatus("BOX_NUM",1);
      out->SetBranchStatus("N_Jets",1);
      
      out->SetBranchStatus("Jet_PT",1);
      out->SetBranchStatus("Jet_Phi",1);
      out->SetBranchStatus("Jet_Eta",1);
            
      out->SetBranchStatus("pTHem1",1);
      out->SetBranchStatus("pTHem2",1);
      out->SetBranchStatus("etaHem1",1);
      out->SetBranchStatus("etaHem2",1);
      out->SetBranchStatus("phiHem1",1);
      out->SetBranchStatus("phiHem2",1);
      
      //pdf info
      out->SetBranchStatus("nCTEQ66", 1);
      out->SetBranchStatus("wCTEQ66", 1);
      out->SetBranchStatus("nMRST2006NNLO", 1);
      out->SetBranchStatus("wMRST2006NNLO", 1);
      out->SetBranchStatus("nNNPDF10100", 1);
      out->SetBranchStatus("wNNPDF10100", 1);
      

      ///////////////////////////////
      ///////////Addresses///////////
      ///////////////////////////////
      
      out->SetBranchAddress("pu_w", &pu_w);
      out->SetBranchAddress("ISR", &ISR);
      
      out->SetBranchAddress("MR", mr);
      out->SetBranchAddress("RSQ", rsq);
      out->SetBranchAddress("nBtag", &btag);
      out->SetBranchAddress("BOX_NUM", &box);
      out->SetBranchAddress("N_Jets", &N_Jets);
      
      out->SetBranchAddress("Jet_PT", Jet_PT);
      out->SetBranchAddress("Jet_Phi", Jet_Phi);
      out->SetBranchAddress("Jet_Eta", Jet_Eta);
            
      out->SetBranchAddress("pTHem1", &pTHem1);
      out->SetBranchAddress("pTHem2", &pTHem2);
      out->SetBranchAddress("etaHem1", &etaHem1);
      out->SetBranchAddress("etaHem2", &etaHem2);
      out->SetBranchAddress("phiHem1", &phiHem1);
      out->SetBranchAddress("phiHem2", &phiHem2);
      
      //pdf info      
      out->SetBranchAddress("nCTEQ66", &nCTEQ66);
      out->SetBranchAddress("wCTEQ66", wCTEQ66);
      out->SetBranchAddress("nMRST2006NNLO", &nMRST2006NNLO);
      out->SetBranchAddress("wMRST2006NNLO", wMRST2006NNLO);
      out->SetBranchAddress("nNNPDF10100", &nNNPDF10100);
      out->SetBranchAddress("wNNPDF10100", wNNPDF10100);
      
      int N_out = out->GetEntries();
            
      double N_passed = 0.0;
      double N_passed_ISR = 0.0;
      
      for(int j = 0; j < N_out; j++){
	out->GetEntry(j);
	double hlt_w = HLTscale(mr[2], rsq[2], hlt);
	//Nominal
	TLorentzVector j1;
	TLorentzVector j2;
	j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
	j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
	double Dphi = j1.DeltaPhi(j2);
	if(mr[2] >= 200.0 && rsq[2] >= 0.5 && btag == 0 && box == 0 && fabs(Dphi) < 2.5){
	  
	  double ff_w = hlt_w*pu_w*ISR;
	  if(mr[2] > 200.0 && mr[2] <= 300.0 ){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    for(int l = 0; l < nCTEQ66; l++){
	      h_acc_pdf_CTEQ[xs_counter][l][0]->Fill(rsq[2], ff_w*wCTEQ66[l]/N_pdf_CTEQ66_isr_Tot[l]);
	    }
	    for(int l = 0; l < nMRST2006NNLO; l++){
	      h_acc_pdf_MRST[xs_counter][l][0]->Fill(rsq[2], ff_w*wMRST2006NNLO[l]/N_pdf_MRST2006NNLO_isr_Tot[l]);
	    }
	    for(int l = 0; l < nNNPDF10100; l++){
	      h_acc_pdf_NNPDF[xs_counter][l][0]->Fill(rsq[2], ff_w*wNNPDF10100[l]/N_pdf_NNPDF10100_isr_Tot[l]);
	    }
	  }else if(mr[2] > 300.0 && mr[2] <= 400.0){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    for(int l = 0; l < nCTEQ66; l++){
	      h_acc_pdf_CTEQ[xs_counter][l][1]->Fill(rsq[2], ff_w*wCTEQ66[l]/N_pdf_CTEQ66_isr_Tot[l]);
	    }
	    for(int l = 0; l < nMRST2006NNLO; l++){
	      h_acc_pdf_MRST[xs_counter][l][1]->Fill(rsq[2], ff_w*wMRST2006NNLO[l]/N_pdf_MRST2006NNLO_isr_Tot[l]);
	    }
	    for(int l = 0; l < nNNPDF10100; l++){
	      h_acc_pdf_NNPDF[xs_counter][l][1]->Fill(rsq[2], ff_w*wNNPDF10100[l]/N_pdf_NNPDF10100_isr_Tot[l]);
	    }
	  }else if(mr[2] > 400.0 && mr[2] <= 600.0){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    for(int l = 0; l < nCTEQ66; l++){
	      h_acc_pdf_CTEQ[xs_counter][l][2]->Fill(rsq[2], ff_w*wCTEQ66[l]/N_pdf_CTEQ66_isr_Tot[l]);
	    }
	    for(int l = 0; l < nMRST2006NNLO; l++){
	      h_acc_pdf_MRST[xs_counter][l][2]->Fill(rsq[2], ff_w*wMRST2006NNLO[l]/N_pdf_MRST2006NNLO_isr_Tot[l]);
	    }
	    for(int l = 0; l < nNNPDF10100; l++){
	      h_acc_pdf_NNPDF[xs_counter][l][2]->Fill(rsq[2], ff_w*wNNPDF10100[l]/N_pdf_NNPDF10100_isr_Tot[l]);
	    }
	  }else if(mr[2] > 600.0 && mr[2] <= 3500.0){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    for(int l = 0; l < nCTEQ66; l++){
	      h_acc_pdf_CTEQ[xs_counter][l][3]->Fill(rsq[2], ff_w*wCTEQ66[l]/N_pdf_CTEQ66_isr_Tot[l]);
	    }
	    for(int l = 0; l < nMRST2006NNLO; l++){
	      h_acc_pdf_MRST[xs_counter][l][3]->Fill(rsq[2], ff_w*wMRST2006NNLO[l]/N_pdf_MRST2006NNLO_isr_Tot[l]);
	    }
	    for(int l = 0; l < nNNPDF10100; l++){
	      h_acc_pdf_NNPDF[xs_counter][l][3]->Fill(rsq[2], ff_w*wNNPDF10100[l]/N_pdf_NNPDF10100_isr_Tot[l]);
	    }
	  }
	  N_passed += hlt_w*pu_w;
	  N_passed_ISR += hlt_w*pu_w*ISR;
	}
	
      }
      
      TString data_card_name1 = dm_sample.c_str();
      TString s1 = data_card_name1;
      
      //CTEQ66 Error Calculation
      int npair = 22;
      for(int cat = 0; cat < 4; cat++){
	for(int bin = 1; bin <= h_acc_pdf_CTEQ[xs_counter][0][cat]->GetNbinsX(); bin++){
	  double wplus = 0.0;
	  double wminus = 0.0;
	  for(int l = 0; l < npair; l++){
	    double wa = 0.0;
	    if(h_acc_pdf_CTEQ[xs_counter][2*l+1][cat]->GetBinContent(bin)){
	      wa = h_acc_pdf_CTEQ[xs_counter][2*l+1][cat]->GetBinContent(bin)/h_acc_pdf_CTEQ[xs_counter][0][cat]->GetBinContent(bin) - 1.0;
	    }
	    double wb = 0.0;
	    if(h_acc_pdf_CTEQ[xs_counter][2*l+2][cat]->GetBinContent(bin)){
	      wb = h_acc_pdf_CTEQ[xs_counter][2*l+2][cat]->GetBinContent(bin)/h_acc_pdf_CTEQ[xs_counter][0][cat]->GetBinContent(bin) - 1.0;
	    }
	    
	    if(!std::isfinite(wa) || !std::isfinite(wb))continue;
	    
	    if (wa>wb) {
	      if (wa<0.) wa = 0.;
	      if (wb>0.) wb = 0.;
	      wplus += wa*wa;
	      wminus += wb*wb;
	    } else {
	      if (wb<0.) wb = 0.;
	      if (wa>0.) wa = 0.;
	      wplus += wb*wb;
	      wminus += wa*wa;
	    }
	    
	  }
	  wplus = sqrt(wplus);
	  wminus = sqrt(wminus);
	  std::cout << "CAT: " << cat+1 << " Bin: " << bin << " DeltaUp: " << wplus*100.0 << " % " <<
	    " DeltaDown: " << wminus*100.0 << " %"<< std::endl;
	  double acc_up = h_acc_pdf_CTEQ[xs_counter][0][cat]->GetBinContent(bin)+wplus*h_acc_pdf_CTEQ[xs_counter][0][cat]->GetBinContent(bin);
	  double acc_down = h_acc_pdf_CTEQ[xs_counter][0][cat]->GetBinContent(bin)-wminus*h_acc_pdf_CTEQ[xs_counter][0][cat]->GetBinContent(bin);
	  h_acc_pdf_CTEQ_up[xs_counter][cat]->SetBinContent(bin, acc_up);
	  h_acc_pdf_CTEQ_down[xs_counter][cat]->SetBinContent(bin, acc_down);
	}
      }
      
      //MRST Error Calculation
      npair = 15;
      for(int cat = 0; cat < 4; cat++){
	for(int bin = 1; bin <= h_acc_pdf_MRST[xs_counter][0][cat]->GetNbinsX(); bin++){
	  double wplus = 0.0;
	  double wminus = 0.0;
	  for(int l = 0; l < npair; l++){
	    double wa = 0.0;
	    if(h_acc_pdf_MRST[xs_counter][2*l+1][cat]->GetBinContent(bin)){
	      wa = h_acc_pdf_MRST[xs_counter][2*l+1][cat]->GetBinContent(bin)/h_acc_pdf_MRST[xs_counter][0][cat]->GetBinContent(bin) - 1.0;
	    }
	    double wb = 0.0;
	    if(h_acc_pdf_MRST[xs_counter][2*l+2][cat]->GetBinContent(bin)){
	      wb = h_acc_pdf_MRST[xs_counter][2*l+2][cat]->GetBinContent(bin)/h_acc_pdf_MRST[xs_counter][0][cat]->GetBinContent(bin) - 1.0;
	    }
	    
	    if(!std::isfinite(wa) || !std::isfinite(wb))continue;
	    
	    if (wa>wb) {
	      if (wa<0.) wa = 0.;
	      if (wb>0.) wb = 0.;
	      wplus += wa*wa;
	      wminus += wb*wb;
	    } else {
	      if (wb<0.) wb = 0.;
	      if (wa>0.) wa = 0.;
	      wplus += wb*wb;
	      wminus += wa*wa;
	    }
	    
	  }
	  wplus = sqrt(wplus);
	  wminus = sqrt(wminus);
	  std::cout << "CAT: " << cat+1 << " Bin: " << bin << " DeltaUp: " << wplus*100.0 << " % " <<
	    " DeltaDown: " << wminus*100.0 << " %"<< std::endl;
	  double acc_up = h_acc_pdf_MRST[xs_counter][0][cat]->GetBinContent(bin)+wplus*h_acc_pdf_MRST[xs_counter][0][cat]->GetBinContent(bin);
	  double acc_down = h_acc_pdf_MRST[xs_counter][0][cat]->GetBinContent(bin)-wminus*h_acc_pdf_MRST[xs_counter][0][cat]->GetBinContent(bin);
	  h_acc_pdf_MRST_up[xs_counter][cat]->SetBinContent(bin, acc_up);
	  h_acc_pdf_MRST_down[xs_counter][cat]->SetBinContent(bin, acc_down);
	}
      }
      
      for(int cat = 0; cat < 4; cat++){
	for(int bin = 1; bin <= h_acc_pdf_MRST[xs_counter][0][cat]->GetNbinsX(); bin++){
	  double cteqUp = h_acc_pdf_CTEQ_up[xs_counter][cat]->GetBinContent(bin);
	  double cteqDown = h_acc_pdf_CTEQ_down[xs_counter][cat]->GetBinContent(bin);
	  double mrstUp = h_acc_pdf_MRST_up[xs_counter][cat]->GetBinContent(bin);
	  double mrstDown = h_acc_pdf_MRST_down[xs_counter][cat]->GetBinContent(bin);
	  
	  double MAX = -99.0;
	  double MIN = -99.0;
	  if(cteqUp > mrstUp){
	    MAX = cteqUp;
	  }else{
	    MAX = mrstUp;
	  }
	  
	  if(cteqUp < mrstUp){
	    MIN = cteqDown;
	  }else{
	    MIN = mrstDown;
	  }
	  
	  h_acc_CENTRAL[xs_counter][cat]->SetBinContent(bin, (MIN+MAX)/2.0);
	  h_acc_PDF_up[xs_counter][cat]->SetBinContent(bin, (MAX - (MIN+MAX)/2.0)/((MIN+MAX)/2.0));
	  h_acc_PDF_down[xs_counter][cat]->SetBinContent(bin, ((MIN+MAX)/2.0 - MIN)/((MIN+MAX)/2.0));
	}
      }
      

      TFile* fo = new TFile("AccBIS/"+s1+"_Acc.root", "RECREATE");
      TString dummy;
      for(int i = 0; i < 4; i++){
	dummy = Form("PDF_SYS_cat%d",i+1);
	h_acc_PDF_up[xs_counter][i]->Write(dummy);
      }
    }
    std::cout << " deb 1" << std::endl;
    xs_counter++;
  }else{
    std::cout << "Unable to open the file" << std::endl;
  }
  std::cout << " deb 2" << std::endl;
  mfile0.close();
  std::cout << " deb 2.1" << std::endl;
  return 0;
  std::cout << " deb 2.2" << std::endl;
}
