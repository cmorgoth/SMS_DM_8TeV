/*
This code obatains the acceptance 
including the PDF systematic,
Calculates the systematic error band
for the pdf acceptance
it outpus a ROOT file with a histogram containing
the error band
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
    
  std::map<std::string, pdf_and_ISR> sms_mass_pair;//map to store all mass point pair and total generate events;
  int ncteq66 = 45;
  int nmrst2006nnlo = 31;
  int nnnpdf10100 = 101;
  
  TH2F* h_N_cteq66_isr[ncteq66];
  TH2F* h_N_mrst2006nnlo_isr[nmrst2006nnlo];
  TH2F* h_N_nnpdf10100_isr[nnnpdf10100];
  TString ss_h_aux;//aux TString to get histograms from TFile
  for(int j = 0; j < ncteq66; j++){
    ss_h_aux = TString(Form("h_N_pdf_CTEQ66_isr_%d",j));
    h_N_cteq66_isr[j] = (TH2F*)fin->Get(ss_h_aux);
  }
  for(int j = 0; j < nmrst2006nnlo; j++){
    ss_h_aux = TString(Form("h_N_pdf_MRST2006NNLO_isr_%d",j));
    h_N_mrst2006nnlo_isr[j] = (TH2F*)fin->Get(ss_h_aux);
  }
  for(int j = 0; j < nnnpdf10100; j++){
    ss_h_aux = TString(Form("h_N_pdf_NNPDF10100_isr_%d",j));
    h_N_nnpdf10100_isr[j] = (TH2F*)fin->Get(ss_h_aux);
  }
    
    
  for( auto tmp : sms_aux ){
    std::stringstream tmp_ss;
    tmp_ss.str(tmp.first);
    int mass1, mass2;
    tmp_ss >> mass1 >> mass2;
    pdf_and_ISR aux_struct;
    for(int j = 0; j < ncteq66; j++){
      aux_struct.N_pdf_CTEQ66_isr_Tot[j] = GetGenEvents(mass1, mass2, h_N_cteq66_isr[j]);
    }
    for(int j = 0; j < nmrst2006nnlo; j++){
      aux_struct.N_pdf_MRST2006NNLO_isr_Tot[j] = GetGenEvents(mass1, mass2, h_N_mrst2006nnlo_isr[j]);
    }
    for(int j = 0; j < nnnpdf10100; j++){
      aux_struct.N_pdf_NNPDF10100_isr_Tot[j] = GetGenEvents(mass1, mass2, h_N_nnpdf10100_isr[j]);
    }
    CreateStructHistos(aux_struct, tmp.first);
    sms_mass_pair[tmp.first] = aux_struct;
  }

  
  double mr[4], rsq[4], Jet_PT[20], Jet_Eta[20],
    Jet_Phi[20];
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2;
  int btag, box, N_Jets;
  double pu_w, ISR;
    
  double Npassed_ISR;
  int nCTEQ66, nMRST2006NNLO, nNNPDF10100;
  double wCTEQ66[100], wMRST2006NNLO[100], wNNPDF10100[150];
  double N_CTEQ66_isr[100], N_MRST2006NNLO_isr[100], N_NNPDF10100_isr[150];
      
  //Defining needed variables
  outT->SetBranchStatus("*", 0);
  outT->SetBranchStatus("pu_w", 1);
  outT->SetBranchStatus("ISR", 1);
            
  outT->SetBranchStatus("MR", 1);
  outT->SetBranchStatus("RSQ",1);
  outT->SetBranchStatus("nBtag", 1);
  outT->SetBranchStatus("BOX_NUM",1);
  outT->SetBranchStatus("N_Jets",1);
      
  outT->SetBranchStatus("Jet_PT",1);
  outT->SetBranchStatus("Jet_Phi",1);
  outT->SetBranchStatus("Jet_Eta",1);
  
  outT->SetBranchStatus("pTHem1",1);
  outT->SetBranchStatus("pTHem2",1);
  outT->SetBranchStatus("etaHem1",1);
  outT->SetBranchStatus("etaHem2",1);
  outT->SetBranchStatus("phiHem1",1);
  outT->SetBranchStatus("phiHem2",1);
  
  //SMS info
  outT->SetBranchStatus("mssm", 1);
  outT->SetBranchAddress("mssm", mssm);
  
  //pdf info
  outT->SetBranchStatus("nCTEQ66", 1);
  outT->SetBranchStatus("wCTEQ66", 1);
  outT->SetBranchStatus("nMRST2006NNLO", 1);
  outT->SetBranchStatus("wMRST2006NNLO", 1);
  outT->SetBranchStatus("nNNPDF10100", 1);
  outT->SetBranchStatus("wNNPDF10100", 1);
      

  ///////////////////////////////
  ///////////Addresses///////////
  ///////////////////////////////
  
  outT->SetBranchAddress("pu_w", &pu_w);
  outT->SetBranchAddress("ISR", &ISR);
      
  outT->SetBranchAddress("MR", mr);
  outT->SetBranchAddress("RSQ", rsq);
  outT->SetBranchAddress("nBtag", &btag);
  outT->SetBranchAddress("BOX_NUM", &box);
  outT->SetBranchAddress("N_Jets", &N_Jets);
      
  outT->SetBranchAddress("Jet_PT", Jet_PT);
  outT->SetBranchAddress("Jet_Phi", Jet_Phi);
  outT->SetBranchAddress("Jet_Eta", Jet_Eta);
  
  outT->SetBranchAddress("pTHem1", &pTHem1);
  outT->SetBranchAddress("pTHem2", &pTHem2);
  outT->SetBranchAddress("etaHem1", &etaHem1);
  outT->SetBranchAddress("etaHem2", &etaHem2);
  outT->SetBranchAddress("phiHem1", &phiHem1);
  outT->SetBranchAddress("phiHem2", &phiHem2);
    
  //pdf info      
  outT->SetBranchAddress("nCTEQ66", &nCTEQ66);
  outT->SetBranchAddress("wCTEQ66", wCTEQ66);
  outT->SetBranchAddress("nMRST2006NNLO", &nMRST2006NNLO);
  outT->SetBranchAddress("wMRST2006NNLO", wMRST2006NNLO);
  outT->SetBranchAddress("nNNPDF10100", &nNNPDF10100);
  outT->SetBranchAddress("wNNPDF10100", wNNPDF10100);
  
  int N_out = outT->GetEntries();
  
  for( auto tmp : sms_mass_pair ){
    std::cout << "NCteq: " << tmp.second.N_pdf_CTEQ66_isr_Tot[0] << std::endl;
  }


  TCanvas* ca = new TCanvas("c","c", 640, 640);
  TFile* f2 = new 
    TFile("~/Software/git/DM_Signal/trigger/hlt_eff_HTMHT_0mu_BOX_0bTag_Final.root");
  TEfficiency* hlt = (TEfficiency*)f2->Get("Eff2d");
  
  double N_passed = 0.0;
  double N_passed_ISR = 0.0;
  
  for(int j = 0; j < N_out; j++){
    outT->GetEntry(j);
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
    if(mr[2] >= 200.0 && rsq[2] >= 0.5 && btag == 0 && box == 0 && fabs(Dphi) < 2.5){
      double ff_w = hlt_w*pu_w*ISR;
      if(mr[2] > 200.0 && mr[2] <= 300.0 ){
	if(rsq[2] > 1.2)rsq[2] = 1.1;
	for(int l = 0; l < nCTEQ66; l++){
	  sms_mass_pair[s_aux].h_acc_pdf_CTEQ[l][0]->Fill(rsq[2], ff_w*wCTEQ66[l]/sms_mass_pair[s_aux].N_pdf_CTEQ66_isr_Tot[l]);
	}
	for(int l = 0; l < nMRST2006NNLO; l++){
	  sms_mass_pair[s_aux].h_acc_pdf_MRST[l][0]->Fill(rsq[2], ff_w*wMRST2006NNLO[l]/sms_mass_pair[s_aux].N_pdf_MRST2006NNLO_isr_Tot[l]);
	}
	for(int l = 0; l < nNNPDF10100; l++){
	  sms_mass_pair[s_aux].h_acc_pdf_NNPDF[l][0]->Fill(rsq[2], ff_w*wNNPDF10100[l]/sms_mass_pair[s_aux].N_pdf_NNPDF10100_isr_Tot[l]);
	}
      }else if(mr[2] > 300.0 && mr[2] <= 400.0){
	if(rsq[2] > 1.2)rsq[2] = 1.1;
	for(int l = 0; l < nCTEQ66; l++){
	  sms_mass_pair[s_aux].h_acc_pdf_CTEQ[l][1]->Fill(rsq[2], ff_w*wCTEQ66[l]/sms_mass_pair[s_aux].N_pdf_CTEQ66_isr_Tot[l]);
	}
	for(int l = 0; l < nMRST2006NNLO; l++){
	  sms_mass_pair[s_aux].h_acc_pdf_MRST[l][1]->Fill(rsq[2], ff_w*wMRST2006NNLO[l]/sms_mass_pair[s_aux].N_pdf_MRST2006NNLO_isr_Tot[l]);
	}
	for(int l = 0; l < nNNPDF10100; l++){
	  sms_mass_pair[s_aux].h_acc_pdf_NNPDF[l][1]->Fill(rsq[2], ff_w*wNNPDF10100[l]/sms_mass_pair[s_aux].N_pdf_NNPDF10100_isr_Tot[l]);
	}
      }else if(mr[2] > 400.0 && mr[2] <= 600.0){
	if(rsq[2] > 1.2)rsq[2] = 1.1;
	for(int l = 0; l < nCTEQ66; l++){
	  sms_mass_pair[s_aux].h_acc_pdf_CTEQ[l][2]->Fill(rsq[2], ff_w*wCTEQ66[l]/sms_mass_pair[s_aux].N_pdf_CTEQ66_isr_Tot[l]);
	}
	for(int l = 0; l < nMRST2006NNLO; l++){
	  sms_mass_pair[s_aux].h_acc_pdf_MRST[l][2]->Fill(rsq[2], ff_w*wMRST2006NNLO[l]/sms_mass_pair[s_aux].N_pdf_MRST2006NNLO_isr_Tot[l]);
	}
	for(int l = 0; l < nNNPDF10100; l++){
	  sms_mass_pair[s_aux].h_acc_pdf_NNPDF[l][2]->Fill(rsq[2], ff_w*wNNPDF10100[l]/sms_mass_pair[s_aux].N_pdf_NNPDF10100_isr_Tot[l]);
	}
      }else if(mr[2] > 600.0 && mr[2] <= 3500.0){
	if(rsq[2] > 1.2)rsq[2] = 1.1;
	for(int l = 0; l < nCTEQ66; l++){
	  sms_mass_pair[s_aux].h_acc_pdf_CTEQ[l][3]->Fill(rsq[2], ff_w*wCTEQ66[l]/sms_mass_pair[s_aux].N_pdf_CTEQ66_isr_Tot[l]);
	}
	for(int l = 0; l < nMRST2006NNLO; l++){
	  sms_mass_pair[s_aux].h_acc_pdf_MRST[l][3]->Fill(rsq[2], ff_w*wMRST2006NNLO[l]/sms_mass_pair[s_aux].N_pdf_MRST2006NNLO_isr_Tot[l]);
	}
	for(int l = 0; l < nNNPDF10100; l++){
	  sms_mass_pair[s_aux].h_acc_pdf_NNPDF[l][3]->Fill(rsq[2], ff_w*wNNPDF10100[l]/sms_mass_pair[s_aux].N_pdf_NNPDF10100_isr_Tot[l]);
	}
      }
      N_passed += hlt_w*pu_w;
      N_passed_ISR += hlt_w*pu_w*ISR;
    }
    
  }
     
  for( auto tmp : sms_mass_pair){
    std::stringstream ss_aux;
    std::stringstream ss_aux2;
    ss_aux.str(tmp.first);
    int mass1, mass2;
    ss_aux >> mass1 >> mass2;
    ss_aux.str("");
    ss_aux2 << mass1 << "_" << mass2;
    TString sms_s = ss_aux2.str().c_str();
    //CTEQ66 Error Calculation
    int npair = 22;
    for(int cat = 0; cat < 4; cat++){
      for(int bin = 1; bin <= tmp.second.h_acc_pdf_CTEQ[0][cat]->GetNbinsX(); bin++){
	double wplus = 0.0;
	double wminus = 0.0;
	for(int l = 0; l < npair; l++){
	  double wa = 0.0;
	  if(tmp.second.h_acc_pdf_CTEQ[2*l+1][cat]->GetBinContent(bin)){
	    wa = tmp.second.h_acc_pdf_CTEQ[2*l+1][cat]->GetBinContent(bin)/tmp.second.h_acc_pdf_CTEQ[0][cat]->GetBinContent(bin) - 1.0;
	  }
	  double wb = 0.0;
	  if(tmp.second.h_acc_pdf_CTEQ[2*l+2][cat]->GetBinContent(bin)){
	    wb = tmp.second.h_acc_pdf_CTEQ[2*l+2][cat]->GetBinContent(bin)/tmp.second.h_acc_pdf_CTEQ[0][cat]->GetBinContent(bin) - 1.0;
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
	/*
	  std::cout << "CAT: " << cat+1 << " Bin: " << bin << " DeltaUp: " << wplus*100.0 << " % " <<
	  " DeltaDown: " << wminus*100.0 << " %"<< std::endl;
	*/
	double acc_up = tmp.second.h_acc_pdf_CTEQ[0][cat]->GetBinContent(bin)+wplus*tmp.second.h_acc_pdf_CTEQ[0][cat]->GetBinContent(bin);
	double acc_down = tmp.second.h_acc_pdf_CTEQ[0][cat]->GetBinContent(bin)-wminus*tmp.second.h_acc_pdf_CTEQ[0][cat]->GetBinContent(bin);
	tmp.second.h_acc_pdf_CTEQ_up[cat]->SetBinContent(bin, acc_up);
	tmp.second.h_acc_pdf_CTEQ_down[cat]->SetBinContent(bin, acc_down);
      }
    }
    
    //MRST Error Calculation
    npair = 15;
    for(int cat = 0; cat < 4; cat++){
      for(int bin = 1; bin <= tmp.second.h_acc_pdf_MRST[0][cat]->GetNbinsX(); bin++){
	double wplus = 0.0;
	double wminus = 0.0;
	for(int l = 0; l < npair; l++){
	  double wa = 0.0;
	  if(tmp.second.h_acc_pdf_MRST[2*l+1][cat]->GetBinContent(bin)){
	    wa = tmp.second.h_acc_pdf_MRST[2*l+1][cat]->GetBinContent(bin)/tmp.second.h_acc_pdf_MRST[0][cat]->GetBinContent(bin) - 1.0;
	  }
	  double wb = 0.0;
	  if(tmp.second.h_acc_pdf_MRST[2*l+2][cat]->GetBinContent(bin)){
	    wb = tmp.second.h_acc_pdf_MRST[2*l+2][cat]->GetBinContent(bin)/tmp.second.h_acc_pdf_MRST[0][cat]->GetBinContent(bin) - 1.0;
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
	/*
	  std::cout << "CAT: " << cat+1 << " Bin: " << bin << " DeltaUp: " << wplus*100.0 << " % " <<
	  " DeltaDown: " << wminus*100.0 << " %"<< std::endl;
	*/
	double acc_up = tmp.second.h_acc_pdf_MRST[0][cat]->GetBinContent(bin)+wplus*tmp.second.h_acc_pdf_MRST[0][cat]->GetBinContent(bin);
	double acc_down = tmp.second.h_acc_pdf_MRST[0][cat]->GetBinContent(bin)-wminus*tmp.second.h_acc_pdf_MRST[0][cat]->GetBinContent(bin);
	tmp.second.h_acc_pdf_MRST_up[cat]->SetBinContent(bin, acc_up);
	tmp.second.h_acc_pdf_MRST_down[cat]->SetBinContent(bin, acc_down);
      }
    }
    
    for(int cat = 0; cat < 4; cat++){
      for(int bin = 1; bin <= tmp.second.h_acc_pdf_MRST[0][cat]->GetNbinsX(); bin++){
	double cteqUp = tmp.second.h_acc_pdf_CTEQ_up[cat]->GetBinContent(bin);
	double cteqDown = tmp.second.h_acc_pdf_CTEQ_down[cat]->GetBinContent(bin);
	double mrstUp = tmp.second.h_acc_pdf_MRST_up[cat]->GetBinContent(bin);
	double mrstDown = tmp.second.h_acc_pdf_MRST_down[cat]->GetBinContent(bin);
	
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
	
	double CENTRAL_VAL = (MIN+MAX)/2.0;
	tmp.second.h_acc_CENTRAL[cat]->SetBinContent(bin, CENTRAL_VAL);
	if(CENTRAL_VAL > 0){
	  tmp.second.h_acc_PDF_up[cat]->SetBinContent(bin, (MAX - CENTRAL_VAL)/CENTRAL_VAL);
	  tmp.second.h_acc_PDF_down[cat]->SetBinContent(bin, (CENTRAL_VAL - MIN)/CENTRAL_VAL);
	}else{
	  tmp.second.h_acc_PDF_up[cat]->SetBinContent(bin, 0.0);
	  tmp.second.h_acc_PDF_down[cat]->SetBinContent(bin, 0.0);
	}
      }
    }
      

    
    TFile* fo = new TFile("AccBIS/SMS_"+sms_s+".root", "RECREATE");
    TString dummy;
    for(int i = 0; i < 4; i++){
      dummy = Form("PDF_SYS_cat%d",i+1);
      dummy = sms_s+"_"+dummy;
      tmp.second.h_acc_PDF_up[i]->Write(dummy);
      //tmp.second.h_acc_CENTRAL[i]->Write(dummy+"_central");
    }
  }
  return 0;
}
