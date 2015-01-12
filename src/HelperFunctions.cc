//ROOT INCLUDES
#include "TLegend.h"
//LOCAL INCLUDES
#include "HelperFunctions.hh"
#include "StructDefinition.hh"

double GetGenEvents(double mass1, double mass2, TH2F* h){
  
  if(mass1 < 0.0 || mass2 < 0.0){
    std::cerr << "Error message : Negative Masses!\nTerminating Program." << std::endl;
    exit(1);
  }
  
  return h->GetBinContent(h->FindBin(mass1, mass2));
};

bool CreateStructHistos(pdf_and_ISR& tmp, std::string s){
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
  TString SMS = s.c_str();
  TString sn;
  for(int i = 0; i < 4; i++){
    for(int l = 0; l < 45; l++){
      sn = TString(Form("_acc_pdfCTEQ_cat%d_index%d",i+1,l));
      tmp.h_acc_pdf_CTEQ[l][i] = new TH1F(SMS+sn, SMS+sn, r2B[i], v.at(i));
    }
    sn = TString(Form("_acc_pdfCTEQ_up_cat%d",i+1));
    tmp.h_acc_pdf_CTEQ_up[i] = new TH1F(SMS+sn, SMS+sn, r2B[i], v.at(i));
    sn = TString(Form("_acc_pdfCTEQ_down_cat%d",i+1));
    tmp.h_acc_pdf_CTEQ_down[i] = new TH1F(SMS+sn, SMS+sn, r2B[i], v.at(i));
    
    for(int l = 0; l < 31; l++){
      sn = TString(Form("_acc_pdfMRST_cat%d_index%d",i+1, l));
      tmp.h_acc_pdf_MRST[l][i] = new TH1F(SMS+sn, SMS+sn, r2B[i], v.at(i));
    }
    sn = TString(Form("_acc_pdfMRST_up_cat%d",i+1));
    tmp.h_acc_pdf_MRST_up[i] = new TH1F(SMS+sn, SMS+sn, r2B[i], v.at(i));
    sn = TString(Form("_acc_pdfMRST_down_cat%d",i+1));
    tmp.h_acc_pdf_MRST_down[i] = new TH1F(SMS+sn, SMS+sn, r2B[i], v.at(i));
    
    for(int l = 0; l < 101; l++){
      sn = TString(Form("_acc_pdfNNPDF_cat%d_index%d",i+1, l));
      tmp.h_acc_pdf_NNPDF[l][i] = new TH1F(SMS+sn, SMS+sn, r2B[i], v.at(i));
    }
    sn = TString(Form("_acc_pdfNNPDF_up_cat%d",i+1));
    tmp.h_acc_pdf_NNPDF_up[i] = new TH1F(SMS+sn, SMS+sn, r2B[i], v.at(i));
    sn = TString(Form("_acc_pdfNNPDF_down_cat%d",i+1));
    tmp.h_acc_pdf_NNPDF_down[i] = new TH1F(SMS+sn, SMS+sn, r2B[i], v.at(i));
    
    sn = TString(Form("_acc_CENTRAL_cat%d",i+1));
    tmp.h_acc_CENTRAL[i] = new TH1F(SMS+sn, SMS+sn, r2B[i], v.at(i));
    sn = TString(Form("_acc_TOTAL_PDF_up_cat%d",i+1));
    tmp.h_acc_PDF_up[i] = new TH1F(SMS+sn, SMS+sn, r2B[i], v.at(i));
    sn = TString(Form("_acc_TOTAL_PDF_down_cat%d",i+1));
    tmp.h_acc_PDF_down[i] = new TH1F(SMS+sn, SMS+sn, r2B[i], v.at(i));  
  }
  
  return true;
};

bool CreateStructHistos(SignalStruct& tmp, std::string s){
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
  
  TString s_aux = s.c_str();
  TString sn;
  int nb_met = 20;
  int nb_njets = 9;
  float met_xlow = 0.0;
  float met_xhigh = 1500.0;
  float njets_xlow = 1.0;
  float njets_xhigh = 10.0;
  //SYSTEMATIC HISTOS
  for(int i = 0; i < 4; i++){
      sn = TString(Form("signal_cat%d_met",i+1));
      tmp.h_met[i] = new TH1F(sn+"_"+s_aux, sn+"_"+s_aux, nb_met, met_xlow, met_xhigh);
      sn = TString(Form("signal_cat%d_njets",i+1));
      tmp.h_njets[i] = new TH1F(sn+"_"+s_aux, sn+"_"+s_aux, nb_njets, njets_xlow, njets_xhigh);
      
      sn = TString(Form("signal_cat%d",i));
      tmp.h_rsq[i] = new TH1F(sn+"_"+s_aux, sn+"_"+s_aux, r2B[i], v.at(i));
      
      tmp.h_rsq_ISR_up[i] = new TH1F(sn+"_"+s_aux+"_ISR_up", sn+"_"+s_aux+"ISR_up", r2B[i], v.at(i));
      tmp.h_rsq_ISR_down[i] = new TH1F(sn+"_"+s_aux+"_ISR_down", sn+"_"+s_aux+"ISR_down", r2B[i], v.at(i));
      
      tmp.h_rsq_JES_up[i] = new TH1F(sn+"_"+s_aux+"_JES_up", sn+"_"+s_aux+"JES_up", r2B[i], v.at(i));
      tmp.h_rsq_JES_down[i] = new TH1F(sn+"_"+s_aux+"_JES_down", sn+"_"+s_aux+"JES_down", r2B[i], v.at(i));
      
      tmp.s_up[i] = new TH1F(sn+"_"+s_aux+"_up", sn+"_"+s_aux+"_up", r2B[i], v.at(i));
      tmp.s_down[i] = new TH1F(sn+"_"+s_aux+"_down", sn+"_"+s_aux+"_down", r2B[i], v.at(i));
  }
  
  return true;
};

bool GetMasses(std::string aux, float& mass1, float& mass2){
  int length = aux.size();
  int first = aux.find_first_of("_");
  int second = aux.find_last_of("_");
  
  std::string s1 = aux.substr(0, first);
  std::string s2 = aux.substr(first+1, second-(first+1));
  //std::cout << s1 << " " << s2 << std::endl;
  mass1 = atoi(s1.c_str());
  mass2 = atoi(s2.c_str());
  return true;
};

bool PlotLimit(TH2F* h, TGraph* obs, TGraph* exp, TGraph* mSigma, TGraph* pSigma, TString name){
  TCanvas* C = new TCanvas("C", "C", 640, 640);
  
  obs->SetLineColor(kGreen);
  obs->SetLineWidth(4);
  obs->SetLineStyle(1);
  
  exp->SetLineColor(kRed);
  exp->SetLineWidth(2);
  exp->SetLineStyle(2);
  
  mSigma->SetLineColor(kGreen-3);
  mSigma->SetLineWidth(1);
  mSigma->SetLineStyle(2);

  pSigma->SetLineColor(kGreen-3);
  pSigma->SetLineWidth(1);
  pSigma->SetLineStyle(2);

   obs->GetXaxis()->SetTitle("stop mass [GeV]");
  //obs->GetXaxis()->CenterTitle(kTRUE);
  //obs->GetYaxis()->CenterTitle(kTRUE);
  obs->GetYaxis()->SetTitleOffset(1.3);
  obs->GetYaxis()->SetTitle("LSP mass [GeV]");
  obs->SetTitle("");
  
  obs->Draw("AL");
  exp->Draw("same");
  pSigma->Draw("same");
  mSigma->Draw("same");
  h->Draw("colzsame");
  
  TLegend* leg = new TLegend(0.15, 0.7, 0.55, 0.85);//(xmin, ymin, xmax, ymax)
  leg->AddEntry(obs, "obs","l");
  leg->AddEntry(exp, "Exp.","l");
  leg->AddEntry(pSigma, "#pm 1#sigma","l");
  
  leg->SetTextSize(.022);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->Draw();
     
  C->SaveAs("SMS_LimitsFile/T2ccLimit.pdf");
  C->SaveAs("SMS_LimitsFile/T2ccLimit.C");
  return true;
};
