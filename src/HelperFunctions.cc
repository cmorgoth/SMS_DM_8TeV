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
}
