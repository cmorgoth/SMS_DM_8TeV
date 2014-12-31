#ifndef STRUCTDEFINITION_HH
#define STRUCTDEFINITION_HH 1
//ROOT INCLUDES
#include "TH1F.h"

struct pdf_and_ISR{

  double N_pdf_CTEQ66_isr_Tot[45];
  double N_pdf_MRST2006NNLO_isr_Tot[31];
  double N_pdf_NNPDF10100_isr_Tot[101];
    
  TH1F* h_acc_pdf_CTEQ[45][4];
  TH1F* h_acc_pdf_MRST[31][4];
  TH1F* h_acc_pdf_NNPDF[101][4];
  
  TH1F* h_acc_pdf_CTEQ_up[4];
  TH1F* h_acc_pdf_MRST_up[4];
  TH1F* h_acc_pdf_NNPDF_up[4];
  
  TH1F* h_acc_pdf_CTEQ_down[4];
  TH1F* h_acc_pdf_MRST_down[4];
  TH1F* h_acc_pdf_NNPDF_down[4];
  
  TH1F* h_acc_CENTRAL[4];
  TH1F* h_acc_PDF_up[4];
  TH1F* h_acc_PDF_down[4];
};

#endif
