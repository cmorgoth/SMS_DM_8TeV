#include "DM_1DRatio.hh"

int RatioPlots(TH1F* h1, TH1F* h2, TString h1Name = "h1Name", TString h2Name = "h2Name", TString fname = "default_name", TString type = "defaultType"){
	
  TCanvas* C = new TCanvas("C", "C	", 400, 500);
  C->cd();  
  TLegend* leg;
  
  TH1F*  RATIO;
  TString label;
  if(type == "MR"){
    RATIO = new TH1F("RATIO", fname + "_" + type , BaseDM::MR_Bins, BaseDM::MR_BinArr);
    label = "M_{R}";
    //h1->GetXaxis()->SetRangeUser(200,3500);
    //h2->GetXaxis()->SetRangeUser(200,3500);
    //RATIO->GetXaxis()->SetRangeUser(200,3500);

  }else if(type == "RSQ" ){
    RATIO = new TH1F("RATIO", fname + "_" + type , BaseDM::RSQ_Bins, BaseDM::RSQ_BinArr);
    label = "R^{2}";
    //h1->GetXaxis()->SetRangeUser(0.5, 2.5);
    //h2->GetXaxis()->SetRangeUser(0.5, 2.5);
    //RATIO->GetXaxis()->SetRangeUser(0.5, 2.5);

  }else{
    delete RATIO;
    delete C;
    std::cout << "Unknown Type, please use: MR or RSQ" << std::endl;
    return -1;
  }
  
  RATIO->Divide(h1, h2, 1, 1, "");
  RATIO->GetYaxis()->SetRangeUser(.0, 2.05);

  h1->SetLineColor(2);
  h1->SetMarkerSize(1.);
  h1->SetMarkerColor(2);
  h1->SetMarkerStyle(20);
  h1->SetFillColor(2);
  
  h2->SetLineColor(4);
  h2->SetMarkerSize(1.);
  h2->SetMarkerColor(4);
  h2->SetMarkerStyle(20);
  h2->SetFillColor(4);
  
  h1->SetStats(0);
  h1->SetTitle("");
  h2->SetTitle("");
  h2->SetStats(0);
  h2->SetXTitle( type );
  h1->SetXTitle( type );
  
  std::cout << "GET X Title Size:  " << RATIO->GetYaxis()->GetTitleSize() << std::endl;
   
  TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();
  h1->DrawCopy();
  h2->Draw("same");
  C->cd();
  
  leg = new TLegend(0.55, 0.7, 0.89, 0.9);//(xmin, ymin, xmax, ymax)
  leg->AddEntry(h1, label + " " + h1Name ,"f");
  leg->AddEntry(h2, label + " " + h2Name ,"f");
  leg->SetTextSize(.022);
  leg->SetFillColor(0);
  leg->Draw();
  pad1->SetLogy();
  C->Update();
  
  TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.25);
  pad2->SetTopMargin(0);
  pad2->SetTopMargin(0.008);
  pad2->SetBottomMargin(0.25);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();

  RATIO->SetLineColor(4);
  RATIO->SetStats(0);
  RATIO->SetTitle("");
  RATIO->GetXaxis()->SetLabelSize(0.1);
  RATIO->GetYaxis()->SetLabelSize(0.08);
  RATIO->GetYaxis()->SetTitleOffset(0.35);
  RATIO->GetXaxis()->SetTitleOffset(0.88);
  RATIO->GetYaxis()->SetTitleSize(0.11);
  RATIO->GetXaxis()->SetTitleSize(0.11);
  RATIO->SetXTitle( label );
  RATIO->SetYTitle("Ratio");
  RATIO->SetLineColor(4);
  RATIO->SetMarkerSize(.7);
  RATIO->SetMarkerColor(4);
  RATIO->SetMarkerStyle(20);
  RATIO->SetFillColor(4);

  RATIO->Draw();
  C->cd();
  C->SaveAs(fname + ".pdf");
  C->SaveAs(fname + ".png");
  
  delete leg;
  delete C;
  delete RATIO;
  
  return 0;
  
};

int RatioPlotsBandV2(TH1F* h1, TH1F* h2, TString h1Name = "h1Name", TString h2Name = "h2Name", TString fname = "default_name", TString type = "defaultType", int nbins = 0, float* bins = NULL, int color = 0){
	
  TCanvas* C = new TCanvas("C", "C	", 400, 500);
  C->cd();  
  TLegend* leg;
  
  TH1F*  RATIO;// = new TH1F("RATIO","Data_to_Prediction",BaseDM::MR_Bins, BaseDM::MR_BinArr);
  TH1F*  RATIO2;
  TString label;
  if(type == "MR"){
    RATIO = new TH1F("RATIO", fname + "_" + type , nbins, bins);
    RATIO2 = new TH1F("RATIO2", fname + "_" + type , nbins, bins);
    label = "M_{R}";
    h1->GetXaxis()->SetRangeUser(200,3500);
    h2->GetXaxis()->SetRangeUser(200,3500);
    RATIO->GetXaxis()->SetRangeUser(200,3500);

  }else if(type == "RSQ" ){
    RATIO = new TH1F("RATIO", fname + "_" + type , nbins, bins);
    RATIO2 = new TH1F("RATIO2", fname + "_" + type , nbins, bins);
    label = "R^{2}";
    h1->GetXaxis()->SetRangeUser(0.5, 1.2);
    h2->GetXaxis()->SetRangeUser(0.5, 1.2);
    RATIO->GetXaxis()->SetRangeUser(0.5, 1.2);
  }else{
    delete RATIO;
    delete RATIO2;
    delete C;
    std::cout << "Unknown Type, please use: MR or RSQ" << std::endl;
    return -1;
  }
  
  //RATIO->Divide(h1, h2, 1, 1, "");
  
  for(int i = 1; i <= h1->GetNbinsX(); i++){
    for(int j = 1; j <=  h1->GetNbinsY(); j++){
      double r = .0, r2 = .0;
      double err = .0, err2 = .0;
      
      if(h2->GetBinContent(i,j) !=0 ){
	r = h1->GetBinContent(i,j)/h2->GetBinContent(i,j);
	err = h1->GetBinError(i,j)/h2->GetBinContent(i,j);
	r2 = h2->GetBinContent(i,j)/h2->GetBinContent(i,j);
	err2 = h2->GetBinError(i,j)/h2->GetBinContent(i,j);
      }
      RATIO->SetBinContent(i,j,r);
      RATIO->SetBinError(i,j,err);
      RATIO2->SetBinContent(i,j,r2);
      RATIO2->SetBinError(i,j,err2);
      //std::cout << "Error: " << h2->GetBinError(i,j) << " " << err << std::endl;
      //std::cout << "BinContent: " << h1->GetBinContent(i,j) << std::endl;
    }
  }
  RATIO->GetYaxis()->SetRangeUser(.0, 3.05);
  RATIO2->GetYaxis()->SetRangeUser(.0, 3.05);
  //RATIO->GetYaxis()->SetRangeUser(.7, 1.3);
  //RATIO2->GetYaxis()->SetRangeUser(.7, 1.3);
  
 
  //RATIO->Sumw2();
  //RATIO2->Sumw2();
  h1->SetLineColor(1);//data
  h1->SetMarkerSize(.3);
  h1->SetMarkerColor(1);
  h1->SetMarkerStyle(20);
  h1->SetFillColor(1);
  h1->SetLineWidth(1.7);
  h2->SetLineWidth(2);
  TH1F* h2clone = (TH1F*)h2->Clone("h2clone");
  h2clone->SetFillColor(0);
  switch(color){
  case 0:
    h2->SetLineColor(kGreen-10);
    h2->SetFillColor(kGreen-10);
    h2clone->SetLineColor(kGreen);
    break;
  case 1:
    h2->SetLineColor(kRed-10);
    h2->SetFillColor(kRed-10);
    h2clone->SetLineColor(kRed);
    break;
  case 2:
    h2->SetLineColor(kBlue-10);
    h2->SetFillColor(kBlue-10);
    h2clone->SetLineColor(kBlue);
    break;
  default:
    std::cout << "---Default---" << std::endl;
    break;
  }
  h1->SetStats(0);
  h1->SetTitle("");
  h2->SetTitle("");
  h2->SetStats(0);
  h2->SetXTitle( type );
  h1->SetXTitle( type );
  //std::cout << "GET X Title Size:  " << RATIO->GetYaxis()->GetTitleSize() << std::endl;
  TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();
  if(h1->GetBinContent(nbins) != 0.0){
    if(h2->GetBinContent(nbins)-h2->GetBinError(nbins) >0){
      h1->GetYaxis()->SetRangeUser(0.2*(h2->GetBinContent(nbins)-h2->GetBinError(nbins)), 2.0*h1->GetBinContent(1));
    }else{
      h1->GetYaxis()->SetRangeUser(0.1, 2*h1->GetBinContent(1));
    }
  }else{
    h1->GetYaxis()->SetRangeUser(0.2, 2*h1->GetBinContent(1));
  }
  h1->Draw("pe");
  h2->DrawCopy("E2same");
  h2clone->DrawCopy("hist same");
  h1->Draw("pesame");
  pad1->Update();
  C->cd();
  
  switch(color){
  case 0:
    h2->SetLineColor(kGreen);
    break;
  case 1:
    h2->SetLineColor(kRed);
    break;
  case 2:
    h2->SetLineColor(kBlue);
    break;
  default:
    std::cout << "---Default---" << std::endl;
    break;
  }
  leg = new TLegend(0.675, 0.83, 0.89, 0.92);//(xmin, ymin, xmax, ymax)
  leg->AddEntry(h1, label + " " + h1Name ,"lep");
  leg->AddEntry(h2, label + " " + h2Name ,"lf");
  leg->SetTextSize(.018);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->Draw();
  //pad1->SetLogy();
  C->Update();

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  t->SetTextSize(0.03);
  t->DrawLatex(0.22,0.95,"CMS Preliminary:");
  t->DrawLatex(0.42,0.95,"#sqrt{s} = 8 TeV,");
  t->DrawLatex(0.62,0.95,"#int L dt = 18.836 fb^{-1}");
  
  TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.25);
  pad2->SetTopMargin(0);
  pad2->SetTopMargin(0.008);
  pad2->SetBottomMargin(0.25);
  pad2->SetGridy();
  //pad2->Draw();
  //pad2->cd();

  RATIO->SetLineColor(1);
  RATIO->SetStats(0);
  RATIO->SetTitle("");
  RATIO->GetXaxis()->SetLabelSize(0.1);
  RATIO->GetYaxis()->SetLabelSize(0.08);
  RATIO->GetYaxis()->SetTitleOffset(0.35);
  RATIO->GetXaxis()->SetTitleOffset(0.88);
  RATIO->GetYaxis()->SetTitleSize(0.11);
  RATIO->GetXaxis()->SetTitleSize(0.11);
  RATIO->SetXTitle( label );
  RATIO->SetYTitle("Ratio");
  RATIO->SetLineColor(1);
  RATIO->SetMarkerSize(.3);
  RATIO->SetMarkerColor(1);
  RATIO->SetMarkerStyle(20);
  switch(color){
  case 0:
    RATIO->SetFillColor(kGreen-10);
    break;
  case 1:
    RATIO->SetFillColor(kRed-10);
    break;
  case 2:
    RATIO->SetFillColor(kBlue-10);
    break;
  default:
    std::cout << "---Default---" << std::endl;
    break;
  }
  //RATIO->SetFillColor(kGreen-10);

  RATIO2->SetStats(0);
  RATIO2->SetTitle("");
  RATIO2->GetXaxis()->SetLabelSize(0.1);
  RATIO2->GetYaxis()->SetLabelSize(0.08);
  RATIO2->GetYaxis()->SetTitleOffset(0.35);
  RATIO2->GetXaxis()->SetTitleOffset(0.88);
  RATIO2->GetYaxis()->SetTitleSize(0.11);
  RATIO2->GetXaxis()->SetTitleSize(0.11);
  RATIO2->SetXTitle( label );
  RATIO2->SetYTitle("Ratio");
  RATIO2->SetLineWidth(2);
  TH1F* RATIO2clone = (TH1F*)RATIO2->Clone("ratio2clone");
  RATIO2clone->SetFillColor(0);
  RATIO2clone->SetLineColor(kGreen);
  switch(color){
  case 0:
    RATIO2->SetLineColor(kGreen-10);
    RATIO2->SetFillColor(kGreen-10);
    RATIO2clone->SetLineColor(kGreen);
    break;
  case 1:
    RATIO2->SetLineColor(kRed-10);
    RATIO2->SetFillColor(kRed-10);
    RATIO2clone->SetLineColor(kRed);
    break;
  case 2:
    RATIO2->SetLineColor(kBlue-10);
    RATIO2->SetFillColor(kBlue-10);
    RATIO2clone->SetLineColor(kBlue);
    break;
  default:
    std::cout << "---Default---" << std::endl;
    break;
  }
  
  //RATIO2->SetFillColor(kGreen-10);
  //RATIO2->SetLineColor(kGreen-10);
  

  //RATIO->SetFillColor(0);
  RATIO->SetLineWidth(2);
  pad2->Draw();
  pad2->cd();
  RATIO->Draw("pe");
  RATIO2->DrawCopy("e2 same");
  RATIO2clone->DrawCopy("hist same");
  RATIO->Draw("pesame");
  
  C->cd();
  C->SaveAs(fname + ".C");
  C->SaveAs(fname + ".root");
  C->SaveAs(fname + ".pdf");
  C->SaveAs(fname + ".png");
  
  
  delete leg;
  delete C;
  delete RATIO;
  
  return 0;
  
};

int RatioPlotsBandMC(TH1F* h1, TH1F* h2, TString h1Name = "h1Name", TString h2Name = "h2Name", TString fname = "default_name", TString type = "defaultType", int nbins = 0, float* bins = NULL, int color = 0){
	
  TCanvas* C = new TCanvas("C", "C	", 400, 500);
  C->cd();  
  TLegend* leg;
  
  TH1F*  RATIO;// = new TH1F("RATIO","Data_to_Prediction",BaseDM::MR_Bins, BaseDM::MR_BinArr);
  TH1F*  RATIO2;
  TString label;
  if(type == "MR"){
    RATIO = new TH1F("RATIO", fname + "_" + type , nbins, bins);
    RATIO2 = new TH1F("RATIO2", fname + "_" + type , nbins, bins);
    label = "M_{R}";
    h1->GetXaxis()->SetRangeUser(200,3500);
    h2->GetXaxis()->SetRangeUser(200,3500);
    RATIO->GetXaxis()->SetRangeUser(200,3500);

  }else if(type == "RSQ" ){
    RATIO = new TH1F("RATIO", fname + "_" + type , nbins, bins);
    RATIO2 = new TH1F("RATIO2", fname + "_" + type , nbins, bins);
    label = "R^{2}";
    h1->GetXaxis()->SetRangeUser(0.5, 1.2);
    h2->GetXaxis()->SetRangeUser(0.5, 1.2);
    RATIO->GetXaxis()->SetRangeUser(0.5, 1.2);
  }else{
    delete RATIO;
    delete RATIO2;
    delete C;
    std::cout << "Unknown Type, please use: MR or RSQ" << std::endl;
    return -1;
  }
  
  //RATIO->Divide(h1, h2, 1, 1, "");
  
  for(int i = 1; i <= h1->GetNbinsX(); i++){
    for(int j = 1; j <=  h1->GetNbinsY(); j++){
      double r = .0, r2 = .0;
      double err = .0, err2 = .0;
      
      if(h2->GetBinContent(i,j) !=0 ){
	r = h1->GetBinContent(i,j)/h2->GetBinContent(i,j);
	err = h1->GetBinError(i,j)/h2->GetBinContent(i,j);
	r2 = h2->GetBinContent(i,j)/h2->GetBinContent(i,j);
	err2 = h2->GetBinError(i,j)/h2->GetBinContent(i,j);
      }
      RATIO->SetBinContent(i,j,r);
      RATIO->SetBinError(i,j,err);
      RATIO2->SetBinContent(i,j,r2);
      RATIO2->SetBinError(i,j,err2);
      std::cout << "Error: " << h2->GetBinError(i,j) << " " << err << std::endl;
      std::cout << "BinContent: " << h1->GetBinContent(i,j) << std::endl;
    }
  }
  RATIO->GetYaxis()->SetRangeUser(.95, 1.05);
  RATIO2->GetYaxis()->SetRangeUser(.95, 1.05);
   
  RATIO->Sumw2();
  RATIO2->Sumw2();
  h1->SetLineColor(1);//data
  h1->SetMarkerSize(.3);
  h1->SetMarkerColor(1);
  //h1->SetMarkerStyle(20);
  //h1->SetFillColor(1);
  h1->SetLineWidth(1.0);
  //h1->SetLineStyle(2);
  h2->SetLineWidth(4.0);
  TH1F* h2clone = (TH1F*)h2->Clone("h2clone");
  h2clone->SetFillColor(0);
  TH1F* h1clone = (TH1F*)h1->Clone("h1clone");
  h1clone->SetFillColor(0);
  switch(color){
  case 0:
    h2->SetLineColor(kGreen-10);
    h2->SetFillColor(kGreen-10);
    h2clone->SetLineColor(kGreen);
    break;
  case 1:
    h2->SetLineColor(kRed-10);
    h2->SetFillColor(kRed-10);
    h2clone->SetLineColor(kRed);
    break;
  case 2:
    h2->SetLineColor(kBlue-10);
    h2->SetFillColor(kBlue-10);
    h2clone->SetLineColor(kBlue);
    break;
  default:
    std::cout << "---Default---" << std::endl;
    break;
  }
  h1->SetStats(0);
  h1->SetTitle("");
  h2->SetTitle("");
  h2->SetStats(0);
  h2->SetXTitle( type );
  h1->SetXTitle( type );
  std::cout << "GET X Title Size:  " << RATIO->GetYaxis()->GetTitleSize() << std::endl;
  TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();
  if(h1->GetBinContent(nbins) != 0.0){
    if(h2->GetBinContent(nbins)-h2->GetBinError(nbins) >0){
      h1->GetYaxis()->SetRangeUser(0.2*(h2->GetBinContent(nbins)-h2->GetBinError(nbins)), 1.3*h1->GetBinContent(1));
    }else{
      h1->GetYaxis()->SetRangeUser(0.1, 1.3*h1->GetBinContent(1));
    }
  }else{
    h1->GetYaxis()->SetRangeUser(0.2, 2*h1->GetBinContent(1));
  }
  h1->Draw("hist p");
  h2->DrawCopy("E2same");
  h2clone->DrawCopy("hist same");
  h1clone->Draw("hist same");
  pad1->Update();
  C->cd();
  
  switch(color){
  case 0:
    h2->SetLineColor(kGreen);
    break;
  case 1:
    h2->SetLineColor(kRed);
    break;
  case 2:
    h2->SetLineColor(kBlue);
    break;
  default:
    std::cout << "---Default---" << std::endl;
    break;
  }
  leg = new TLegend(0.6, 0.83, 0.85, 0.92);//(xmin, ymin, xmax, ymax)
  leg->AddEntry(h1, label + " " + h1Name ,"l");
  leg->AddEntry(h2, label + " " + h2Name ,"lf");
  leg->SetTextSize(.018);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->Draw();
  //pad1->SetLogy();
  C->Update();

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  t->SetTextSize(0.03);
  t->DrawLatex(0.22,0.95,"CMS Preliminary:");
  t->DrawLatex(0.42,0.95,"#sqrt{s} = 8 TeV,");
  t->DrawLatex(0.62,0.95,"#int L dt = 18.836 fb^{-1}");
  
  TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.25);
  pad2->SetTopMargin(0);
  pad2->SetTopMargin(0.008);
  pad2->SetBottomMargin(0.25);
  pad2->SetGridy();
  //pad2->Draw();
  //pad2->cd();

  RATIO->SetLineColor(1);
  RATIO->SetStats(0);
  RATIO->SetTitle("");
  RATIO->GetXaxis()->SetLabelSize(0.1);
  RATIO->GetYaxis()->SetLabelSize(0.08);
  RATIO->GetYaxis()->SetTitleOffset(0.35);
  RATIO->GetXaxis()->SetTitleOffset(0.88);
  RATIO->GetYaxis()->SetTitleSize(0.11);
  RATIO->GetXaxis()->SetTitleSize(0.11);
  RATIO->SetXTitle( label );
  RATIO->SetYTitle("Ratio");
  RATIO->SetLineColor(1);
  RATIO->SetMarkerSize(.3);
  RATIO->SetMarkerColor(1);
  RATIO->SetMarkerStyle(20);
  switch(color){
  case 0:
    RATIO->SetFillColor(kGreen-10);
    break;
  case 1:
    RATIO->SetFillColor(kRed-10);
    break;
  case 2:
    RATIO->SetFillColor(kBlue-10);
    break;
  default:
    std::cout << "---Default---" << std::endl;
    break;
  }
  //RATIO->SetFillColor(kGreen-10);

  RATIO2->SetStats(0);
  RATIO2->SetTitle("");
  RATIO2->GetXaxis()->SetLabelSize(0.1);
  RATIO2->GetYaxis()->SetLabelSize(0.08);
  RATIO2->GetYaxis()->SetTitleOffset(0.35);
  RATIO2->GetXaxis()->SetTitleOffset(0.88);
  RATIO2->GetYaxis()->SetTitleSize(0.11);
  RATIO2->GetXaxis()->SetTitleSize(0.11);
  RATIO2->SetXTitle( label );
  RATIO2->SetYTitle("Ratio");
  RATIO2->SetLineWidth(2);
  TH1F* RATIO2clone = (TH1F*)RATIO2->Clone("ratio2clone");
  RATIO2clone->SetFillColor(0);
  RATIO2clone->SetLineColor(kGreen);
  switch(color){
  case 0:
    RATIO2->SetLineColor(kGreen-10);
    RATIO2->SetFillColor(kGreen-10);
    RATIO2clone->SetLineColor(kGreen);
    break;
  case 1:
    RATIO2->SetLineColor(kRed-10);
    RATIO2->SetFillColor(kRed-10);
    RATIO2clone->SetLineColor(kRed);
    break;
  case 2:
    RATIO2->SetLineColor(kBlue-10);
    RATIO2->SetFillColor(kBlue-10);
    RATIO2clone->SetLineColor(kBlue);
    break;
  default:
    std::cout << "---Default---" << std::endl;
    break;
  }
  
  //RATIO2->SetFillColor(kGreen-10);
  //RATIO2->SetLineColor(kGreen-10);
  

  //RATIO->SetFillColor(0);
  RATIO->SetLineWidth(1);
  pad2->Draw();
  pad2->cd();
  RATIO->Draw("pe");
  RATIO2->DrawCopy("e2 same");
  RATIO2clone->DrawCopy("hist same");
  RATIO->Draw("pesame");
  
  C->cd();
  C->SaveAs(fname + ".C");
  C->SaveAs(fname + ".root");
  C->SaveAs(fname + ".pdf");
  C->SaveAs(fname + ".png");
  
  
  delete leg;
  delete C;
  delete RATIO;
  
  return 0;
  
};

int BandMC_TGraph(TGraphAsymmErrors* h1, TH1F* h2, TString h1Name = "h1Name", TString h2Name = "h2Name", TString fname = "default_name", TString type = "defaultType", int nbins = 0, float* bins = NULL, int color = 0){
	
  TCanvas* C = new TCanvas("C", "C", 400, 500);
  C->cd();  
  TLegend* leg;
  
  TString label;
  
  h2->SetLineWidth(1.1);
  TH1F* h2clone = (TH1F*)h2->Clone("h2clone");
  h2clone->SetFillColor(0);
  switch(color){
  case 0:
    h2->SetLineColor(kGreen-10);
    h2->SetFillColor(kGreen-10);
    h2clone->SetLineColor(kGreen);
    break;
  case 1:
    h2->SetLineColor(kRed-10);
    h2->SetFillColor(kRed-10);
    h2clone->SetLineColor(kRed);
    break;
  case 2:
    h2->SetLineColor(kBlue-10);
    h2->SetFillColor(kBlue-10);
    h2clone->SetLineColor(kBlue);
    break;
  default:
    std::cout << "---Default---" << std::endl;
    break;
  }
  
  h1->GetXaxis()->SetRangeUser(0.5, 1.2);
  h1->GetYaxis()->SetRangeUser(0.01*h2->GetBinContent(1), 1.5*h2->GetBinContent(h2->GetNbinsX()));
  h2->SetXTitle( type );
  h1->GetXaxis()->SetTitle( "R^{2}" );
  h1->GetXaxis()->CenterTitle(1);
  
  TPad *pad1 = new TPad("pad1","pad1",0,0.0,1,1);
  //pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();
  h1->SetFillColor(2);
  h1->SetFillStyle(3001);
  h1->SetTitle("");
    
  h1->Draw("a2same");
  //h2->DrawCopy("E2same][");
  h2clone->DrawCopy("hist same][");
  h1->Draw("2");
  pad1->Update();
  C->cd();
  
  switch(color){
  case 0:
    h2->SetLineColor(kGreen);
    break;
  case 1:
    h2->SetLineColor(kRed);
    break;
  case 2:
    h2->SetLineColor(kBlue);
    break;
  default:
    std::cout << "---Default---" << std::endl;
    break;
  }
  leg = new TLegend(0.68, 0.78, 0.89, 0.89);//(xmin, ymin, xmax, ymax)
  leg->AddEntry(h1, label + " " + h1Name ,"lf");
  leg->AddEntry(h2, label + " " + h2Name ,"lf");
  leg->SetTextSize(.018);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->Draw();
  //pad1->SetLogy();
  C->Update();

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  t->SetTextSize(0.03);
  t->DrawLatex(0.22,0.95,"CMS Preliminary:");
  t->DrawLatex(0.42,0.95,"#sqrt{s} = 8 TeV,");
  t->DrawLatex(0.62,0.95,"#int L dt = 18.836 fb^{-1}");
  
  
  C->cd();
  C->SaveAs(fname + ".C");
  C->SaveAs(fname + ".root");
  C->SaveAs(fname + ".pdf");
  C->SaveAs(fname + ".png");
  
  
  delete leg;
  delete C;
  
  return 0;
  
};

int BandMC_TGraph(TGraphAsymmErrors* h1, TH1F* h2, TGraphAsymmErrors* h3, TH1F* h4, TString h1Name = "h1Name", TString h2Name = "h2Name", TString fname = "default_name", TString type = "defaultType", int nbins = 0, float* bins = NULL, int color = 0){
	
  TCanvas* C = new TCanvas("C", "C", 400, 500);
  C->cd();  
  TLegend* leg;
  
  TString label;
  
  h2->SetLineWidth(1.2);
  TH1F* h2clone = (TH1F*)h2->Clone("h2clone");
  h2clone->SetFillColor(0);
  h4->SetLineWidth(1.2);
  TH1F* h4clone = (TH1F*)h4->Clone("h4clone");
  h4clone->SetFillColor(0);
  switch(color){
  case 0:
    h2->SetLineColor(kGreen-10);
    h2->SetFillColor(kGreen-10);
    h2clone->SetLineColor(kGreen);
    h4->SetLineColor(kGreen-10);
    h4->SetFillColor(kGreen-10);
    h4clone->SetLineColor(kGreen);
    h1->SetFillColor(kGreen-10);
    h1->SetLineColor(kGreen-10);
    h3->SetFillColor(kGreen-10);
    h3->SetLineColor(kGreen-10);
    break;
  case 1:
    h2->SetLineColor(kRed-10);
    h2->SetFillColor(kRed-10);
    h2clone->SetLineColor(kRed);
    h4->SetLineColor(kRed-10);
    h4->SetFillColor(kRed-10);
    h4clone->SetLineColor(kRed);
    h1->SetFillColor(kRed-10);
    h1->SetLineColor(kRed-10);
    h3->SetFillColor(kRed-10);
    h3->SetLineColor(kRed-10);
    break;
  case 2:
    h2->SetLineColor(kBlue-10);
    h2->SetFillColor(kBlue-10);
    h2clone->SetLineColor(kBlue);
    h4->SetLineColor(kBlue-10);
    h4->SetFillColor(kBlue-10);
    h4clone->SetLineColor(kBlue);
    h1->SetFillColor(kBlue-10);
    h1->SetLineColor(kBlue-10);
    h3->SetFillColor(kBlue-10);
    h3->SetLineColor(kBlue-10);
    break;
  default:
    std::cout << "---Default---" << std::endl;
    break;
  }
  
  h1->GetXaxis()->SetRangeUser(0.5, 1.2);
  if(h2->GetBinContent(1) < h2->GetBinContent(h2->GetNbinsX())){
    h1->GetYaxis()->SetRangeUser(0.01*h2->GetBinContent(1), 1.5*h2->GetBinContent(h2->GetNbinsX()));
  }else{
    h1->GetYaxis()->SetRangeUser(0.01*h2->GetBinContent(h2->GetNbinsX()), 1.5*h2->GetBinContent(3));
  }
  h2->SetXTitle( type );
  h1->GetXaxis()->SetTitle( "R^{2}" );
  h1->GetXaxis()->CenterTitle(1);
  
  TPad *pad1 = new TPad("pad1","pad1",0.0,0.25,1,1);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();
  h1->SetFillStyle(3001);
  h1->SetTitle("");
    
  h1->Draw("a2same");
  //h2->DrawCopy("E2same][");
  h2clone->DrawCopy("hist same][");
  h1->Draw("2");
  pad1->Update();
  C->cd();
  
  switch(color){
  case 0:
    h2->SetLineColor(kGreen);
    h4->SetLineColor(kGreen);
    break;
  case 1:
    h2->SetLineColor(kRed);
    h4->SetLineColor(kRed);
    break;
  case 2:
    h2->SetLineColor(kBlue);
    h4->SetLineColor(kBlue);
    break;
  default:
    std::cout << "---Default---" << std::endl;
    break;
  }
  leg = new TLegend(0.68, 0.78, 0.89, 0.89);//(xmin, ymin, xmax, ymax)
  leg->AddEntry(h1, label + " " + h1Name ,"lf");
  leg->AddEntry(h2clone, label + " " + h2Name ,"lf");
  leg->SetTextSize(.018);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->Draw();
  //pad1->SetLogy();
  C->Update();

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  t->SetTextSize(0.03);
  t->DrawLatex(0.22,0.95,"CMS Preliminary:");
  t->DrawLatex(0.42,0.95,"#sqrt{s} = 8 TeV,");
  t->DrawLatex(0.62,0.95,"#int L dt = 18.836 fb^{-1}");
  
  TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.25);
  pad2->SetTopMargin(0);
  pad2->SetTopMargin(0.008);
  pad2->SetBottomMargin(0.25);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
  
  h3->GetXaxis()->SetTitle( "R^{2}" );
  h3->GetXaxis()->CenterTitle(1);
  h3->GetXaxis()->SetLabelSize(0.1);
  h3->GetYaxis()->SetLabelSize(0.08);
  h3->GetYaxis()->SetTitleOffset(0.35);
  h3->GetXaxis()->SetTitleOffset(0.88);
  h3->GetYaxis()->SetTitleSize(0.11);
  h3->GetXaxis()->SetTitleSize(0.11);
  
  //h3->SetFillColor(2);
  h3->SetFillStyle(3001);
  h3->SetTitle("");
  h3->GetXaxis()->SetRangeUser(0.5, 1.2);
  h3->GetYaxis()->SetRangeUser(0.8, 1.2);
  h3->Draw("a2same");
  //h2->DrawCopy("E2same][");
  h4clone->DrawCopy("hist same][");
  h3->Draw("2");
  pad2->Update();
  
  C->cd();
  C->SaveAs(fname + ".C");
  C->SaveAs(fname + ".root");
  C->SaveAs(fname + ".pdf");
  C->SaveAs(fname + ".png");
  
  
  delete leg;
  delete C;
  
  return 0;
  
};

int RatioPlotsBand(TH1F* h1, TH1F* h2, TString h1Name = "h1Name", TString h2Name = "h2Name", TString fname = "default_name", TString type = "defaultType"){
	
  TCanvas* C = new TCanvas("C", "C	", 400, 500);
  C->cd();  
  TLegend* leg;
  
  TH1F*  RATIO;// = new TH1F("RATIO","Data_to_Prediction",BaseDM::MR_Bins, BaseDM::MR_BinArr);
  TH1F*  RATIO2;
  TString label;
  if(type == "MR"){
    RATIO = new TH1F("RATIO", fname + "_" + type , BaseDM::MR_Bins, BaseDM::MR_BinArr);
    RATIO2 = new TH1F("RATIO2", fname + "_" + type , BaseDM::MR_Bins, BaseDM::MR_BinArr);
    label = "M_{R}";
    h1->GetXaxis()->SetRangeUser(200,3500);
    h2->GetXaxis()->SetRangeUser(200,3500);
    RATIO->GetXaxis()->SetRangeUser(200,3500);

  }else if(type == "RSQ" ){
    RATIO = new TH1F("RATIO", fname + "_" + type , BaseDM::RSQ_Bins, BaseDM::RSQ_BinArr);
    RATIO2 = new TH1F("RATIO2", fname + "_" + type , BaseDM::RSQ_Bins, BaseDM::RSQ_BinArr);
    label = "R^{2}";
    h1->GetXaxis()->SetRangeUser(0.5, 2.5);
    h2->GetXaxis()->SetRangeUser(0.5, 2.5);
    RATIO->GetXaxis()->SetRangeUser(0.5, 2.5);
  }else{
    delete RATIO;
    delete RATIO2;
    delete C;
    std::cout << "Unknown Type, please use: MR or RSQ" << std::endl;
    return -1;
  }
  
  //RATIO->Divide(h1, h2, 1, 1, "");
  
  for(int i = 1; i <= h1->GetNbinsX(); i++){
    for(int j = 1; j <=  h1->GetNbinsY(); j++){
      double r = .0, r2 = .0;
      double err = .0, err2 = .0;
      
      if(h2->GetBinContent(i,j) !=0 ){
	r = h1->GetBinContent(i,j)/h2->GetBinContent(i,j);
	err = h1->GetBinError(i,j)/h2->GetBinContent(i,j);
	r2 = h2->GetBinContent(i,j)/h2->GetBinContent(i,j);
	err2 = h2->GetBinError(i,j)/h2->GetBinContent(i,j);
      }
      RATIO->SetBinContent(i,j,r);
      RATIO->SetBinError(i,j,err);
      RATIO2->SetBinContent(i,j,r2);
      RATIO2->SetBinError(i,j,err2);
      std::cout << "Error: " << h2->GetBinError(i,j) << " " << err << std::endl;
      std::cout << "BinContent: " << h1->GetBinContent(i,j) << std::endl;
    }
  }
  RATIO->GetYaxis()->SetRangeUser(.0, 2.05);
  RATIO2->GetYaxis()->SetRangeUser(.0, 2.05);
  //RATIO->GetYaxis()->SetRangeUser(.7, 1.3);
  //RATIO2->GetYaxis()->SetRangeUser(.7, 1.3);
  
 
  RATIO->Sumw2();
  RATIO2->Sumw2();
  h1->SetLineColor(1);//data
  h1->SetMarkerSize(.3);
  h1->SetMarkerColor(1);
  h1->SetMarkerStyle(20);
  h1->SetFillColor(1);
  h1->SetLineWidth(1.7);
  
  h2->SetLineColor(kGreen-10);
  h2->SetFillColor(kGreen-10);
  h2->SetLineWidth(2);
  TH1F* h2clone = (TH1F*)h2->Clone("h2clone");
  h2clone->SetFillColor(0);
  h2clone->SetLineColor(kGreen);
  
  h1->SetStats(0);
  h1->SetTitle("");
  h2->SetTitle("");
  h2->SetStats(0);
  h2->SetXTitle( type );
  h1->SetXTitle( type );
  
  std::cout << "GET X Title Size:  " << RATIO->GetYaxis()->GetTitleSize() << std::endl;
   
  TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();
  h1->GetYaxis()->SetRangeUser(80.5, 100000);
  h1->Draw("pe");
  h2->DrawCopy("E2same");
  h2clone->DrawCopy("hist same");
  h1->Draw("pesame");
  pad1->Update();
  C->cd();
  
  h2->SetLineColor(kGreen);
  leg = new TLegend(0.65, 0.7, 0.89, 0.9);//(xmin, ymin, xmax, ymax)
  leg->AddEntry(h1, label + " " + h1Name ,"lep");
  leg->AddEntry(h2, label + " " + h2Name ,"lf");
  leg->SetTextSize(.02);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->Draw();
  pad1->SetLogy();
  C->Update();

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  t->SetTextSize(0.03);
  t->DrawLatex(0.22,0.95,"CMS Preliminary:");
  t->DrawLatex(0.42,0.95,"#sqrt{s} = 8 TeV,");
  t->DrawLatex(0.62,0.95,"#int L dt = 18.836 fb^{-1}");
  
  TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.25);
  pad2->SetTopMargin(0);
  pad2->SetTopMargin(0.008);
  pad2->SetBottomMargin(0.25);
  pad2->SetGridy();
  //pad2->Draw();
  //pad2->cd();

  RATIO->SetLineColor(1);
  RATIO->SetStats(0);
  RATIO->SetTitle("");
  RATIO->GetXaxis()->SetLabelSize(0.1);
  RATIO->GetYaxis()->SetLabelSize(0.08);
  RATIO->GetYaxis()->SetTitleOffset(0.35);
  RATIO->GetXaxis()->SetTitleOffset(0.88);
  RATIO->GetYaxis()->SetTitleSize(0.11);
  RATIO->GetXaxis()->SetTitleSize(0.11);
  RATIO->SetXTitle( label );
  RATIO->SetYTitle("Ratio");
  RATIO->SetLineColor(1);
  RATIO->SetMarkerSize(.3);
  RATIO->SetMarkerColor(1);
  RATIO->SetMarkerStyle(20);
  RATIO->SetFillColor(kGreen-10);

  RATIO2->SetStats(0);
  RATIO2->SetTitle("");
  RATIO2->GetXaxis()->SetLabelSize(0.1);
  RATIO2->GetYaxis()->SetLabelSize(0.08);
  RATIO2->GetYaxis()->SetTitleOffset(0.35);
  RATIO2->GetXaxis()->SetTitleOffset(0.88);
  RATIO2->GetYaxis()->SetTitleSize(0.11);
  RATIO2->GetXaxis()->SetTitleSize(0.11);
  RATIO2->SetXTitle( label );
  RATIO2->SetYTitle("Ratio");
  RATIO2->SetFillColor(kGreen-10);
  RATIO2->SetLineColor(kGreen-10);
  RATIO2->SetLineWidth(2);
  
  TH1F* RATIO2clone = (TH1F*)RATIO2->Clone("ratio2clone");
  RATIO2clone->SetFillColor(0);
  RATIO2clone->SetLineColor(kGreen);

  //RATIO->SetFillColor(0);
  RATIO->SetLineWidth(2);
  pad2->Draw();
  pad2->cd();
  RATIO->Draw("pe");
  RATIO2->DrawCopy("e2 same");
  RATIO2clone->DrawCopy("hist same");
  RATIO->Draw("pesame");
  
  C->cd();
  C->SaveAs(fname + ".pdf");
  C->SaveAs(fname + ".png");
  
  delete leg;
  delete C;
  delete RATIO;
  
  return 0;
  
};



int RatioPlotsV2(THStack* s, TH1F* h1, TH1F* h2, TString h1Name = "h1Name", TString h2Name = "h2Name", TString fname = "default_name", TString type = "defaultType", int nbins = 0, float* bins = NULL, TLegend* le = NULL){

  TCanvas* C = new TCanvas("C", "C", 400, 500);
  C->cd();

  TH1F*  RATIO;
  TString label;
  std::cout << "debug 0" << std::endl;
  if(type == "MR"){
    RATIO = new TH1F("RATIO", fname + "_" + type , BaseDM::MR_Bins, BaseDM::MR_BinArr);
    label = "M_{R}";
    h1->GetXaxis()->SetRangeUser(200.,3500.);
    h2->GetXaxis()->SetRangeUser(200.,3500.);
    RATIO->GetXaxis()->SetRangeUser(200.,3500.);
    RATIO->GetYaxis()->SetRangeUser(.0, 2.0);
    s->SetMaximum(100000.);
  }else if(type == "RSQ" ){
    std::cout << "debug 1" << std::endl;
    RATIO = new TH1F("RATIO", fname + "_" + type , nbins, bins);
    //RATIO2 = new TH1F("RATIO2", fname + "_" + type , nbins, bins);
    label = "R^{2}";
    std::cout << "debug 2" << std::endl;
    h1->GetXaxis()->SetRangeUser(0.5, 1.2);
    h2->GetXaxis()->SetRangeUser(0.5, 1.2);
    RATIO->GetXaxis()->SetRangeUser(0.5, 1.2);
    std::cout << "debug 3" << std::endl;
  }else if(type == "MET"){
    RATIO = new TH1F("RATIO", fname + "_" + type , 50, 0, 1000);
    label = "#slash{E}_{T}  GeV";
    s->SetMaximum(10000.);
  }else if(type == "NJETS"){
    RATIO = new TH1F("RATIO", fname + "_" + type , 9, 1, 10);
    label = "Jet Multiplicity";
    s->SetMaximum(100000.);
  }else{
    delete RATIO;
    delete C;
    std::cout << "Unknown Type, please use: MR or RSQ" << std::endl;
    return -1;
  }
  
  std::cout << "=====================Dividing Histograms=====================" << std::endl;
  RATIO->Divide(h1, h2, 1, 1, "");
  RATIO->GetYaxis()->SetRangeUser(.0, 2.05);
  h1->SetMarkerSize(.7);
  h1->SetStats(0);
  s->SetMinimum(1.);
  TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
  pad1->SetBottomMargin(0.0);
  pad1->Draw();
  pad1->cd();
  if(h1->GetBinContent(nbins) != 0.0){
    s->SetMinimum(0.2*h1->GetBinContent(nbins));
    s->SetMaximum(4*h1->GetBinContent(1));
  }else if(h1->GetBinContent(1) != 0.0){
    s->SetMinimum(0.2);
    s->SetMaximum(4*h1->GetBinContent(1));
  }else{
    s->SetMinimum(0.2);
    s->SetMaximum(10.0);
  }
  s->SetTitle("");
  s->Draw();
  if(type == "MET"){
    s->GetYaxis()->SetTitle("Events/20 GeV");
  }else{
    s->GetYaxis()->SetTitle("Events");
  }
  s->GetYaxis()->SetTitleOffset(1.25);
  gPad->Modified();
  h1->SetStats(0);
  h1->Draw("same");
  C->cd();
  
  le->SetFillColor(0);
  le->SetBorderSize(0);
  le->SetTextSize(0.02);
  le->Draw();
  pad1->SetLogy();
  C->Update();
  
  std::cout << "debug 5" << std::endl;
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  t->SetTextSize(0.03);
  t->DrawLatex(0.22,0.95,"CMS Preliminary:");
  t->DrawLatex(0.42,0.95,"#sqrt{s} = 8 TeV,");
  t->DrawLatex(0.62,0.95,"#int L dt = 18.84 fb^{-1}");
  
  TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.25);
  pad2->SetTopMargin(0.008);
  pad2->SetBottomMargin(0.25);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
  RATIO->SetLineColor(4);
  RATIO->SetStats(0);
  RATIO->SetTitle("");
  RATIO->GetXaxis()->SetLabelSize(0.1); 
  RATIO->GetYaxis()->SetLabelSize(0.1); 
  RATIO->GetYaxis()->SetTitleOffset(0.35);
  RATIO->GetXaxis()->SetTitleOffset(0.88);
  RATIO->GetYaxis()->SetTitleSize(0.11);
  RATIO->GetXaxis()->SetTitleSize(0.11);
  RATIO->SetXTitle( label );
  RATIO->SetYTitle("Data/MC");
  RATIO->SetLineColor(4);
  RATIO->SetMarkerSize(.7);
  RATIO->SetMarkerColor(4);
  RATIO->SetMarkerStyle(20);
  RATIO->SetFillColor(4);
  RATIO->Draw();
  C->cd();
  C->SaveAs(fname + ".pdf");
  C->SaveAs(fname + ".png");
  C->SaveAs(fname + ".C");
 
  delete C;
  delete RATIO;

  return 0;
  
};

int PlotCosmetics(TH1F* h, TString fname = "defaultFN" , TString xlabel = "X", TString ylabel = "Y", TString gopt = "", TString mASS = "", TString current = ""){

  TCanvas* C = new TCanvas("C", "C", 400, 500);
  C->cd();
  
  h->SetLineColor(kGreen+2);
  h->SetLineStyle(2);
  h->SetLineWidth(2);
  h->SetFillColor(kGreen+2);
  h->SetFillStyle(3002);
  
  TLegend* leg = new TLegend(0.54, 0.84, 0.89, 0.88);//(xmin, ymin, xmax, ymax)
  leg->SetTextAlign(22);
  leg->AddEntry(h, " DM, m = "+mASS+" GeV, "+current ,"l");
  leg->SetTextSize(.028);
  leg->SetFillColor(0);
  leg->SetLineColor(0);


  h->SetStats(0);
  h->SetTitle("");
  h->GetXaxis()->SetLabelSize(0.03); 
  h->GetYaxis()->SetLabelSize(0.03);
  h->GetYaxis()->CenterTitle(1);
  h->GetXaxis()->CenterTitle(1);
  h->GetYaxis()->SetTitleOffset(1.5);
  h->GetXaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetTitleSize(0.035);
  h->GetXaxis()->SetTitleSize(0.035);
  h->SetXTitle( xlabel );
  h->SetYTitle( ylabel );
  h->Draw();
  
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  t->SetTextSize(0.03);
  t->DrawLatex(0.22,0.95,"CMS Preliminary:");
  t->DrawLatex(0.42,0.95,"#sqrt{s} = 8 TeV,");
  t->DrawLatex(0.62,0.95,"#int L dt = 18.84 fb^{-1}");
  
  
  leg->Draw();

  C->SaveAs(fname + ".pdf");
  C->SaveAs(fname + ".png");
  C->SaveAs(fname + ".C");
  
  return 0;

};

int PlotCosmetics2D(TH2F* h, TString fname = "defaultFN" , TString xlabel = "X", TString ylabel = "Y", TString gopt = "", TString mASS = "", TString current = ""){

  TCanvas* C = new TCanvas("C", "C", 400, 500);
  C->cd();

  h->SetStats(0);
  h->SetTitle("");
  h->GetXaxis()->SetLabelSize(0.03); 
  h->GetYaxis()->SetLabelSize(0.03);
  h->GetYaxis()->CenterTitle(1);
  h->GetXaxis()->CenterTitle(1);
  h->GetYaxis()->SetTitleOffset(1.5);
  h->GetXaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetTitleSize(0.035);
  h->GetXaxis()->SetTitleSize(0.035);
  h->SetXTitle( xlabel );
  h->SetYTitle( ylabel );
  h->Draw(gopt);
  
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  t->SetTextSize(0.03);
  t->DrawLatex(0.22,0.95,"CMS Preliminary:");
  t->DrawLatex(0.42,0.95,"#sqrt{s} = 8 TeV,");
  t->DrawLatex(0.62,0.95,"#int L dt = 18.84 fb^{-1}");
  
  
  C->SaveAs(fname + ".pdf");
  C->SaveAs(fname + ".png");
  C->SaveAs(fname + ".C");
  
  return 0;

};
