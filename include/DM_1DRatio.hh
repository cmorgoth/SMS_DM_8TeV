#ifndef DM_1DRATIOS
#define DM_1DRATIOS 1

#include <iostream>
#include <sstream>
#include "TCanvas.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TPad.h"
#include "TString.h"
#include "DM_Base.hh"
#include "THStack.h"
#include "TGraphAsymmErrors.h"

int RatioPlots( TH1F*, TH1F*, TString, TString, TString, TString );
int RatioPlotsBand( TH1F*, TH1F*, TString, TString, TString, TString );
int RatioPlotsBandV2( TH1F*, TH1F*, TString, TString, TString, TString, int nbins, float* bins, int color);
int RatioPlotsV2(THStack*, TH1F*, TH1F*, TString, TString, TString, TString, int, float*, TLegend*);
int RatioPlotsBandMC( TH1F*, TH1F*, TString, TString, TString, TString, int nbins, float* bins, int color);
int BandMC_TGraph(TGraphAsymmErrors*, TH1F*, TString, TString, TString, TString, int, float*, int);
int BandMC_TGraph(TGraphAsymmErrors*, TH1F*, TGraphAsymmErrors*, TH1F*, TString, TString, TString, TString, int, float*, int);

int PlotCosmetics(TH1F* h, TString, TString, TString, TString, TString, TString);
int PlotCosmetics2D(TH2F* h, TString, TString, TString, TString, TString, TString);
#endif
