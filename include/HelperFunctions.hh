#ifndef HELPER_FUNCTIONS_HH
#define HELPER_FUCTIONS_HH 1
//C++ INCLUDES
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <stdlib.h>

//ROOT INCLUDES
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"

//LOCAL INCLUDES
#include "StructDefinition.hh"

double GetGenEvents(double mass1, double mass2, TH2F* h);
bool CreateStructHistos(pdf_and_ISR& tmp, std::string s);
bool CreateStructHistos(SignalStruct& tmp, std::string s);
bool GetMasses(std::string aux, float& mass1, float& mass2);

bool PlotLimit(TH2F* h, TGraph* obs, TGraph* exp, TGraph* mSigma, TGraph* pSigma, TString name);

#endif
