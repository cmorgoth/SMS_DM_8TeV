#ifndef HELPER_FUNCTIONS_HH
#define HELPER_FUCTIONS_HH 1
//C++ INCLUDES
#include <iostream>
#include <stdlib.h>
#include <vector>

//ROOT INCLUDES
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TString.h"

//LOCAL INCLUDES
#include "StructDefinition.hh"

double GetGenEvents(double mass1, double mass2, TH2F* h);
bool CreateStructHistos(pdf_and_ISR& tmp, std::string s);

#endif
