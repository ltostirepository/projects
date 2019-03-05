//
//  _smoothing_and_DLED
//  
//
//  Created by Luca Tosti on 05/03/2019.
//

#ifndef functions_sipm_smoothing_and_DLED_h
#define functions_sipm_smoothing_and_DLED_h

#include "../include/sipm/general_lib.h"
#include "../include/sipm/Peaks.h"

// smoothing and DLED
void TGraphRunningMean( TGraph *g, TGraph **gm, Int_t windowsize);
void gDLEDpoint( const TGraph* gorig, TGraph **gdelayed, Int_t pdelay);
void SetHTime(TH1* h);
TH2F* NormalizeX(TH2F* hh);

#endif /* functions_sipm_smoothing_and_DLED.h */
