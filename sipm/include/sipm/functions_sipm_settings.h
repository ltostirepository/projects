//
//  sipm_settings_h
//  
//
//  Created by Luca Tosti on 05/03/2019.
//

#ifndef functions_sipm_settings_h
#define functions_sipm_settings_h

#include "../include/sipm/general_lib.h"
#include "../include/sipm/Peaks.h"

// settings
void ScaleGraph( TGraph* g, Double_t scale);
Double_t* GenerateLogBinning(int nbins, float min, float max);
Double_t* GenerateBinning(int nbins, float min, float max);

#endif /* functions_sipm_settings.h */
