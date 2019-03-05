//
//  _Peaks_research
//  
//
//  Created by Luca Tosti on 05/03/2019.
//

#ifndef functions_sipm_Peaks_research_h
#define functions_sipm_Peaks_research_h

#include "../include/sipm/general_lib.h"
#include "../include/sipm/Peaks.h"

//Peaks research
void initpeaks( Peaks *p );
void FindPeaks( TGraph* g, Peaks *peaks, Double_t threshold, bool average, bool fixleftside );
void FindNextPeaks( TGraph *g, Peaks* peaksdled, Peaks *peaks, Double_t twindow, bool average );
void RemovePeaks( TGraph *g, TGraph **gnopeak, Peaks *peaks, Double_t tleft=40e-9, Double_t tright=300e-9);


#endif /* functions_sipm_Peaks_research.h */
