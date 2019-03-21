//
//  _Peaks_charge
//  
//
//  Created by Luca Tosti on 05/03/2019.
//

#ifndef functions_sipm_Peaks_charge_h
#define functions_sipm_Peaks_charge_h


//Peaks charge
void integra_picco( TGraph *g, Peaks *peaksdled,Peaks *peaks, Double_t ped_signal , Double_t time_before, Double_t time_after, Double_t sampling_time_ns);
Bool_t GoodCharge( Peaks* p, Short_t index);

#endif /* functions_sipm_Peaks_charge.h */
