#include "../include/sipm/general_lib.h"
#include "../include/sipm/Peaks.h"
#include "../include/sipm/functions_sipm_Peaks_charge.h"

using namespace std;

void integra_picco( TGraph *g, Peaks *peaksdled,Peaks *peaks, Double_t ped_signal , Double_t time_before, Double_t time_after, Double_t sampling_time_ns){
    for(int i_carica = 0 ; i_carica < peaksdled->N;i_carica++){
        Double_t carica_picco=0;
        int bin_before = (int)(time_before/sampling_time_ns);
        int bin_after =  (int)(time_after/sampling_time_ns);
        Double_t tempo = peaksdled->t.at(i_carica);
        //int bin_picco = (int)(tempo/sampling_time_ns);
        int bin_picco = peaksdled->xt.at(i_carica);
        //int running_bin;
        for(int i_carica_tempo = -bin_before ; i_carica_tempo < bin_after  ; i_carica_tempo++){
            if(bin_picco+i_carica_tempo>=0 && bin_picco+i_carica_tempo<g->GetN()){
                carica_picco += g->GetY()[bin_picco+i_carica_tempo];
            }
            else{
                continue;
            }
        }
        carica_picco -= (bin_after+bin_before+1)*ped_signal;
        carica_picco *= (1e5/8)*sampling_time_ns;//10e-3*dt*10e-9/500*1.6*10e-19  //tensioni in mV tempo in secondi resistenza in ohm (500) -> carica in C -> C/e = #elettroni
        peaksdled->Q.at(i_carica)= carica_picco;
        peaks->Q.at(i_carica) = carica_picco;
    }
}



Bool_t GoodCharge( Peaks* p, Short_t index){
    if(index>=p->N) return false; //out of index
    
    else
        if(p->N<2)
        { return true; } //one peak only, good
    
        else
        {
            static Double_t mindeltat = 300e-9;
            if(index == 0){ //first peak, only check right side
                
                if((p->t.at(1)-p->t.at(0))> mindeltat) return true; else return false;
            }
            else
                if( index == (p->N)-1 ) //last peak, only check left side
                    
                {
                    if((p->t.at(index)-p->t.at(index-1))> mindeltat) return true; else return false;
                }
                else
                {
                    if( ( (p->t.at(index)-p->t.at(index-1))> mindeltat) && ((p->t.at(index+1)-p->t.at(index))> mindeltat) ) return true; else return false;
                }
        }
    return false;
}
