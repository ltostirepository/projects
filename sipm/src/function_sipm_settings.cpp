#include "../include/sipm/general_lib.h"
#include "../include/sipm/Peaks.h"
#include "../include/sipm/functions_sipm_settings.h"

using namespace std;



void ScaleGraph( TGraph* g, Double_t scale){
    for(int ip=0; ip<g->GetN(); ip++){
        g->SetPoint(ip,g->GetX()[ip],scale*g->GetY()[ip]);
        
    }
}



Double_t* GenerateLogBinning(int nbins, float min, float max){
    
    Double_t *X = new Double_t[nbins+1];
    Double_t log_interval = ( TMath::Log10(max) - TMath::Log10(min) ) / nbins;
    for(int i=0; i<=nbins; i++) X[i] = TMath::Power( 10, TMath::Log10(min) + i * log_interval );
    return X;
    
}



Double_t* GenerateBinning(int nbins, float min, float max){
    
    Double_t *X = new Double_t[nbins+1];
    Double_t interval = ( max - min ) / nbins;
    for(int i=0; i<=nbins; i++) X[i] = min + i*interval;
    return X;
    
}
