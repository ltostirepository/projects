//
//  Peaks.h
//  
//
//  Created by Luca Tosti on 05/03/2019.
//

#ifndef functions_sipm_h
#define functions_sipm_h

#include "../include/sipm/general_lib.h"
#include "../include/sipm/Peaks.h"

// settings
void ScaleGraph( TGraph* g, Double_t scale);
Double_t* GenerateLogBinning(int nbins, float min, float max);
Double_t* GenerateBinning(int nbins, float min, float max);

// smoothing and DLED
void TGraphRunningMean( TGraph *g, TGraph **gm, Int_t windowsize);
void gDLEDpoint( const TGraph* gorig, TGraph **gdelayed, Int_t pdelay);
void SetHTime(TH1* h);
TH2F* NormalizeX(TH2F* hh);

//Peaks research
void initpeaks( Peaks *p );
void FindPeaks( TGraph* g, Peaks *peaks, Double_t threshold, bool average, bool fixleftside );
void FindNextPeaks( TGraph *g, Peaks* peaksdled, Peaks *peaks, Double_t twindow, bool average );
void RemovePeaks( TGraph *g, TGraph **gnopeak, Peaks *peaks, Double_t tleft=40e-9, Double_t tright=300e-9);


//Peaks charge
void integra_picco( TGraph *g, Peaks *peaksdled,Peaks *peaks, Double_t ped_signal , Double_t time_before, Double_t time_after, Double_t sampling_time_ns)
Bool_t GoodCharge( Peaks* p, Short_t index);

//useless
Double_t media_colonna(TH2F *histo,int colonna);
void picco_medio(TH2F *histo2d,TGraph *graph);

// PDF printing
void PrintCanvas( TObjArray* arr, string path, string pdfname, bool pdfonly=false );

// principal function
void drawnolaser( const  char* treefilename = "f_in_only_tgraphs.root",int ich=0, bool invert=false, int imin=-999, int imax=-999,bool print=false ,Double_t dt_shift = 8e-10);


#endif /* functions_sipm.h */
