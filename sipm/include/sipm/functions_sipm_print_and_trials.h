//
//  print_and_trials
//  
//
//  Created by Luca Tosti on 05/03/2019.
//

#ifndef functions_sipm_print_and_trials_h
#define functions_sipm_print_and_trials_h

#include "../include/sipm/general_lib.h"
#include "../include/sipm/Peaks.h"

//useless trials
Double_t media_colonna(TH2F *histo,int colonna);
void picco_medio(TH2F *histo2d,TGraph *graph);

// PDF printing
void PrintCanvas( TObjArray* arr, string path, string pdfname, bool pdfonly=false );

// principal function
void drawnolaser( const  char* treefilename = "f_in_only_tgraphs.root",int ich=0, bool invert=false, int imin=-999, int imax=-999,bool print=false ,Double_t dt_shift = 8e-10);


#endif /* functions_sipm_print_and_trials.h */
