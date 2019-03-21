//
//  print_and_trials
//  
//
//  Created by Luca Tosti on 05/03/2019.
//

#ifndef functions_sipm_print_and_trials_h
#define functions_sipm_print_and_trials_h



//useless trials
Double_t media_colonna(TH2F *histo,int colonna);
void picco_medio(TH2F *histo2d,TGraph *graph);

// PDF printing
void PrintCanvas( TObjArray* arr, std::string path, std::string pdfname, bool pdfonly=false );


#endif /* functions_sipm_print_and_trials.h */
