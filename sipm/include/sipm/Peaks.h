//
//  Peaks.h
//  
//
//  Created by Luca Tosti on 05/03/2019.
//

#ifndef Peaks_h
#define Peaks_h

typedef struct Peaks{
    Short_t N;
    vector<Int_t> xt;
    vector<Double_t> t;
    vector<Double_t> V;
    vector<Double_t> dt; //from previous peack
    vector<Double_t> Q;
} Peaks;

#endif /* Peaks.h */
