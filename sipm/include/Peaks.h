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
    std::vector<Int_t> xt;
    std::vector<Double_t> t;
    std::vector<Double_t> V;
    std::vector<Double_t> dt; //from previous peak
    std::vector<Double_t> Q;
} Peaks;

#endif /* Peaks.h */
