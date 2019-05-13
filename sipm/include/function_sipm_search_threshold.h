

#ifndef drawnolaser_h
#define drawnolaser_h


#include "general_lib.h"
#include "Peaks.h"
#include "function_sipm_search_threshold.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h"


using namespace std;

/*
 Double_t Search_threshold();

 Double_t get_max_V(TGraph *g);

 void fill_histo(TH1F *h, TGraph *g);

 Double_t probability_p(Double_t h_i, Double_t V_i, Double_t V_max);

 Double_t probability_q(Double_t p_i);

 void create_r(Double_t *p, Double_t *q, Int_t N);

 void create_R(Double_t *p, Double_t *q, Int_t N, Double_t *R, Double_t *r);

 void create_Q(Double_t *q, Int_t N, Double_t *Q);

 void generate_data(Double_t Vmax);
*/

int factorial(int n);
Double_t poisson(Double_t lambda, int k);
void generate_data(Double_t lambda, Double_t V_1,Double_t sigma_sgn, Double_t sigma_noise, Int_t N);
void testRandom(Int_t nrEvents=500000000);

#endif
