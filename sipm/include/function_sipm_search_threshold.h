

#ifndef drawnolaser_h
#define drawnolaser_h


#include "general_lib.h"
#include "Peaks.h"
#include "function_sipm_search_threshold.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TVirtualFFT.h"
#include "TAxis.h"
#include "TH1.h"
#include "TArrayD.h"

using namespace std;


Double_t Search_threshold(const  char* file);

Double_t get_max_V(TGraph *g);

void fill_histo(TH1F *h, TGraph *g);

Double_t probability_p(TH1D *h,Int_t i,Double_t h_i, Double_t V_i, Double_t V_max);

Double_t probability_q(Double_t p_i);

void create_r(Double_t *p, Double_t *q, Int_t N);

void create_R(Double_t *p, Double_t *q, Int_t N, Double_t *R, Double_t *r, Double_t *e);

void create_Q(Double_t *q, Int_t N, Double_t *Q,Double_t *r);

Double_t alpha(Double_t lambda,Double_t ct);

Double_t payoff(TH1D *h,Int_t i,Double_t h_i, Double_t V_i, Double_t V_max);

Double_t sigma_2_extimation(Double_t V, Double_t E_x,Double_t E_x_2,Double_t E_x_3);
int factorial(int n);
Double_t probability(Double_t lambda, int k,Double_t alpha0, Double_t ct);
Double_t norm(Double_t lambda, Double_t alpha0, Double_t ct);
void generate_data(Double_t lambda, Double_t V_1,Double_t sigma_sgn, Double_t sigma_noise, Double_t alpha0, Double_t ct, Int_t N);
void testRandom(Int_t nrEvents=500000000);

Double_t ScaleX(Double_t x,Double_t min, Double_t max);

Double_t ScaleY(Double_t y, Double_t N);

void ScaleAxis(TAxis *a, Double_t (*Scale)(Double_t));

void ScaleXaxis(TH1 *h, Double_t (*Scale)(Double_t));

void ScaleYaxis(TH1 *h, Double_t (*Scale)(Double_t));






#endif
