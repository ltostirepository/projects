#include "general_lib.h"
#include "Peaks.h"
#include "function_sipm_search_threshold.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h"



using namespace std;

/*
Double_t Search_threshold(){
    
    Double_t th;

    // definisci variabili dei punti successivi
    
    
    // TEST genera finto DLED per prove
    
    // prendi grafico DLED e costruisci histo punti (no cerca picchi) escludi valori negativi e nulli
    
    
    // costruisci funz prob h_i*(2/V_max^2)*V_i
    
    // crea valori p_i e q_i=1-p_i per ogni tensione V_i = [0;Vmax]
    
    // crea valori ri
    
    // crea R_s = sum(r_i, from N to s) e Q_s= prod(q_N to q_s) con s in[0;N]
    
    // plotta p, q, R, Q
    
    // printa a schermo valore soglia scelta
    
    return th ;
    
}


Double_t get_max_V(TGraph *g){
    
    Double_t max=0;
    Int_t N = g->GetN();
    Double_t *Y = g->GetY();
    
    for(Int_t i=0;i<N;i++){
        if(Y[i]>max){
            max=Y[i];
        }
    }
    return max;
}


void fill_histo(TH1F *h, TGraph *g){
    
    Double_t *X = g->GetX();
    Double_t *Y = g->GetY();
    Int_t N = g->GetN();
    
    for(Int_t i=0; i<N; i++){
        if(Y[i]<=0){continue;}
        h->Fill(Y[i]);
    }
}



Double_t probability_p(Double_t h_i, Double_t V_i, Double_t V_max){

    return h_i*(2/(V_max*V_max))*V_i ;    //<<<<<<<------ devo aggiungere il dV ? non sembra
}



Double_t probability_q(Double_t p_i){
    
    return 1 - p_i ;
    
}



void create_r(Double_t *p, Double_t *q, Int_t N){
    
    for(int i=0;i<N;i++){
        if (q[i]=0){
            cout<<"q_"<<i<<" = 0, fix it in search threshold\n";
            exit(1);
        }
        r[i] = p[i]/q[i];
    }
    
}



void create_R(Double_t *p, Double_t *q, Int_t N, Double_t *R, Double_t *r){
    
    create_r(p,q,N,r);
    
    R[N] = r[n];
    Double_t sum = r[N];
    
    for(int s=N-1;s>-1;s--){
        
        sum += r[s];
        R[s] = sum;
    
    }
    
    
}




void create_Q(Double_t *q, Int_t N, Double_t *Q){
    
    Q[N] = q[N];
    Double_t prod = q[N];
    
    for(int s=N-1;s<-1;s--){
        
        prod = *= r[s];
        Q[s] = prod;
        
    }
    
}

void generate_data(Double_t Vmax){
    
    // genera U[0, 1]
    
    // dividi in sezioni (prob di ogni picco) ; se x € sezione -> aumenta conteggio sezione
    
    // genera conteggi gaussiana = conteggi sezione (picco o rumore)
    
    // mischia posizione picchi (vettore "tempi" -> shuffle -> associa picchi sezioni)
    
    //genera un documento da cui ripescare i dati senza simularli ogni volta
    
}

*/

void testRandom(Int_t N) {
    
    TFile *f = new TFile("test_random.root","RECREATE");
    
    TRandom *r1=new TRandom();
    TRandom2 *r2=new TRandom2();
    TRandom3 *r3=new TRandom3();
    
    TH1F *h1=new TH1F("h1","TRandom",500,0,1);
    TH1F *h2=new TH1F("h2","TRandom2",500,0,1);
    TH1F *h3=new TH1F("h3","TRandom3",500,0,1);
    
    
    for (Int_t i=0; i<N; i++) { h1->Fill(r1->Uniform(0,1)); }
    cout << "Random: done" << endl;
    
    
    for (Int_t i=0; i<N; i++) { h2->Fill(r2->Uniform(0,1)); }
    cout << "Random2: done"<< endl;
    
    
    for (Int_t i=0; i<N; i++) { h3->Fill(r3->Uniform(0,1)); }
    cout << "Random3: done" << endl;
    
    f->Write();
    
    cout<<"END TEST_RANDOM\n";
}
