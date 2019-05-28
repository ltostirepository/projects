#include "general_lib.h"
#include "Peaks.h"
#include "function_sipm_search_threshold.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h"



using namespace std;


Double_t Search_threshold(const  char* file){
 
    TFile *f = new TFile("test_threshold.root","RECREATE");
    TFile A(file);
    TH1D *h = (TH1D*)A.Get("h1");
    
    Double_t th = 0;
    
    // definisci variabili dei punti successivi
    
    
    
 // ####### (per ora prendi direttamente l'histo) prendi grafico DLED e costruisci histo punti (no cerca picchi) escludi valori negativi e nulli
    Int_t entries_h=h->GetEntries();
    Int_t entries_h_1= 0;
    entries_h_1= entries_h +1;
    
    
    Double_t p[501]; // imposta dim dall'histo
    Double_t q[501]; // imposta dim dall'histo
    Double_t r[501]; // imposta dim dall'histo
    
    Double_t R[501]; // imposta dim dall'histo
    Double_t Q[501]; // imposta dim dall'histo
    
    Double_t V_histo[501]; // imposta dim dall'histo
    Double_t Vmax = h->GetBinCenter(500);
    // costruisci funz prob h_i*(2/V_max^2)*V_i
    /*Double_t p[entries_h_1]; // imposta dim dall'histo
    Double_t q[entries_h_1]; // imposta dim dall'histo
    Double_t r[entries_h_1]; // imposta dim dall'histo
    cout<<"in4\n";
    Double_t R[entries_h_1]; // imposta dim dall'histo
    Double_t Q[entries_h_1]; // imposta dim dall'histo
    cout<<"in5\n";
    Double_t V_histo[entries_h_1]; // imposta dim dall'histo
    Double_t Vmax = h->GetBinCenter(entries_h); //così V_max dipende dall'histo (devi scorrere il vettore e trovare il più grande utile)
    */
    
    //Double_t scale = 1.0/h->Integral();
    //h->Scale(scale);// normalizzazione histo
    for (int i=1; i<502; i++){ //parto da 1 per non contare l'indice dell'underflow (che comunque indicherebbe un bin vuoto)
        //arrivo a 502 per contare l'overflow
        
        // crea valori p_i e q_i=1-p_i per ogni tensione V_i = [0;Vmax]
        V_histo[i] = h->GetBinCenter((Int_t)i);
        p[i] = probability_p(h,(Int_t)i,h->GetBinContent((Int_t)i), h->GetBinCenter((Int_t)i) , Vmax);//TH1D *h,Int_t i, Double_t V_i, Double_t V_max
        q[i] = probability_q(p[i]);
        
    }
    
    // crea R_s = sum(r_i, from N to s) e Q_s= prod(q_N to q_s) con s in[0;N]
    create_R(p,q,500,R,r);
    //cout<<"erererer\n";
    create_Q(q,500,Q,r);
    
    
    // plotta p, q, R, Q
    TGraph *gp = new TGraph(501);
    gp->SetName("gp");
    TGraph *gq = new TGraph(501);
    gq->SetName("gq");
    TGraph *gR = new TGraph(501);
    gR->SetName("gR");
    TGraph *gQ = new TGraph(501);
    gQ->SetName("gQ");
    TGraph *gRQ = new TGraph(501);
    gRQ->SetName("gRQ");
    
    
    for(Int_t d = 0; d<501;d++){
        
        gp->SetPoint(d, V_histo[d], p[d]);
        gq->SetPoint(d, V_histo[d], q[d]);
        gR->SetPoint(d, V_histo[d], R[d]);
        gQ->SetPoint(d, V_histo[d], Q[d]);
        gRQ->SetPoint(d, V_histo[d], R[d]*Q[d]);

    }
    
    
 
    f->Write();
    f->WriteTObject(h);
    f->WriteTObject(gp);
    f->WriteTObject(gq);
    f->WriteTObject(gR);
    f->WriteTObject(gQ);
    f->WriteTObject(gRQ);
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
    
    //Double_t *X = g->GetX();
    Double_t *Y = g->GetY();
    Int_t N = g->GetN();
    
    for(Int_t i=0; i<N; i++){
        if(Y[i]<=0){continue;}
        h->Fill(Y[i]);
    }
}



Double_t probability_p(TH1D *h,Int_t i,Double_t h_i, Double_t V_i, Double_t V_max){

    double integral = h->Integral(i,501);
    //return (integral*h_i*(2/(V_max*V_max))*V_i)/(h->GetEntries()) ;//<<<<<<<------ devo aggiungere il dV ? non sembra
    return integral/(h->GetEntries());
}



Double_t probability_q(Double_t p_i){
    
    return 1 - p_i ;
    
}



void create_r(Double_t *p, Double_t *q, Int_t N, Double_t *r){
    
    for(int i=0;i<N;i++){
        if (q[i]==0){
            cout<<"q_"<<i<<" = 0, fix it in search threshold\n";
            r[i] = 1;
            q[i] = 0.1;
            //exit(1);
        }
        else{r[i] = p[i]/q[i];}
        
    }
    
}



void create_R(Double_t *p, Double_t *q, Int_t N, Double_t *R, Double_t *r){
    
    create_r(p,q,N,r);
    
    R[N] = r[N];
    Double_t sum = r[N];
    
    for(int s=N-1;s>-1;s--){
        
        sum += r[s];
        R[s] = sum;
    
    }
    
    
}




void create_Q(Double_t *q, Int_t N, Double_t *Q, Double_t *r ){
    
    Q[N] = q[N];
    Double_t prod = q[N];
    
    if (prod == 0 || q[N]==0 ){cout<<"AAAAHHHHHHHHHHHHH!!!!!\n";}
    
    for(int s=N-1;s>(-1);s--){
        
        prod *= q[s];
        Q[s] = prod;
        
    }
    
}


Double_t sigma_2_extimation(Double_t V, Double_t E_x,Double_t E_x_2,Double_t E_x_3){
    Double_t sigma_2;
    
    sigma_2=((pow(V,3)*E_x) + (3*pow(V,2)*pow(E_x,2)) + ((pow(E_x,3)-E_x_3)*V))/(3*E_x_2);
    
    return sigma_2;
}

int factorial(int n){
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

Double_t probability(Double_t lambda, int k,Double_t alpha0, Double_t ct){
    
    Double_t a = (Double_t)(exp(-((double)lambda))*pow(((double)lambda),k)/(double)factorial(k));
    
    if(k==0){
        a*=alpha0;
    }
    else if(k==1){
        a*=1;
    }
    else{
        //a*=alpha(lambda,ct);
        //cout<<"alpha "<<alpha(lambda,ct)<<"\n";
        a*= (1+ct) ;
        
    }
    
    //cout<<"N = "<<norm(lambda, alpha0, ct)<<"\n";
    a = a/norm(lambda, alpha0, ct);
    
    return a;
    
}

Double_t norm(Double_t lambda, Double_t alpha0, Double_t ct){
    
    Double_t N = exp(-((double)lambda))*( alpha0 + lambda )+ (/*alpha(lambda,ct)*/(1+ct))*(1-exp(-((double)lambda))-(exp(-((double)lambda))*lambda));
    
    return N;
}

Double_t alpha(Double_t lambda,Double_t ct){
    
    Double_t e = (Double_t)exp(-((double)lambda));
    Double_t a = (ct*e)/((1-e-e*lambda)*(1-ct));
    
    return a;
    
}


void generate_data(Double_t lambda, Double_t V_1,Double_t sigma_sgn, Double_t sigma_noise,Double_t alpha0, Double_t ct, Int_t N){
    
    TFile *f = new TFile("test_random_generate_data.root","RECREATE");
    // genera U[0, 1]
    int N_picchi = 10;
    cout<<" average value (if negligible noise ) "<<lambda*V_1<<"\n";
    Double_t sections_prob_value[N_picchi];
    Double_t sections_th[N_picchi];
    Int_t N_histo=500;
    int sections[N_picchi];
    int sum = 0;
    
    TRandom3 *r1=new TRandom3();
    Double_t r=0;
    TH1D *h = new TH1D("h1","TRandom",N_histo,0,V_1*N_picchi);
    
    
    Int_t maxxx=0;
    
    Double_t normaliz=0;
    for(int i=0; i<100; i++){
        if(probability(lambda,i, alpha0, ct)<1e-20){cout<<"negligible prob , stop at "<<i<<" peak\n"; maxxx=(Int_t)i ;break;}
        normaliz += probability(lambda,i, alpha0, ct);
        cout<<"probability "<<i<<" "<<probability(lambda,i, alpha0, ct)<<"\n";
    }
    cout<<"normalization ----> "<<normaliz<<"\n";
    
    
    
    for(int i=0; i<N_picchi-1; i++){
        
        sections_prob_value[i] = probability(lambda,i, alpha0, ct);
        if (i==0){
            sections_prob_value[N_picchi-1] = probability(lambda,N_picchi-1,alpha0, ct);
            sections_th[N_picchi-1]=1;
            sections_th[0] = sections_prob_value[0];
            
        }
        else{
            sections_th[i] = sections_prob_value[i] + sections_th[i-1];
        }
        
    }
    
    for(int i=0;i<N_picchi;i++){
        sections[i]=0;
    }
    
    for(int i=0;i<N;i++){
        r=r1->Uniform(0,1);
        
        for(int k=0;k<N_picchi;k++){
            if(r<sections_th[k]){
                sections[k] += 1 ;
                break;
            }
        }
    }
    
    Double_t conteggio=0;
    for(int k=0;k<N_picchi;k++){
        conteggio += sections[k];
    }
    cout<<"conteggio ---> "<< conteggio <<"\n";
    
    if(sections[N_picchi-1]>N*sections_prob_value[N_picchi-1]*(1+0.05)){
        
        cout<<"\n\nwarning: too much signal in last gaussian function -> renormalization last peak\n";
        sections[N_picchi-1]=N*sections_prob_value[N_picchi-1];
        for(int i=0;i<N_picchi;i++){
            sum+=sections[i];
        }
        cout<<"total number of events "<< sum <<" / "<<N<<" -> "<< ((N-sum)/(double)N)*100 <<"% error in distribution \n";
        cout<<"NEW number of events : "<< sum <<" ; "<< (N-sum) <<" events lost\n\n";
        h->SetBinContent(N_histo+1,(Double_t)(N-sum));
        
    }
    else{cout<<"no cut on peaks -> N=sum\n\n";sum=N;}
    
    Double_t sigma;
    /*
    Double_t E_x=0;
    Double_t E_x_2=0;
    Double_t E_x_3=0;
    Double_t Var_x=0;
    Double_t k_var=0;
    Double_t V_plus=0;
    Double_t V_minus=0;
    Double_t sigma_extim_plus=0;
    Double_t sigma_extim_minus=0;
    Double_t lambda_extim_plus=0;
    Double_t lambda_extim_minus=0;
    */
    //Double_t conta=1;
    for(int i = 0; i<N_picchi;i++){
        
        if (i==0){
            sigma=sigma_noise;
        }
        else{
            sigma=sigma_sgn;
        }
        
        cout<<"i*V_1,sigma "<<i*V_1<<" "<<sigma<<"\n";
        cout<<"sections"<<"["<<i<<"]"<< sections[i]<<"\n";
        cout<<"sections_prob_value"<<"["<<i<<"]"<< sections_prob_value[i]<<"\n";
        cout<<"sections_th"<<"["<<i<<"]"<< sections_th[i]<<"\n";
        
        
        
        for(int k = 0;k<sections[i];k++){
            
            r=0.0;
            
            if(i==0){
                r=r1->Gaus(0,sigma);
            }
            else if(i==1){
                r=r1->Gaus(V_1,sigma);
                r+= r1->Gaus(0,sigma_noise);
                
            }
            else{
                for(int j=0;j<i;j++){
                    r+=r1->Gaus(V_1,sigma);
                    r+= r1->Gaus(0,sigma_noise);
                }
            }
            
            
            if(r<0){
                k--;
            }
            else{
                //cout<<"r_1 "<<r<<"\n";
                h->Fill(r);
                //E_x = (conta-1/conta)*E_x + (r/conta);
                //E_x_2 = (conta-1/conta)*E_x_2 + (r*r/conta);
                //E_x_3 = (conta-1/conta)*E_x_3 + (r*r*r/conta);
                //conta+=1;
                //cout<<"r_2 "<<r<<"\n";
            }
        }
    }
    
    //E_x/=sum;
    //E_x_2/=sum;
    //E_x_3/=sum;
    /*cout<<"conta "<<conta<<" sum "<<sum<<"\n";
    
    TGraph *g = new TGraph(N_histo);
    g->SetName("g");
    
    Double_t V_d=0;
    Double_t s =0;
    
    for(Int_t d = 0; d<N_histo*10;d++){
    
        V_d = ((V_1*N_picchi)/(N_histo*10))*d;
        
        s=sigma_2_extimation(V_d,E_x,E_x_2,E_x_3);
        if (s>0){s= pow(s,0.5);}
        else{s=0;}
        
        g->SetPoint(d, V_d ,s);
                
    }
    */
    TGraph *g = new TGraph(maxxx);
    g->SetName("g");
    Double_t s =0;
    
    for(Int_t d = 0; d<maxxx;d++){
        
        s+=probability(lambda, (int)d ,alpha0, ct);
        g->SetPoint(d, (Double_t)d ,s);
        
    }
    //////result Nan -> check it
    /*
    cout<<"E_x "<<E_x<<"\n";
    cout<<"E_x_2 "<<E_x_2<<"\n";
    cout<<"E_x_3 "<<E_x_3<<"\n";
    
    
    Var_x=E_x_2 - (E_x*E_x);
    k_var= ((E_x*E_x*E_x) - E_x_3) + (3*E_x_2*((Var_x/E_x)+E_x));
    
    cout<<"Var_x "<<Var_x<<"\n";
    cout<<"k_var "<<k_var<<"\n\n";
    
    V_plus=(3*Var_x+pow((9*Var_x*Var_x-4*E_x*k_var),0.5))/(2*E_x);
    V_minus=(3*Var_x-pow((9*Var_x*Var_x-4*E_x*k_var),0.5))/(2*E_x);
    
    sigma_extim_plus=V_plus*((Var_x/E_x)+E_x) - V_plus*V_plus ;
    sigma_extim_minus= V_minus*((Var_x/E_x)+E_x) - V_minus*V_minus ;
    lambda_extim_plus=E_x/V_plus;
    lambda_extim_minus=E_x/V_minus;
    
    cout<<"V "<<V_1<<"\n";
    cout<<"V_plus "<<V_plus<<"\n";
    cout<<"V_minus "<<V_minus<<"\n\n";
    
    cout<<"sigma_sgn "<<sigma_sgn*sigma_sgn<<"\n";
    cout<<"sigma_extim_plus "<<sigma_extim_plus<<"\n";
    cout<<"sigma_extim_minus "<<sigma_extim_minus<<"\n\n";
    
    
    cout<<"lambda "<<lambda<<"\n";
    cout<<"lambda_extim_plus "<<lambda_extim_plus<<"\n";
    cout<<"lambda_extim_minus "<<lambda_extim_minus<<"\n\n";
     */
    
    g->Write();
    f->Write();
    //f->Close();
    cout<<"END CREATE_DATA\n";
    
    // dividi in sezioni (prob di ogni picco) ; se x € sezione -> aumenta conteggio sezione
    
    // genera conteggi gaussiana = conteggi sezione (picco o rumore)
    
    //ho generato direttamente l'histo salvato in test_random_generate_data.root
    
}



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
