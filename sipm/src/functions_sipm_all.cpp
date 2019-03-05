#include "../include/sipm/general_lib.h"
#include "../include/sipm/Peaks.h"
#include "../include/sipm/functions_sipm_all.h"


using namespace std;


//settings

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


//smoothing and DLED
void TGraphRunningMean( TGraph *g, TGraph **gm, Int_t windowsize){
    
    (*gm) = new TGraph();
    
    Int_t nn=g->GetN();
    vector<Double_t> xx;
    for(int i=0;i<nn;i++) xx.push_back( g->GetX()[i] );
    vector<Double_t> yy;
    for(int i=0;i<nn;i++) yy.push_back( g->GetY()[i] );
    vector<Double_t> xxtemp;
    vector<Double_t> yytemp;
    
    vector<Double_t> v;
    Int_t vsize=windowsize;
    Int_t iv=vsize-1;
    Int_t ix=(iv/2);
    
    Double_t sum=0, mean=0;
    v.push_back(0); //dummy, to synchronize with the algorithm
    for(int i=0; i<(iv); i++) { v.push_back(yy[i]); sum+=yy[i]; }
    while( iv < (int)xx.size() ){
        sum -= v.at(0);
        v.erase( v.begin() );
        v.push_back( yy[iv] );
        sum += yy[iv];
        mean = sum / vsize;
        xxtemp.push_back(xx[ix]);
        yytemp.push_back(mean);
        iv++; ix++;
    }
    //if(debug) cout<< xxtemp.size() <<" "<< yytemp.size() << endl;
    
    for(int ii=0; ii<(int)xxtemp.size(); ii++) (*gm)->SetPoint(ii,xxtemp[ii],yytemp[ii]);
    //if(debug) cout<<"(*gm)->GetN() "<<(*gm)->GetN()<<endl;
    xxtemp.clear();
    yytemp.clear();
    xx.clear();
    yy.clear();
    v.clear();
}




void gDLEDpoint( const TGraph* gorig, TGraph **gdelayed, Int_t pdelay){
    
    (*gdelayed) = new TGraph(); (*gdelayed)->SetName( Form("%s_delay",gorig->GetName()) );
    Double_t *xaxis = gorig->GetX();
    Double_t *yaxis = gorig->GetY();
    //if(debug)  cout<<"xaxis "<< xaxis<< " gorig->GetX()"<< gorig->GetX()<<" gorig "<< gorig <<endl;
    //if(debug)  cout<<"yaxis "<< yaxis<< " gorig->GetY()"<< gorig->GetY()<<" gorig "<< gorig <<endl;
    Double_t minx = xaxis[0];
    Double_t maxx = xaxis[gorig->GetN()-1];
    for(int ip=0; ip<=gorig->GetN(); ip++){
        if( ip+pdelay<0 || ip+pdelay>gorig->GetN()-1) continue;
        Double_t x=xaxis[ip];
        Double_t xdel=xaxis[ip+pdelay];
        if( xdel<minx || xdel>maxx ) continue;
        else
        {
            Double_t toset = yaxis[ip] - yaxis[ip+pdelay];
            (*gdelayed)->SetPoint( (*gdelayed)->GetN(), x, toset );
        }
    }
}



void SetHTime(TH1* h){
    h->GetXaxis()->SetTimeFormat("%H:%m%F1970-01-01 00:00:00");
    h->GetXaxis()->SetTimeDisplay(true);
    return;
}



TH2F* NormalizeX(TH2F* hh){     //Get pointer to TH2F normalized in X slices
    std::string hname; hname.assign( Form("%s_normX", hh->GetName() ) );
    if( gDirectory->Get( hname.c_str() ) )
    {
        std::cout<<Form("NormalizeX::Deleting old %s",hname.c_str() )<<std::endl;
        gDirectory->Get( hname.c_str() ) ->Delete();
    }
    TH2F *hNorm = (TH2F*)hh->Clone( Form("%s_normX", hh->GetName() ) );
    hNorm->Sumw2();
    Double_t nentries;
    Double_t value;
    for(int ibinx=1; ibinx<=hNorm->GetNbinsX(); ibinx++)
    {
        nentries=0;
        for(int ibiny=1; ibiny<=hNorm->GetNbinsY(); ibiny++)
        {
            nentries+=hNorm->GetBinContent(ibinx,ibiny);
        }
        for(int ibiny=1; ibiny<=hNorm->GetNbinsY(); ibiny++)
        {
            if( nentries<=0 ) continue;
            value = hNorm->GetBinContent(ibinx,ibiny);
            value /= nentries;
            hNorm->SetBinContent(ibinx,ibiny,value);
        }
    }
    return hNorm;
}


//Peaks reseach
void initpeaks( Peaks *p ){
    p->N=0;
    p->t.clear();
    p->V.clear();
    p->dt.clear();
    p->xt.clear();
    p->Q.clear();
}



void FindPeaks( TGraph* g, Peaks *peaks, Double_t threshold, bool average, bool fixleftside ){
    initpeaks(peaks);
    Int_t ip=-1;
    Double_t *Y = g->GetY();
    Double_t *X = g->GetX();
    while( (ip+1)<g->GetN() ){
        ip++;
        if( Y[ip]<threshold ){
            continue;
        }
        else{
            //first above threshold. Find when it goes below again
            Int_t minip=ip;
            while((ip+1)<g->GetN() && Y[ip]>threshold  ){ ip++;}
            if((ip+1)-(int)g->GetN()>0){continue;}
            Int_t maxip=ip;
            //find local maximum
            Double_t localVmax=-999, localtmax=-999;
            Int_t localxmax=-999;
            for( Int_t jj=minip; jj<maxip; jj++ ){
                if( fabs(Y[jj])>localVmax ) {
                    localVmax=Y[jj];
                    localtmax=X[jj];
                    localxmax=jj;
                }
            }
            if(average)
            {
                int nn=1;
                if( (localxmax-1)>=0 )        { localVmax += fabs(Y[localxmax-1]); nn++; }
                if( (localxmax+1)<g->GetN() ) { localVmax += fabs(Y[localxmax+1]); nn++; }
                localVmax /= nn;
            }
            //if(debug) cout<<"prova2\n";
            int xminside=localxmax-1;
            while( xminside>1 && Y[xminside]>Y[xminside-1] ) xminside--;
            //while( xminside>1 && Y[xminside]>Y[xminside-2] ) xminside--;//aggiunto io
            Double_t minside=Y[xminside];
            int xmaxside=localxmax+11;//+11
            while( xmaxside<g->GetN()-1 && Y[xmaxside]>Y[xmaxside+1] ) xmaxside++;
            Double_t maxside=Y[xmaxside];
            Double_t ave = 0.5*(minside+maxside);
            //if(debug) cout<<"prova3\n";
            //fill peaks structure
            peaks->N++;
            peaks->xt.push_back(localxmax);
            peaks->t.push_back(localtmax);
            if( fixleftside ) { peaks->V.push_back(localVmax-minside); }//remove only left side to avoi bias introduce by undershoot on the right
            else { peaks->V.push_back(localVmax); }
            peaks->Q.push_back(-666);
            if( peaks->N>1 ){ peaks->dt.push_back( fabs(localtmax - peaks->t.at( peaks->N-2 ) )); }
            else peaks->dt.push_back(-999);//valore fisso per non confondersi durante la rimozione
        }
    }
    printf("%s::%d-peaks-found\n",__FUNCTION__,peaks->N);
    return;
}


void FindNextPeaks( TGraph *g, Peaks* peaksdled, Peaks *peaks, Double_t twindow, bool average ){
    initpeaks(peaks);
    Double_t *Y = g->GetY();
    Double_t *X = g->GetX();
    //Bool_t *found = new Bool_t[g->GetN()]; for(int i=0; i<g->GetN(); i++) found[i]=false;
    Bool_t found[g->GetN()]; for(int i=0; i<g->GetN(); i++) found[i]=false;
    
    //if(debug) cout<<"peaksdled->N [FindNextPeaks]"<<peaksdled->N<<endl;
    for(int ip=0; ip<peaksdled->N; ip++)
    {
        //if(debug) {if(ip%1000==0) printf("ip= %d N= %d\n", ip, peaksdled->N);}
        Double_t localVmax=-999, localtmax=-999;
        Int_t localxmax = peaksdled->xt.at(ip);
        Int_t ix = peaksdled->xt.at(ip);
        Double_t deltat=0;
        while( ix<g->GetN() && deltat<twindow ){
            if( Y[ix] >localVmax && found[ix]==false) {
                localVmax=Y[ix];
                localtmax=X[ix];
                localxmax=ix;
            }
            deltat = (X[ix]-peaksdled->t.at(ip));
            ix++;
        }
        if(average)
        {
            int nn=1;
            if( (localxmax-1)>=0 )        { localVmax += fabs(Y[localxmax-1]); nn++; }
            if( (localxmax+1)<g->GetN() ) { localVmax += fabs(Y[localxmax+1]); nn++; }
            localVmax /= nn;
        }
        //fill peaks structure
        //cout<<"fill peak structure\n";
        peaks->N++;
        peaks->xt.push_back(localxmax);
        peaks->t.push_back(localtmax);
        peaks->V.push_back(localVmax);
        peaks->Q.push_back(-666);
        if( peaks->N>1 ){ peaks->dt.push_back( fabs(localtmax - peaks->t.at( peaks->N-2 ) )); }
        else peaks->dt.push_back(-999);//valore fisso per facilitare la rimozione
        
        found[localxmax]=true;
    }
    
    
}



void RemovePeaks( TGraph *g, TGraph **gnopeak, Peaks *peaks, Double_t tleft, Double_t tright ){
    
    Int_t n_intt = g->GetN();
    //if (debug) cout<<"RemovePeaks - 0.0\n";
    int n = (int)n_intt;
    //if (debug) cout<<"n = "<<n<<" n_intt = "<<n_intt<<"\n";
    //if (debug) {cout<<"try fake n  -> 50 instead of "<<n<<endl; n=50;}
    //if (debug) cout<<"cast done\n";
    Double_t* x = new Double_t[n];/////n is too big for a normal allocation of an array , now it is in dynamic location
    Double_t* y = new Double_t[n];
    //if (debug) cout<<"before loop\n";
    for(int i=0; i<n; i++) {
        x[i]=g->GetX()[i];
        y[i]=g->GetY()[i];
    }
    //if (debug) cout<<"RemovePeaks - 0.1\n";
    Double_t dx = (y[n-1]-y[0])/(n-1);
    //if (debug) cout<<"RemovePeaks - 0.2\n";
    vector<Double_t> minxtoremove, maxxtoremove;
    minxtoremove.push_back(-1); maxxtoremove.push_back(tright); //always remove the first tright interval to avoid tails from peaks
    //if (debug) cout<<"RemovePeaks - 0.3\n";
    for(int ip=0; ip<peaks->N; ip++){
        minxtoremove.push_back(peaks->t.at(ip)-tleft);
        maxxtoremove.push_back(peaks->t.at(ip)+tright);
    }
    //for(int i=0; i<(int)minxtoremove.size();i++){ cout<<Form("%.3e\t%.3e\n", minxtoremove.at(i), maxxtoremove.at(i)); }
    //if (debug) cout<<"RemovePeaks - 1\n";
    
    Int_t thispeak=0;
    bool jumppeak=true;
    static vector<Double_t> newx, newy;
    Int_t newn=0;
    for(int ii=0; ii<n; ii++)
    {
        if(x[ii]>=minxtoremove.at(thispeak) && x[ii]<=maxxtoremove.at(thispeak)) { jumppeak=true; continue; }
        else{
            if(jumppeak) {  if( (thispeak+1)<((int)minxtoremove.size())){thispeak++;} jumppeak=false;}  //the second if is to protect the case in which the last part of the waveform has to \be kept, so thispeak is never bigger than peaks->N
            else // we skip one point. nevermind, to be conservative. otherwise there is conflicts with neighbouring intervals
            {
                newx.push_back(x[ii]); newy.push_back(y[ii]); newn++;
            }
        }
    }
    //if (debug) cout<<"RemovePeaks - 2\n";
    (*gnopeak) = new TGraph(); (*gnopeak)->SetNameTitle( Form("%s_nopeak",g->GetName()), g->GetTitle() );
    for(int i=0; i<newn;i++) {  (*gnopeak)->SetPoint(i,newx.at(i),newy.at(i)); }
    newx.clear();
    newy.clear();
    //if (debug) cout<<"RemovePeaks - 3\n";
    delete[] x;
    delete[] y;
}





//charge
void integra_picco( TGraph *g, Peaks *peaksdled,Peaks *peaks, Double_t ped_signal , Double_t time_before, Double_t time_after, Double_t sampling_time_ns){
    for(int i_carica = 0 ; i_carica < peaksdled->N;i_carica++){
        Double_t carica_picco=0;
        int bin_before = (int)(time_before/sampling_time_ns);
        int bin_after =  (int)(time_after/sampling_time_ns);
        Double_t tempo = peaksdled->t.at(i_carica);
        //int bin_picco = (int)(tempo/sampling_time_ns);
        int bin_picco = peaksdled->xt.at(i_carica);
        //int running_bin;
        for(int i_carica_tempo = -bin_before ; i_carica_tempo < bin_after  ; i_carica_tempo++){
            if(bin_picco+i_carica_tempo>=0 && bin_picco+i_carica_tempo<g->GetN()){
                carica_picco += g->GetY()[bin_picco+i_carica_tempo];
            }
            else{
                continue;
            }
        }
        carica_picco -= (bin_after+bin_before+1)*ped_signal;
        carica_picco *= (1e5/8)*sampling_time_ns;//10e-3*dt*10e-9/500*1.6*10e-19  //tensioni in mV tempo in secondi resistenza in ohm (500) -> carica in C -> C/e = #elettroni
        peaksdled->Q.at(i_carica)= carica_picco;
        peaks->Q.at(i_carica) = carica_picco;
    }
}



Bool_t GoodCharge( Peaks* p, Short_t index){
    if(index>=p->N) return false; //out of index
    
    else
        if(p->N<2)
        { return true; } //one peak only, good
    
        else
        {
            static Double_t mindeltat = 300e-9;
            if(index == 0){ //first peak, only check right side
                
                if((p->t.at(1)-p->t.at(0))> mindeltat) return true; else return false;
            }
            else
                if( index == (p->N)-1 ) //last peak, only check left side
                    
                {
                    if((p->t.at(index)-p->t.at(index-1))> mindeltat) return true; else return false;
                }
                else
                {
                    if( ( (p->t.at(index)-p->t.at(index-1))> mindeltat) && ((p->t.at(index+1)-p->t.at(index))> mindeltat) ) return true; else return false;
                }
        }
    return false;
}



//useless trials
Double_t media_colonna(TH2F *histo,int colonna){
    Double_t media=0;
    int cont=0;
    for(int i =1;i<histo->GetNbinsY();i++){
        cont += histo->GetBinContent(colonna,i);
        media += (histo->GetYaxis()->GetBinCenter(i))*(histo->GetBinContent(colonna,i));
    }
    media = (Double_t)(media/(Double_t)cont);
    return media;
}



void picco_medio(TH2F *histo2d,TGraph *graph){
    Double_t res=0;
    for(int i=0; i<histo2d->GetNbinsX();i++){
        res = media_colonna(histo2d,i);
        graph->SetPoint(i,histo2d->GetXaxis()->GetBinCenter(i),res);
    }
}


// PDF printing


void PrintCanvas( TObjArray* arr, string path, string pdfname, bool pdfonly ){
    if( !arr) return;
    TCanvas *c = new TCanvas();
    c->SaveAs( Form("%s/%s[", path.c_str(), pdfname.c_str()) );
    for(int ic=0; ic<(int)arr->GetEntries(); ic++)
    {
        c = (TCanvas*)(arr->At(ic) );
        if( !c ) continue;
        c->SaveAs( Form("%s/%s.pdf", path.c_str(), c->GetName()  ) );
        if( !pdfonly){
            c->SaveAs( Form("%s/%s.root", path.c_str(), c->GetName()  ) );
            c->SaveAs( Form("%s/%s.C", path.c_str(), c->GetName()  ) );
            c->SaveAs( Form("%s/%s.eps", path.c_str(), c->GetName()  ) );
        }
        c->SaveAs( Form("%s/%s", path.c_str(), pdfname.c_str()) );
    }
    c = new TCanvas();
    c->SaveAs( Form("%s/%s]", path.c_str(), pdfname.c_str()) );
    return;
}


// principal function (se non serve più togli)
/*
 void drawnolaser( const  char* treefilename = "f_in_only_tgraphs.root", int ich=0, bool invert=false, int imin=-999, int imax=-999, bool print=false ,Double_t dt_shift = 8e-10){
 
 cout<<"........START PROGRAM........\n";
 int status;
 Double_t peakthreshold = 0.004;//0.0001;//V
 const static Double_t sipmarea = 36; //mm2
 TObjArray *arrCanvas = new TObjArray();
 TObjArray *arrwfCanvas = new TObjArray();
 
 
 //////////////////attento , sto automatizzando i range -> potrebbero sorgere casini
 
 Double_t vaxismin =-5;//-10;// 0;//mV
 Double_t vaxismax = 50;//100;//25.0;//mV
 Double_t vaxisdledmin = -5;//mV
 Double_t vaxisdledmax = 50;//mV
 Double_t mintime =  -5e5;//ns -500us
 Double_t maxtime =   5e5;//ns  500us
 int run_mean_points = 14;//14
 int point_shift_dled = 7;//1
 
 if (debug) printf("line: %d\n",__LINE__);
 
 TFile *f = TFile::Open(treefilename);
 TFile *fout;
 TTree *tout;
 f->cd();
 string runstr(f->GetName()); runstr=runstr.substr(runstr.find_last_of("/")+1,10);
 char crunstr[2048]; sprintf(crunstr,"%s",runstr.c_str());
 int Npoints = 0;
 Int_t dt_shift_int = (Int_t)(1e12*dt_shift);
 Double_t dt_sampling = 0;
 cout<<"dt_shift_int "<<dt_shift_int<<"ps\n";
 
 if (debug) printf("line: %d\n",__LINE__);
 
 //#### CHECK Ntriggers from f
 Int_t Ntriggers = 0;
 
 for( int ig1=0; ig1<10000000; ig1++){
 TGraph *g1 = (TGraph*)f->Get( Form("gwf_%08d",ig1) );
 if( !g1 ) break;
 Ntriggers++;
 delete g1;
 }
 
 cout<<"Ntriggers -----> "<< Ntriggers  <<"\n";
 if (debug) printf("line: %d\n",__LINE__);
 printf("%s::Set-Triggers-to-%d\n", __FUNCTION__, Ntriggers);
 
 if (debug) printf("line: %d\n",__LINE__);
 
 //################# allocation data* : TGraphs* , Peaks*    evita di portarti dietro vgwf vgwfnopeak vgwfdlednopeak e in generale
 
 vector<TGraph*> vgwf;                    //original waveform
 vector<TGraph*> vgwfsm;                  //shooth waveform
 vector<TGraph*> vgdled;                  //shooth waveform tranformed with DLED
 vector<TGraph*> vgwfnopeak;              //original waveform (without peaks)
 vector<TGraph*> vgwfdlednopeak;          //vgdled (without peaks)
 vector<Peaks*>  vpeaks;                   // original peaks stored
 vector<Peaks*>  vpeaksdled;               // DLED peaks stored
 vgwf.resize(Ntriggers);        for(int it=0; it<Ntriggers; it++) vgwf.at(it) = NULL;
 vgwfsm.resize(Ntriggers);      for(int it=0; it<Ntriggers; it++) vgwfsm.at(it) = NULL;
 vgdled.resize(Ntriggers);      for(int it=0; it<Ntriggers; it++) vgdled.at(it) = NULL;
 vgwfnopeak.resize(Ntriggers);  for(int it=0; it<Ntriggers; it++) vgwfnopeak.at(it) = NULL;
 vgwfdlednopeak.resize(Ntriggers);  for(int it=0; it<Ntriggers; it++) vgwfdlednopeak.at(it) = NULL;
 vpeaks.resize(Ntriggers);      for(int it=0; it<Ntriggers; it++) vpeaks.at(it) = new Peaks();
 vpeaksdled.resize(Ntriggers);  for(int it=0; it<Ntriggers; it++) vpeaksdled.at(it) = new Peaks();
 if (debug) printf("line: %d\n",__LINE__);
 //################# Clone data from file to storage (vectors)
 
 //for( int ig=0; ig<Ntriggers; ig++){
 //TGraph *g = (TGraph*)f->Get( Form("gwf_%08d",ig) );
 //if( !g ) { printf("Cannot find gwf_%08d\n Exit....",ig); exit(78); }
 //vgwf.at(0) = g;//(TGraph*)(g->Clone(Form("gwf_%08d",ig)));
 //printf("%s\n",vgwf.at(0)->GetName());
 //}
 
 Int_t ntbins=200;
 Double_t *Naxis = GenerateBinning(1000,-10,10);
 //// taxis è l'asse dei dt
 Double_t *hz_axis = GenerateBinning(100,0,0.5);
 Double_t *taxis = GenerateBinning(ntbins,-1e+2,1e+3);//-1e+2,5e+3
 Double_t *taxis_log = GenerateLogBinning(ntbins,1e-1,1e+3);//-1e+2,5e+3
 Int_t qbins = 1000;
 Double_t *binaxis = GenerateBinning(1000,-400,1200);
 Double_t *qaxis = GenerateBinning(qbins,-20,200);
 Double_t *ttrgaxis = GenerateBinning(ntbins,90,105);
 Int_t nvbins=1000;
 Double_t *vaxis = GenerateBinning(nvbins,vaxismin,vaxismax);
 Double_t *vaxis_dv = GenerateBinning(1000,-vaxisdledmax,vaxisdledmax);
 Double_t *vaxisdlednoth = GenerateBinning(nvbins,vaxisdledmin,vaxisdledmax);
 Double_t *vaxisdled = GenerateBinning(nvbins,vaxisdledmin,vaxisdledmax);// peakthreshold al posto di vaxisdledmin
 Double_t *vaxislog = GenerateLogBinning(nvbins,0.1,100);
 Double_t *vaxislogdled = GenerateLogBinning(nvbins,0.1,50);
 Int_t ntimebins=200;//20
 Double_t *timeaxis = GenerateBinning(ntimebins,mintime,maxtime);
 Double_t *eventaxis = GenerateBinning(ntbins,0,20); // (int)t->GetEntries() al posto di 20
 Float_t tree_Vmax, tree_Tmax, tree_dt,tree_Vmax_dled, tree_Tmax_dled, tree_dt_dled,tree_Q, tree_deltaV, tree_deltaVdled, tree_deltaQ,tree_Ceff,tree_ped, tree_ped_dled,tree_deltaVdled_norm, tree_deltaVdled_2,tree_Vmax_dled_norm,tree_Vmax_dled_prev, tree_T_after;
 Bool_t  tree_goodQ;
 UShort_t tree_wf;
 
 
 //################# creating output file & tree
 if (debug) printf("line: %d\n",__LINE__);
 status=mkdir(Form("output_%s_%dps",crunstr,dt_shift_int), 0700);//dt_shift_int
 cout<<"directory creation[ "<< Form("output_%s_%dps",crunstr,dt_shift_int)  <<" ](good if 0) --> "<<status<<endl;
 fout = new TFile( Form("output_%s_%dps/fout_drawnolaser_%s_%dps.root",crunstr,dt_shift_int,crunstr,dt_shift_int),"recreate");
 fout->cd();
 tout = new TTree("tree","tree");
 if (debug) printf("line: %d\n",__LINE__);
 
 tout->Branch("Vmax",      &tree_Vmax,      "Vmax/F");
 tout->Branch("Tmax",      &tree_Tmax,      "Tmax/F");
 tout->Branch("dt",        &tree_dt,        "dt/F");
 tout->Branch("wf",        &tree_wf,        "wf/s");
 tout->Branch("Vmax_dled", &tree_Vmax_dled, "Vmax_dled/F");
 tout->Branch("Vmax_dled_norm", &tree_Vmax_dled_norm, "Vmax_dled_norm/F");
 tout->Branch("Vmax_dled_prev", &tree_Vmax_dled_prev, "Vmax_dled_prev/F");
 tout->Branch("Tmax_dled", &tree_Tmax_dled, "Tmax_dled/F");
 tout->Branch("dt_dled",   &tree_dt_dled,   "dt_dled/F");
 tout->Branch("Q",         &tree_Q,         "Q/F");
 tout->Branch("deltaQ",    &tree_deltaQ,         "deltaQ/F");
 tout->Branch("Ceff",      &tree_Ceff,         "Ceff/F");
 tout->Branch("goodQ",     &tree_goodQ,     "goodQ/O");
 tout->Branch("deltaV",    &tree_deltaV,      "deltaV/F");
 tout->Branch("deltaVdled",&tree_deltaVdled,      "deltaVdled/F");
 tout->Branch("deltaVdled_norm",&tree_deltaVdled_norm,      "deltaVdled_norm/F");
 tout->Branch("deltaVdled_2",&tree_deltaVdled_2,      "deltaVdled_2/F");
 tout->Branch("ped",       &tree_ped,      "tree_ped/F");
 tout->Branch("ped_dled",       &tree_ped_dled,      "tree_ped_dled/F");
 tout->Branch("dt_after",   &tree_T_after,   "dt_after/F");
 
 
 
 TH1F *dt = new TH1F("dt","dt",ntbins,taxis_log);
 TH1F *hgwfnopeaks = new TH1F("hgwfnopeaks","hgwfnopeaks",nvbins,vaxisdled);//300,-10,4);
 TH2F *prova_t_dt = new TH2F("prova_t_dt","prova tempi t_dt;t_v;dt_v",ntimebins,timeaxis,ntbins,taxis);
 TH2F *prova_t_dt_dled = new TH2F("prova_t_dt_dled","prova tempi t_dt_dled;t_dled;dt_dled",ntimebins,timeaxis,ntbins,taxis);
 TH2F *prova_tempi = new TH2F("prova_tempi","prova tempi;dt_v;dt_dled",ntbins,taxis,ntbins,taxis);
 TH2F *prova_tempi_t = new TH2F("prova_tempi_t","prova tempi_t;t_v;t_dled",ntbins,taxis,ntbins,taxis);
 TH2F *prova_V = new TH2F("prova_V","prova V;v;vdled",nvbins,vaxis,nvbins,vaxisdled);
 TH2F *prova_N_dt_norm = new TH2F("prova_N_dt_norm","prova_N_dt_norm;dt;N",ntbins,taxis,100,Naxis);
 TH2F *prova_N_hz_norm = new TH2F("prova_N_hz_norm","prova_N_hz_norm;GHz;N",100,hz_axis,100,Naxis);
 TH2F *prova_N_dv_norm = new TH2F("prova_N_dv_norm","prova_N_dv_norm;dled ;N",nvbins,vaxisdled,1000,Naxis);//nvbins,vaxisdled,100,Naxis
 TH2F *prova_N_dv_2 = new TH2F("prova_N_dv_2","prova_N_dv_2;dled mV;N",nvbins,vaxisdled,1000,0,6);//,nvbins,vaxisdled
 TH2F *hVmaxdt = new TH2F("hVmaxdt","Scope waveform;#Delta T from previous peak (ns);Amplitude (mV)",ntbins,taxis_log,nvbins,vaxis);
 TH2F *hVmaxdt_prova = new TH2F("hVmaxdt_prova","Scope waveform;#Delta T from previous peak (ns);Amplitude (mV)",ntbins,0.1,2000,1000,-4,20);
 TH2F *hVmaxdtdled = new TH2F("hVmaxdtdled","DLED waveform;#Delta T from previous peak (ns);Amplitude (mV)",ntbins,taxis_log,nvbins,vaxisdled);
 TH2F *hVmaxdtdled_prova= new TH2F("hVmaxdtdled_prova","DLED waveform;#Delta T from previous peak (ns) time ;Amplitude (mV)",ntbins,2000,4500,1000,0,2);
 TH1F *hVnosig = new TH1F("hVnosig","V samples with no signal;Amplitude (mV)", nvbins,vaxisdled);//200,-10, 10);
 TH1F *hVnosigdled = new TH1F("hVnosigdled","V samples with no signal;DLED Amplitude (mV)",nvbins,-vaxisdledmax/5,vaxisdledmax/5);
 TH1F *hVmaxdled = new TH1F("hVmaxdled","DLED waveform;V peak (mV);Entries",nvbins,vaxisdled);
 TH1F *hVmax = new TH1F("hVmax","hvmax;V peak (mV);Entries",nvbins,vaxis);
 TH1F *hVmaxdledlog = new TH1F("hVmaxdledlog","hvmaxdledlog;V peak (mV);Entries",nvbins,vaxislogdled);
 TH1F *hVmaxlog = new TH1F("hVmaxlog","hvmaxlog;V peak (mV);Entries",nvbins,vaxislog);
 TH1F *hdtdled = new TH1F("hdtdled","hdtdled;#Delta T from previous peak (ns);Entries",ntbins,taxis_log);
 TH1F *hdt = new TH1F("hdt","hvmax;#Delta T from previous peak (ns);Entries",ntbins,taxis_log);
 TH1F *hVprepeak = new TH1F("hVprepeak","hVprepeak;Amplitude (mV);Entries",300,-vaxismax/10,vaxismax/10);
 TH2F *hVnosigtime = new TH2F("hVnosigtime","V samples with no signal;Time;Amplitude (mV)",ntimebins,timeaxis,nvbins,vaxis);
 TH2F *hVmaxtime = new TH2F("hVmaxtime","hVmaxtime;Time;V peak (mV)",ntimebins,timeaxis,nvbins,vaxis);
 TH2F *hVmaxdledtime = new TH2F("hVmaxdledtime","hVmaxdledtime;Time;V peak (mV)",ntimebins,timeaxis,nvbins,vaxisdlednoth);
 TH2F *hVmaxVmaxdled = new TH2F("hVmaxVmaxdled","hVmaxVmaxdled;Vmax;Vmax dled",nvbins, vaxis, nvbins, vaxisdled);
 TH2F *hVmaxtrgtime = new TH2F("hVmaxtrgtime","hVmaxtrgtime;Time after trigger;V peak (mV)",ntbins,ttrgaxis,nvbins,vaxis);
 TH2F *hVmaxdledtrgtime = new TH2F("hVmaxdledtrgtime","hVmaxdledtrgtime;Time after trigger;V peak (mV)",ntbins,ttrgaxis,nvbins,vaxisdlednoth);
 TH2F *hVmaxdledentry = new TH2F("hVmaxdledentry","hVmaxdledentry;Tree Entry;V peak (mV)",ntbins,eventaxis,nvbins,vaxisdlednoth);
 //TH2F *ritardo_zero = new TH2F("ritardo_zero","ritardo_zero;t_v;dt_v",1000,0,10,100,-10,10);
 //TH2F *ritardi_piccoli = new TH2F ("ritardi_piccoli","ritardi_piccoli;t_v;dt_v",1000,0,10,100,-10,10);
 Double_t ped_dled=0;
 Double_t ped_signal=0;
 TH1F *hdled = new TH1F("hdled","hdled",1000,vaxisdledmin,vaxisdledmax);
 TH1F *hsign = new TH1F("hsign","hsign",1000,-20,vaxismax);///hsign
 TGraphErrors *gtdc = new TGraph(); gtdc->SetNameTitle("gtdc","Dark Count Rate;Amplitude (mV);Rate (kHz/mm^{2})");
 TH2F *Q_dt= new TH2F("carica_dt","carica vs_dt;#Delta T from previous peak (ns);# electrons (x 10^{6})",ntbins,taxis,qbins,qaxis);
 TH1F *Q = new TH1F("carica","carica;#electrons (x 10^{6})",qbins,qaxis);
 TH1F *hCeff = new TH1F("Ceff","Ceff; Q/V (unità da definire)",1000,-1,1);
 TH2F *Q_dled= new TH2F("carica_dled","carica vs dled;# electrons (x 10^{6});V peaks (mV)",qbins,qaxis,nvbins,vaxisdled);
 TH2F *Q_Vmax= new TH2F("carica_Vmax","carica vs Vmax;V peaks (mV);# electrons (x 10^{6})",nvbins,vaxis,qbins,qaxis);
 TH2F *Q_histo2d = new TH2F("Q_histo2d","Q_histo2d;bin;volt",ntbins,taxis,1000,vaxis_dv);
 TH1F *delta_V = new TH1F("delta_V","delta_v;mV",1000,-10,10);//-vaxismax,vaxismax);
 TH1F *delta_V_dled = new TH1F("delta_V_dled","delta_vdled;mV",1000,-10,10);//-vaxisdledmax,vaxisdledmax);
 TH1F *delta_q = new TH1F("delta_q","delta_q;unità da definire",1000,-50,50);
 TH1F *hped = new TH1F("ped","ped;mV",1000,-20,40);
 TH1F *hped_dled = new TH1F("ped_dled","ped_dled;mV",1000,-20,40);
 //TH1F *hdled = new TH1F("hdled","hdled",1000,-0.005, 0.008);
 TAxis *xaxis_hsign = hsign->GetXaxis();
 
 /////////////////MEGA__LOOP////////////////
 if (debug) printf("line: %d\n",__LINE__);
 
 
 //################# checks for loop parameters
 if(imin<0) imin=0;
 if(imax<0 || imax > Ntriggers ) imax=Ntriggers;
 Double_t totaltime = 0;
 Int_t Nanaltriggers=0;
 float tot_fract = imax-imin+1;
 int cont = 1;
 if (debug) printf("line: %d\n",__LINE__);
 //################# The loop begin
 for(int i=imin; i<imax;i++){
 f->cd();
 TGraph *g = (TGraph*)f->Get( Form("gwf_%08d",i) );
 if( !g ) break;
 vgwf.at(0) = g;
 dt_sampling = (vgwf.at(0)->GetX()[1001])-(vgwf.at(0)->GetX()[1000]);
 cout<< "dt_sampling "<<dt_sampling*1e9<<" ns ------ dt_sampling "<<dt_sampling*1e12<<" ps\n";
 if( i==0 ) {
 Npoints=(int)g->GetN();
 printf("%s::Set-Waveform-Points-to-%d\n", __FUNCTION__, Npoints);
 }
 
 fout->cd();
 
 
 if( (i/tot_fract) >= ((1./20)*cont)){printf("%.2d%%\n",cont*100/20); cont++;}
 if (debug) printf("line: %d\n",__LINE__);
 TGraph* gwf   = vgwf.at(0); ;
 TGraph* gwfsm = vgwfsm.at(0);
 TGraph* gdled = vgdled.at(0);
 TGraph* gwfdlednopeak = vgwfdlednopeak.at(0);
 TGraph* gwfnopeak = vgwfnopeak.at(0);
 if (gwf == NULL)                     printf("gwf NO GRAFICO_line: %d\n",__LINE__);
 if (gwfsm != NULL )                  printf("gwfsm  GRAFICO_line: %d\n",__LINE__);
 if (gdled != NULL)                   printf("gdled  GRAFICO_line: %d\n",__LINE__);
 if (gwfdlednopeak != NULL)           printf("gwfdlednopeak NO GRAFICO_line: %d\n",__LINE__);
 if (vgwfdlednopeak.at(i) != NULL )   printf("gwfdlednopeak NO elemento_line: %d\n",__LINE__);
 if (gwfnopeak != NULL)               printf(" gwfnopeak NO GRAFICO_line: %d\n",__LINE__);
 if (vgwfnopeak.at(i) != NULL )       printf("gwfnopeak NO elemento_line: %d\n",__LINE__);
 Peaks *peaks = vpeaks.at(i);
 Peaks *peaksdled = vpeaksdled.at(i);
 if (peaks == NULL)                   printf("peaks problem_line: %d\n",__LINE__);
 if (vpeaks.at(i) == NULL )           printf("peaks NO elemnt in vector_line: %d\n",__LINE__);
 if (peaksdled == NULL)               printf("peaksdled problem_line: %d\n",__LINE__);
 if (vpeaksdled.at(i) == NULL )       printf("vpeaksdled NO element in vector_line: %d\n",__LINE__);
 if (debug) printf("line: %d\n",__LINE__);
 
 if(debug) cout<<"prova ----------  gwf->GetN()  "<< gwf->GetN() << endl;
 for (int i_for =0;i_for<gwf->GetN();i_for++) gwf->GetY()[i_for] *= -1 ;
 TGraphRunningMean( gwf, &(gwfsm), run_mean_points);
 cout<<"running mean "<<run_mean_points<<" points == "<<run_mean_points*dt_sampling*1e9<<"ns \n";
 gwfsm->SetName( Form("%s_sm",gwf->GetName() ) );
 if (i==imin) {gwf->Write();gwfsm->Write();}/// saved for visual check
 if (debug) printf("line: %d\n",__LINE__);
 if(debug) cout<<"gwfsm "<<gwfsm << " gdled "<<gdled<<endl;
 if (gwfsm == NULL ){printf("NULL gwfsm line: %d\n",__LINE__);}
 else {if(debug) {printf("good ----- not NULL gwfsm line: %d\n",__LINE__);}}
 if(debug) cout<<"gwfsm_X "<<gwfsm->GetX() <<endl;
 
 for(int ip=1; ip<gwfsm->GetN(); ip++){
 hsign->Fill( (gwfsm->GetY()[ip])*1e3 );
 }
 
 gDLEDpoint( gwfsm, &(gdled), (-1)*point_shift_dled); //-1 point shift 8e-10
 cout<<"dled "<<point_shift_dled<<" points == "<<point_shift_dled*dt_sampling*1e9<<"ns \n";
 if (debug) printf("line: %d\n",__LINE__);
 for(int ip=1; ip<gdled->GetN(); ip++){
 hdled->Fill( (gdled->GetY()[ip])*1e3 );
 }
 if (debug) printf("line: %d\n",__LINE__);
 //ScaleGraph(gdled,1);
 gdled->SetName( Form("g_%d_dled",i) );
 if (i==imin) {gwf->Write();gwfsm->Write();gdled->Write();}/// saved for visual check
 if(debug) gdled->Write();
 if (debug) printf("line: %d\n",__LINE__);
 ScaleGraph(gdled,1);
 if (debug) {fout->cd(); gdled->Write();}
 
 totaltime += (gwf->GetX()[gwf->GetN()-1] - gwf->GetX()[0]) * ( (gwf->GetN()+1)/(gwf->GetN()));
 //total time = (last X - first X) plus "half bin" before first X and after last Xs
 Nanaltriggers++;
 FindPeaks( gdled, peaksdled, peakthreshold,false, false );
 if (debug) printf("line: %d\n",__LINE__);
 FindNextPeaks( gwfsm, peaksdled, peaks, 5e-9, false );
 if (debug) printf("line: %d\n",__LINE__);
 RemovePeaks( gwfsm , &gwfnopeak, peaksdled);
 if (debug) printf("line: %d\n",__LINE__);
 RemovePeaks( gdled, &gwfdlednopeak, peaksdled);
 ped_signal = gwfnopeak->GetMean(2);
 cout<<"ped_signal "<<ped_signal*1e3<<" mV\n";
 ped_signal = (xaxis_hsign->GetBinCenter(hsign->GetMaximumBin()))*1e-3;
 cout<<"ped_signal_hsign "<<ped_signal*1e3<<" mV\n";
 integra_picco( gwfsm, peaksdled, peaks, ped_signal ,10 , 300, dt_sampling*1e9); //gwf, 10 , 300  //all in ns , voltages in mV
 
 
 if (i==imin) {gwfnopeak->Write(); gwfdlednopeak->Write(); cout<<"1-up\n";}
 for(int ip=0; ip<(int)gwfsm->GetN(); ip++){
 if(debug && (i==1)){gwfsm->Write();}
 hVprepeak->Fill((gwfsm->GetY()[ip]-ped_signal)*1e3);
 }
 for(int ip=0; ip<gwfnopeak->GetN(); ip++){
 hVnosig->Fill((gwfnopeak->GetY()[ip]-ped_signal)*1e3);
 //hVnosigtime->Fill(utime,gwfnopeak->GetY()[ip]-ped_signal);////// utime it's pratically the time stamp -> i don't use it so hVnosigtime now has the real peak time
 hVnosigtime->Fill((gwfnopeak->GetX()[ip]),(gwfnopeak->GetY()[ip]-ped_signal)*1e3);
 
 }
 for(int ip=0; ip<gwfdlednopeak->GetN(); ip++){
 hVnosigdled->Fill((gwfdlednopeak->GetY()[ip]-ped_dled)*1e3);
 }
 for(int ip=0;  ip<peaksdled->N; ip++) {
 if(peaksdled->dt.at(ip)!=-999){////////////////////////////qui si aggiunge filtro dt>ns es. 400ns
 hVmaxdled->Fill((peaksdled->V.at(ip)-ped_dled)*1e3);
 }
 if(debug) cout<<"check1\n";
 tree_ped = (Float_t)ped_signal*1e3;
 tree_ped_dled =(Float_t) ped_dled*1e3;
 tree_Vmax_dled = (Float_t)((peaksdled->V.at(ip)-ped_dled)*1e3);//mV
 tree_Tmax_dled = (Float_t)1e+9*peaksdled->t.at(ip);//ns
 tree_dt_dled = (Float_t)1e+9*peaksdled->dt.at(ip);//ns
 if(ip < (peaksdled->N)-1 ){tree_T_after = (Float_t)1e+9*(peaksdled->t.at(ip+1)-peaksdled->t.at(ip));}
 else{tree_T_after = 0;}
 tree_Vmax = (Float_t)(peaks->V.at(ip)-ped_signal)*1e3;//mV
 tree_Tmax = (Float_t)1e+9*peaks->t.at(ip);//ns
 tree_dt = (Float_t)1e+9*peaks->dt.at(ip);//ns
 tree_wf = (UShort_t)i;
 tree_Q = (Float_t)1e3*1e-6*peaksdled->Q.at(ip);//milioni di elettroni (1e3 è per usare le wf in V->mV)
 if(ip>=1){
 tree_Vmax_dled_prev = (Float_t)((peaksdled->V.at(ip-1))*1e3-tree_ped_dled);//mV
 tree_Vmax_dled_norm = (tree_Vmax_dled*tree_Vmax_dled_prev)/pow(tree_Vmax_dled + (peaksdled->V.at(ip-1)-ped_dled)*1e3,2);
 tree_deltaV =  (Float_t) (tree_Vmax-(peaks->V.at(ip-1)-ped_signal)*1e3); ///(Float_t) tree_Vmax-(peaks->V.at(ip-1)-ped_signal)*1e3;
 tree_deltaVdled = (Float_t)(tree_Vmax_dled - (peaksdled->V.at(ip-1)-ped_dled)*1e3); /// tolto TMath::Abs
 tree_deltaVdled_norm = (Float_t) TMath::Abs((tree_Vmax_dled - (peaksdled->V.at(ip-1)-ped_dled)*1e3))/(tree_Vmax_dled + (peaksdled->V.at(ip-1)-ped_dled)*1e3);     //(Float_t) tree_Vmax_dled - (peaksdled->V.at(ip-1)-ped_dled)*1e3;
 tree_deltaVdled_2 = (Float_t) pow((tree_Vmax_dled - (peaksdled->V.at(ip-1)-ped_dled)*1e3),2)/(tree_Vmax_dled + (peaksdled->V.at(ip-1)-ped_dled)*1e3);
 tree_deltaQ = (Float_t)tree_Q-1e3*1e-6*peaksdled->Q.at(ip-1);
 tree_Ceff = (Float_t)(1.6e-19*1e6*1e12*tree_Q/tree_Vmax);/// dovrebbe essere in pico-farad ma le tensioni in mV -> nF ?
 }
 tree_goodQ = GoodCharge( peaksdled, ip );
 tout->Fill();
 }
 if(debug) cout<<"check2\n";
 if(gwf)       delete gwf;
 if(gwfsm)     delete gwfsm;
 if(gdled)     delete gdled;
 if(gwfnopeak) delete gwfnopeak;
 if(peaks)     delete peaks;
 if(peaksdled) delete peaksdled;
 if(gwfdlednopeak) delete gwfdlednopeak;
 //if(g) delete g;
 // cout<<"\n deleted g ->"<< i <<"\n";
 
 }
 
 cout<<"end loop -> filling histograms from tree \n";
 
 hsign->Write();
 ////////MANIPOLA TREE/////////////
 for(int entry_tree = 0; entry_tree < tout->GetEntries();entry_tree++ ){
 tout->GetEntry(entry_tree);
 if(tree_goodQ==1){ ///0.28
 //    if(tree_dt < 300) continue; filtro ritardi
 
 delta_V->Fill(tree_deltaV);
 //dummy_delta_v_dled = tree_Vmax_dled - dummy_voltage_dled;
 //if((tree_deltaVdled > 0.26) && (tree_deltaVdled < 0.36)){delta_V_dled->Fill(tree_Vmax_dled);}
 
 delta_V_dled->Fill(tree_deltaVdled);
 prova_N_dt_norm->Fill(tree_dt_dled,tree_deltaVdled_norm); ///tree_Vmax_dled
 prova_N_dv_norm->Fill(tree_Vmax_dled,(2/(1/tree_deltaVdled_norm)-1));
 prova_N_hz_norm->Fill(tree_dt_dled,tree_deltaVdled_norm);
 prova_N_dv_2->Fill(tree_Vmax_dled,tree_deltaVdled_2);//tree_deltaVdled_2
 //delta_q->Fill(-(Float_t)TMath::ATan((Double_t)(tree_Vmax/tree_Q)));
 ///se è uniforme -> prob distr di cauchy di -acot(Q/V)
 delta_q->Fill(tree_deltaQ);
 //delta_q->Fill(-(Float_t)TMath::ATan((Double_t)(tree_Q/dummy_q)));///1/C = tree_Vmax/tree_Q -> Cauchy
 
 Q_histo2d->Fill(tree_dt_dled,tree_deltaVdled);
 hCeff->Fill(tree_Ceff);
 hped->Fill(tree_ped);
 hped_dled->Fill(tree_ped_dled);
 prova_t_dt->Fill(tree_Tmax,tree_dt);
 prova_t_dt_dled->Fill(tree_Tmax_dled,tree_dt_dled);
 prova_tempi->Fill(tree_dt,tree_dt_dled);
 prova_tempi_t->Fill(tree_Tmax,tree_Tmax_dled);
 if(((tree_deltaVdled > 0.4) && (tree_deltaVdled < 0.6))){ prova_V->Fill(tree_Vmax,tree_Vmax_dled);}
 dt->Fill(tree_dt);
 hVmaxdt->Fill(tree_dt, tree_Vmax);
 if(tree_dt > 400){hVmax->Fill(tree_Vmax);}
 hVmaxlog->Fill(tree_Vmax);
 if(1){hdt->Fill(tree_dt);} // tree_Tmax > 2000
 hVmaxtime->Fill(tree_Tmax_dled,tree_Tmax);// utime ->tree_Tmax_dled, utime è il timestamp ->non lo uso, uso il tempo del picco
 hVmaxtrgtime->Fill(tree_Tmax,tree_Vmax);
 hVmaxdtdled->Fill(tree_dt_dled,tree_Vmax_dled );
 hVmaxdledlog->Fill(tree_Vmax_dled);
 if(1){ hdtdled->Fill(tree_dt_dled);} /// tree_Tmax_dled > 2000
 hVmaxdledtime->Fill(tree_dt_dled,tree_Vmax_dled);// ltime
 hVmaxdledtrgtime->Fill(tree_Tmax_dled,tree_Vmax_dled);
 hVmaxdledentry->Fill(entry_tree,tree_Vmax_dled);
 hVmaxVmaxdled->Fill(tree_Vmax,tree_Vmax_dled);
 
 //Double_t carica_temp_i, V_temp_i, ritardo_i, Vmax_i;
 if(entry_tree > 0 && entry_tree+1 < tout->GetEntries()){
 //carica_temp_i = tree_Q;
 //ritardo_i = tree_dt_dled;
 //V_temp_i = tree_Vmax_dled;
 //Vmax_i = tree_Vmax;
 //tout->GetEntry(entry_tree+1);
 //if(ritardo_i>400 && tree_dt_dled> 400){
 Q->Fill(tree_Q);
 //cout<<"tree_Q example "<<tree_Q<<"\n";
 Q_dt->Fill(tree_dt_dled,tree_Q);
 Q_Vmax->Fill(tree_Vmax,tree_Q);
 Q_dled->Fill(tree_Q,tree_Vmax_dled);
 //}
 
 //tout->GetEntry(entry_tree);
 }
 }
 
 }
 //da rivedere , manca la def di Q_histo2d
 //TGraph *picco_graph = new TGraph();
 //picco_medio(Q_histo2d,picco_graph);
 //picco_graph->SetNameTitle("picco_graph"," picco_graph");
 //picco_graph->Write();
 //
 cout<<"---------------DCR---------------\n";
 printf("%s::Total-Triggers-%d\n",__FUNCTION__,Nanaltriggers);
 printf("%s::Total-Time-%.3e\n",__FUNCTION__,totaltime);
 printf("%s::Average-Time-Per-TriggerT-%.3e\n",__FUNCTION__,totaltime/Nanaltriggers);
 TH1F *h = hVmaxdled;
 for(int ibin=1; ibin<=(h->GetNbinsX()); ibin++){/// aggiungi +1 se sommi overflow
 
 Double_t peaks_above_this_threshold = h->Integral(ibin,h->GetNbinsX()+1);
 Double_t rate_of_dark_counts_above_this_threshold = peaks_above_this_threshold / totaltime;
 gtdc->SetPoint(gtdc->GetN(),h->GetBinCenter(ibin),(rate_of_dark_counts_above_this_threshold/1e+3)/sipmarea);
 cout<<"totaltime "<<totaltime<<"total triggers "<<Nanaltriggers<<"-------------- controlla che sia in ms\n";
 gtdc->SetPointError(gtdc->GetN(),0,TMath:Sqrt((rate_of_dark_counts_above_this_threshold/1e+3)/sipmarea*totaltime));
 }
 TCanvas *cVmaxdt = new TCanvas("cVmaxdt","cVmaxdt"); arrCanvas->Add(cVmaxdt);
 cVmaxdt->Divide(2,2);
 cVmaxdt->cd(1)->SetGrid(); cVmaxdt->cd(1)->SetLogx(); cVmaxdt->cd(1)->SetLogz();
 hVmaxdt->Draw("COLZ");
 cVmaxdt->cd(2);
 hVmax->Draw("");
 cVmaxdt->cd(3)->SetLogz();  cVmaxdt->cd(3)->SetGrid();
 TH2F *hVmaxtime_normx = NormalizeX(hVmaxtime); SetHTime(hVmaxtime_normx); hVmaxtime_normx->Draw("COLZ");
 cVmaxdt->cd(4)->SetLogx();
 hdt->Draw("");
 
 TCanvas *cVmaxdtdled = new TCanvas("cVmaxdtdled","cVmaxdtdled"); arrCanvas->Add(cVmaxdtdled);
 cVmaxdtdled->Divide(2,2);
 cVmaxdtdled->cd(1)->SetGrid(); cVmaxdtdled->cd(1)->SetLogx(); cVmaxdtdled->cd(1)->SetLogz();
 hVmaxdtdled->Draw("COLZ");
 cVmaxdtdled->cd(2);
 hVmaxdled->Draw("");
 cVmaxdtdled->cd(3)->SetLogz();  cVmaxdtdled->cd(3)->SetGrid();
 TH2F *hVmaxdledtime_normx = NormalizeX(hVmaxdledtime); SetHTime(hVmaxdledtime_normx); hVmaxdledtime_normx->Draw("COLZ");
 cVmaxdtdled->cd(4)->SetGrid(); cVmaxdtdled->cd(4)->SetLogz();
 hVmaxdledentry->Draw("COLZ");
 TCanvas *cVmaxVmaxdled = new TCanvas("cVmaxVmaxdled","cVmaxVmaxdled"); arrCanvas->Add(cVmaxVmaxdled);
 cVmaxVmaxdled->Divide(2,2);
 cVmaxVmaxdled->cd(1)->SetGrid();
 hVmaxVmaxdled->Draw("COLZ");
 cVmaxVmaxdled->cd(3)->SetGrid();
 hVmaxtrgtime->Draw("COLZ");
 cVmaxVmaxdled->cd(4)->SetGrid();
 hVmaxdledtrgtime->Draw("COLZ");
 
 TCanvas *cVnosig = new TCanvas("cVnosig","cVnosig");  arrCanvas->Add(cVnosig);
 cVnosig->Divide(2,2);
 cVnosig->cd(1)->SetGrid(); cVnosig->cd(1)->SetLogz();
 TH2F *hVnosigtime_normx = NormalizeX(hVnosigtime); SetHTime(hVnosigtime_normx); hVnosigtime_normx->Draw("COLZ");
 cVnosig->cd(2)->SetLogy();
 TH1F *hVnosigtime_py = (TH1F*)hVnosigtime->ProjectionY(); hVnosigtime_py->Draw("");
 TGraph *gVnosigtime_mean = new TGraph();
 TGraph *gVnosigtime_rms = new TGraph();
 for(int ibin=1; ibin<=hVnosigtime->GetNbinsX(); ibin++){
 TH1F *h = (TH1F*)hVnosigtime->ProjectionY(Form("hVnosigtime_py_%d",ibin),ibin,ibin);
 gVnosigtime_mean->SetPoint(gVnosigtime_mean->GetN(),hVnosigtime->GetXaxis()->GetBinCenter(ibin),h->GetMean());
 gVnosigtime_rms ->SetPoint(gVnosigtime_rms ->GetN(),hVnosigtime->GetXaxis()->GetBinCenter(ibin),h->GetRMS());
 delete h;
 }
 cVnosig->cd(3)->SetGrid(); gVnosigtime_mean->SetNameTitle("gVnosigtime_mean",";Time;Pedestal"); gVnosigtime_mean->Draw("AP");
 cVnosig->cd(4)->SetGrid(); gVnosigtime_rms ->SetNameTitle("gVnosigtime_rms", ";Time;Ped RMS");  gVnosigtime_rms->Draw("AP");
 
 
 
 TCanvas *cgtdc = new TCanvas("cgtdc","cgtdc"); cgtdc->cd()->SetGrid();  cgtdc->cd()->SetLogy(); arrCanvas->Add(cgtdc);
 gtdc->SetLineWidth(3); gtdc->Draw("AL");
 gtdc->Write();//<--------------------------------------------------------------
 Double_t gtdvmax = TMath::MaxElement( gtdc->GetN(), gtdc->GetY() );
 TH1 * hVmaxdledscaled = (TH1*)hVmaxdled->Clone("hVmaxdledscaled"); hVmaxdledscaled->Scale( gtdvmax / hVmaxdled->Integral() );
 hVmaxdledscaled->Draw("same");
 if (debug) printf("line: %d\n",__LINE__);
 
 /// FINE programma (salvataggio e poco altro)
 if(print){
 status=mkdir(Form("output_%s_%dps/waveforms_%s",crunstr,dt_shift_int,crunstr),0700);
 cout<<"directory creation[ "<< Form("output_%s/waveforms_%s",crunstr,crunstr) <<" ](good if 0) --> "<<status<<endl;
 PrintCanvas( arrwfCanvas, Form("output_%s_%dps/waveforms_%s",crunstr,dt_shift_int,crunstr) , "drawnolaser.pdf");
 }
 
 status=mkdir(Form("output_%s_%dps/plots_%s",crunstr,dt_shift_int,crunstr),0700);
 cout<<"directory creation[ "<< Form("output_%s/plots_%s",crunstr,crunstr) <<" ](good if 0) --> "<<status<<endl;
 PrintCanvas( arrCanvas, Form("output_%s_%dps/plots_%s",crunstr,dt_shift_int,crunstr), "drawnolaser.pdf" );
 cout<<"salvataggio file\n";
 
 fout->Write();
 gtdc->Write();
 tout->Write();
 fout->Close();
 
 printf("fine scrittura_line: %d\n",__LINE__);
 
 gROOT->Reset();
 
 printf("COMPLETED \n");
 return;
 
 }
 
 */
