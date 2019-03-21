#include "general_lib.h"
#include "Peaks.h"
#include "functions_sipm_all.h"


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


