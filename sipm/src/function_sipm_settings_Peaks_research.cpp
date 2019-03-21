#include "general_lib.h"
#include "Peaks.h"
#include "functions_sipm_Peaks_research.h"

using namespace std;



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
            //Double_t maxside=Y[xmaxside];
            //Double_t ave = 0.5*(minside+maxside);
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
    //Double_t dx = (y[n-1]-y[0])/(n-1);
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



