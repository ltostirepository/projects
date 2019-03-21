#include "general_lib.h"
#include "Peaks.h"
#include "functions_sipm_smoothing_and_DLED.h"

using namespace std;


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
