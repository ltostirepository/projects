#include "general_lib.h"
#include "Peaks.h"
#include "functions_sipm_print_and_trials.h"

using namespace std;


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


