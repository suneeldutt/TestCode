#include <iostream>
#include "TH1.h"
#include <string>

using namespace std;

//TH1D *CreateHisto(string name, int nbins, double nmin, double nmax, string title, string xtitle, string ytitle, bool CenterTitle, bool rebin);
TH1D *CreateHisto(string name, int nbins, double nmin, double nmax, string title="", string xtitle="", string ytitle="", bool CenterTitle="", bool rebin=1){
TH1D *temp = new TH1D(name.c_str(), title.c_str(), nbins, nmin, nmax);
if(!title.empty()) temp->SetTitle(title.c_str());
if(!xtitle.empty()) temp->SetXTitle(xtitle.c_str());
if(!ytitle.empty()) temp->SetYTitle(ytitle.c_str());
if(CenterTitle){
temp->GetXaxis()->CenterTitle();
temp->GetYaxis()->CenterTitle();
}
if(rebin) temp->Rebin(rebin);
return temp;
}
