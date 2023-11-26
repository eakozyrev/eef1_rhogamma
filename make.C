#include "TChain.h"
#include "yield.C+"


int main(){
 
    TChain *chain= new TChain("h1");
/*    chain->Add("trees/run1.root");
    chain->Add("trees/run2.root");
    chain->Add("trees/run3.root");
    chain->Add("trees/run4.root");
    chain->Add("trees/run5.root");    
    chain->Add("trees/run6.root"); */
    chain->Add("trees/data.root");
    yield as(chain);
    cout << "trees/data.root is opened "  << endl;
    as.Loop("run2.root");
    
    return 1;
}

int main_mc_m0(){
 
    TChain *chain= new TChain("h1");
    chain->Add("trees/mc/m0/mc.root");
    yield as(chain);
    cout << "trees/mc/m0/mc.root is opened "  << endl;
    as.Loop("mc_m0.root");
    
    return 1;
}

int main_mc_m1(){
 
    TChain *chain= new TChain("h1");
    chain->Add("trees/mc/m1/mc.root");
    yield as(chain);
    cout << "trees/mc/m1/mc.root is opened "  << endl;
    as.Loop("mc_m1.root");
    
    return 1;
}


int main_mc_etap(){
 
    TChain *chain= new TChain("h1");
    chain->Add("../etaprime_rhogamma/mc/mc.root");
    yield as(chain);
    cout << "../etaprime_rhogamma/mc/mc.root is opened "  << endl;
    as.Loop("mc_etap.root");
    
    return 1;
}



void make(){

main();
main_mc_m0();
main_mc_m1();
//main_mc_etap();
}


void draw(){

  TLine l(10.26865,1e-8,10.26865,100);
  l.Draw();
  TLine l1(10.25546,1e-8,10.25546,100);
  l1.Draw();
  TLine l2(10.2325,1e-8,10.2325,100);
  l2.Draw();

}

void handle(){

 TFile *_file0 = TFile::Open("hdata.root");;
 TH1D *h2pi_exp = (TH1D*)_file0->Get("h2pi");
 double norm1 = 0;
 for(int i = 0; i < 100; i++){norm1 = norm1 + h2pi_exp->GetBinContent(i);}
 cout << "norm1 = " << norm1 << endl;
 h2pi_exp->Rebin(10);
 
 TFile *_file1 = TFile::Open("hdata_off4S.root");;
 TH1D *h2pi_off = (TH1D*)_file1->Get("h2pi");
 double norm2 = 0;
 for(int i = 0; i < 100; i++){norm2 = norm2 + h2pi_off->GetBinContent(i);}
 cout << "norm2 = " << norm1 << endl;
 h2pi_off->Rebin(10);
 double scale = 27.96/43.92;
 h2pi_exp->Add(h2pi_off,-277./253);
 h2pi_exp->SetLineColor(2);
 h2pi_exp->Draw();
 //h2pi_off->SetNormFactor(h2pi_exp->GetEntries()); 
 //h2pi_off->Draw("same");

}
