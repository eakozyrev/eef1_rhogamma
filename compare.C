#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TProfile.h"
#include "TFile.h"
#include "TMatrixTSym.h"
#include "TGraphErrors.h"
#include "TGaxis.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TVector.h"
#include "TStyle.h"
#include "TText.h"
#include "TBenchmark.h"
#include <iostream>
#include <fstream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "THStack.h"

void compare(string filee, string filemc, string variable,string cutt = "1==1"){

  TCut cut = cutt.c_str();
  TFile *newfile = TFile::Open((filee).c_str());
  TTree* tree = (TTree*)newfile->Get("Tree");
  TH1D *htest = new TH1D("htest","htest",1000000,-10000,10000);
  tree->Draw((variable + " >> htest").c_str(),cut);
  double start = -10000;
  double end = 100000;
  for(int i = 1; i <= 1000000; i++){
	  
	  if(htest->GetBinContent(i) > 0)start = htest->GetBinCenter(i);
	  
  }
  for(int i = 1000000; i > 0; i--){
	  if(htest->GetBinContent(i) > 0)end = htest->GetBinCenter(i);
  }
  TH1D *hm2ph = new TH1D("hm2ph","hm2ph",100,start,end);
  tree->Draw((variable + " >> hm2ph").c_str(),cut);
  double notm = tree->GetEntries(cut);
  cout << notm << endl;
  TFile *newfilemc = TFile::Open(filemc.c_str());
  TTree* treemc = (TTree*)newfilemc->Get("Tree");
  TH1D *hm2ph_mc = new TH1D("hm2ph_mc","hm2ph_mc",100,start,end);
  treemc->Draw((variable + " >> hm2ph_mc").c_str(),cut);
  hm2ph_mc->SetNormFactor(notm);
  hm2ph->Draw();
  hm2ph->SetLineWidth(2.);
  hm2ph->SetLineColor(4);
  hm2ph_mc->SetFillColor(4);
  hm2ph_mc->SetFillStyle(3002);
  hm2ph_mc->Draw("same");
  //compare("scan2011_tr_ph_fc_e750_v5.root", "../4pisel/histograms/mc_750_10156.root")
}















