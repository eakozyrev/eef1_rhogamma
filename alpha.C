#include <TH2.h>
#include <TF1.h>
#include <TH1.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TProfile.h"
#include "TFile.h"
#include "TText.h"
#include "TMatrixTSym.h"
#include "TGraphErrors.h"
#include "TGaxis.h"
#include "TMath.h"
#include "TVector.h"
#include "TStyle.h"
#include "TBenchmark.h"
#include <iostream>
#include <fstream>
#include "RooChi2Var.h"
#include <string>
#include "TCut.h"
#include "TFile.h"
#include "TTree.h"
#include "THStack.h"
#include "RooPolynomial.h"
#include "TEfficiency.h"

#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

#define chi25Ccut 25.
#include <stdlib.h>

using namespace std;


#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooConstVar.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "RooKeysPdf.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
using namespace RooFit;


double paramm[100];

inline bool exists_file(const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}


double fit2g1(string filee, string filemc, double qmin, double qmax){

  cout <<  filee << " " <<  filemc << endl;
  // Observable
  double startt = qmin;//startt, endd
  double endd = qmax;
  RooRealVar x("mf1","m_{#pi^{+}#pi^{-}#gamma}, GeV/c^{2}",startt, endd);
  RooRealVar frac("frac","frac",0.9,0.040,1.0);

  TFile *newfile = TFile::Open((""+filee).c_str());
  TTree* tree = (TTree*)newfile->Get("Tree");

   x.setBins(50); 
   RooRealVar cosa0("cosa0","cosa0",-4,4); 
   RooRealVar Q2("Q2","Q2",0,400); 

   string selection = Form("mf1 > %g && mf1 < %g && Q2 > %g && Q2 < %g",startt,endd, 0., 100.);

   RooDataSet *data = new RooDataSet("data","data",RooArgSet(x,Q2),Import(*tree),Cut(selection.c_str()));
   RooKeysPdf kest_exp("kest_exp","kest_exp",x,*data,RooKeysPdf::MirrorBoth);
   
  //*************************************************************************
  //===========================MC============================================

  TFile *newfilemc = TFile::Open((""+filemc).c_str());
  TTree* treemc = (TTree*)newfilemc->Get("Tree");

  string selectionmc = Form("mf1 > %g && mf1 < %g && Q2 > %g-5 && Q2 < %g+2",startt,endd, 0., 100.); 
  RooDataSet datamc("datamc","datamc",RooArgSet(x,Q2),Import(*treemc),Cut(selectionmc.c_str()));
  RooKeysPdf kest1("kest1","kest1",x,datamc,RooKeysPdf::MirrorBoth) ;
  RooRealVar mg("mg","mg",-0.01,-0.05,0.05); 
  RooRealVar sg("sg","sg",0.003);//,0.0000001,0.1); 
  RooGaussian gauss("gauss","gauss",x,mg,sg);
  cout << "===========================================================================9" << endl;
  RooFFTConvPdf resol("lxg","landau (X) gauss",x,kest1,gauss);
  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++9" << endl;
  //**************************************************************************
  //========================== BCKGR ===========================================
  /*
  TFile *newfileBCKGR = TFile::Open(("../4pisel/histograms/"+filebkgr).c_str());
  TTree* treeBCKGR = (TTree*)newfileBCKGR->Get("Tree");
  RooDataSet dataBCKGR("dataBCKGR","dataBCKGR",RooArgSet(RooArgSet(x,helicity1,helicity2,chi25C,tthg1, tthg2, tthg3, tthg4),RooArgSet(tthpic1,tthpic2,momg1,momg2,momg3,momg4)),Import(*treeBCKGR),Cut(selection.c_str()));
  RooKeysPdf kestBCKGR("kestBCKGR","kestBCKGR",x,dataBCKGR,RooKeysPdf::MirrorBoth);
  */
  //**************************************************************************
  //==========================MODEL===========================================
  RooRealVar c0("c0","c0",10.,0.,30.); 
  RooRealVar c1("c1","c1",0.,-10.,10.); 
  RooRealVar c2("c2","c2",0.0);//,-0.001,0.001); 
  RooGenericPdf P("P","c0 + c1*(mf1-1.28) + c2*(mf1-1.28)*(mf1-1.28)",RooArgSet(x,c0,c1,c2)); 
  RooPolynomial pol0("pol0","pol0",x,RooArgList());
  
  RooAddPdf model("model","model",RooArgList(resol,P),frac) ;
  //model.plotOn(xframe1);
  // Construct unbinned likelihood of model w.r.t. data
  RooAbsReal* nll = model.createNLL(*data) ;

  // I n t e r a c t i v e   m i n i m i z a t i o n ,   e r r o r   a n a l y s i s
  // -------------------------------------------------------------------------------

  // Create MINUIT interface object
  //RooMinuit m(*nll) ;

  // Activate verbose logging of MINUIT parameter space stepping
  //m.setVerbose(kTRUE) ;
  //m.setVerbose(kFALSE) ;

  // Call MIGRAD to minimize the likelihood
  RooFitResult* r = model.fitTo(*data,Save());
  //m.migrad();
  cout << "-log(L) at minimum = " << r->minNll() << endl ;
  cout << " (double)frac.getVal() = " <<  (double)frac.getVal() << endl;
 
  
  paramm[0] = (double)frac.getVal()*tree->GetEntries(TCut(selection.c_str()));
  paramm[1] = frac.getError()*tree->GetEntries(TCut(selection.c_str()));

  paramm[2] = (double)mg.getVal();
  paramm[3] = mg.getError();

  paramm[4] = (double)sg.getVal();
  paramm[5] = sg.getError();
  
  // Print values of all parameters, that reflect values (and error estimates)
  // that are back propagated from MINUIT
  //model.getParameters(x)->Print("s") ;
  RooPlot* xframe1 = x.frame(Title(Form("%g < Q^{2} < %g GeV^{2}",qmin,qmax)));
  data->plotOn(xframe1);
  model.plotOn(xframe1);
  cout << "frame.chiSquare() = " << xframe1->chiSquare() << endl;
  model.plotOn(xframe1,Components(P),LineStyle(kDashed));
  TCanvas *s = new TCanvas();
  TPad*    upperPad = new TPad("upperPad", "upperPad", 0.0,0.05,1.0,1.0);
  TPad*    lowerPad = new TPad("lowerPad", "lowerPad", 0.0,0.1,1.0,0.3);
  upperPad->Draw();
  //lowerPad->Draw();
  upperPad->cd();
  //gPad->SetLeftMargin(0.15) ; 
  xframe1->GetYaxis()->SetTitleOffset(0.7);
  xframe1->SetNdivisions(8,"Y");
  xframe1->GetYaxis()->SetLabelOffset(0.007); 
  xframe1->Draw();
  
  cout << "===============================================" << endl;
  cout << "                     end                       " << endl;
  cout << "===============================================" << endl;
  s->SaveAs(Form("figs/%g_%g.png",qmin,qmax));
  cout << "============================================ " << endl;
  return 1.;

}


void fitall(){


  // paramm[0] +/- paramm[1]
  
  TFile *newfilemc = TFile::Open("mc_m0.root");
  TTree* treemc = (TTree*)newfilemc->Get("Tree");
  TFile *newfilemc0 = TFile::Open("trees/mc/m0/mc.root");
  TTree* treemc0 = (TTree*)newfilemc0->Get("h1");

  double effmc0 = (double)treemc->GetEntries()/treemc0->GetEntries();
  double deffmc0 = sqrt(effmc0*(1.-effmc0)/treemc0->GetEntries());


  TFile *newfilemc1 = TFile::Open("mc_m1.root");
  TTree* treemc1 = (TTree*)newfilemc1->Get("Tree");
  TFile *newfilemc01 = TFile::Open("trees/mc/m1/mc.root");
  TTree* treemc01 = (TTree*)newfilemc01->Get("h1");

  double effmc1 = (double)treemc1->GetEntries()/treemc01->GetEntries();
  double deffmc1 = sqrt(effmc1*(1.-effmc1)/treemc01->GetEntries());
 

  fit2g1("run2.root","mc_m1.root",1.15,1.4);
  double Nf1 = paramm[0];
  double dNf1 = paramm[1];
 
  
  fit2g1("run2.root","mc_etap.root",0.9,1.03);
  double Netap = paramm[0];
  double dNetap = paramm[1];

  
  
  TFile *newfile = TFile::Open("mc_etap.root");
  TTree* tree = (TTree*)newfile->Get("Tree");
  TFile *newfil0 = TFile::Open("../etaprime_rhogamma/mc/mc.root");
  TTree* tree0 = (TTree*)newfil0->Get("h1");
  
  double eff = (double)tree->GetEntries()/tree0->GetEntries();
  double deff = sqrt(eff*(1.-eff)/tree0->GetEntries());
  
  cout << "effmc0 = " << effmc0 << " +/- " << deffmc0 <<  endl;
  cout << "effmc1 = " << effmc1 << " +/- " << deffmc1 <<  endl;
  cout << "Nf1 = " << Nf1 << " +/- " << dNf1 << endl;
  cout << "Netap = " << Netap << " +/- " << dNetap << endl;  
  cout << "eff etap = " << eff << " +/- " << deff << endl;
  
  
}



void effic_draw(string file, string same = "", string title = ""){
  TFile *newfilem0 = TFile::Open(file.c_str());
  TTree* treem0 = (TTree*)newfilem0->Get("Tree");
  TH1D *h_Q2 = new TH1D("h_Q2","h_Q2",40,1,20);
  h_Q2->SetXTitle("Q^{2} (GeV)^{2}");
  TCanvas s;
  treem0->Draw("Q2 >> h_Q2");
  TTree* treem0gen = (TTree*)newfilem0->Get("Treegen");
  TH1D *h_Q2gen = new TH1D("h_Q2gen","h_Q2gen",40,1,20);
  treem0gen->Draw("Q2gen >> h_Q2gen");
  s.Close();
  TEfficiency *hEfficiency = new TEfficiency(*h_Q2,*h_Q2gen);
  hEfficiency->SetTitle((title + " ; Q^{2} (GeV^{2}) ; #varepsilon").c_str());
  hEfficiency->Draw(same.c_str());

}

vector<double*> vector_effic(string file){
  TFile *newfilem0 = TFile::Open(file.c_str());
  TTree* treem0 = (TTree*)newfilem0->Get("Tree");
  TTree* treem0gen = (TTree*)newfilem0->Get("Treegen");
  double val[] = {2,4,5,6,7,10,20};
  vector<double*> d1;
  for(int i = 0; i < sizeof(val)/sizeof(1.)-1; i++){
    double *arr = new double[2];
    TCut cut0 = Form("Q2gen > %g && Q2gen < %g", val[i], val[i+1]);
    TCut cut = Form("Q2 > %g && Q2 < %g", val[i], val[i+1]);
    double eff_ = (double)treem0->GetEntries(cut)/(double)treem0gen->GetEntries(cut0);
    double deff_ = sqrt(eff_*(1.-eff_)/(double)treem0gen->GetEntries(cut0));
    cout << val[i] << " " << val[i+1] << " " << eff_ << endl;
    arr[0] = eff_;
    arr[1] = deff_;
    d1.push_back(arr);
  }

  return d1;
  
}

double fm(double x){
  return ((int)(x*10000))/10000.;
}

double calculate_R(){

  fit2g1("run2.root","mc_m1.root",1.15,1.4);
  
  ifstream stream_cross0("../eef1/results/cross0.dat");
  ifstream stream_cross1("../eef1/results/cross0.dat");
  ifstream stream_Nf1("../eef1/results/nevents.dat");
  double val[] = {2,4,5,6,7,10,20};
  double Nf1_pipieta = 0;
  double dNf1_pipieta = 0;
  vector<double*> eff_rhog_0 = vector_effic("mc_m0.root");
  vector<double*>eff_rhog_1 = vector_effic("mc_m1.root");
  for(int i = 0; i < sizeof(val)/sizeof(1.)-1; i++){

    cout << val[i] << " \\div " << val[i+1] << " & " << fm(eff_rhog_0[i][0]) << " \\pm " <<  fm(eff_rhog_0[i][1]) << " & " << fm(eff_rhog_1[i][0]) << " \\pm " <<  fm(eff_rhog_1[i][1]) << "\\\\ " << endl;

  }

  double sigma0 = 0;
  double sigma1 = 0;

  double eff_pipieta = 0;
  double eff_rho_gamma = 0;

  double *eff_pipieta_0 = new double[sizeof(val)/sizeof(1.)-1];
  double *eff_pipieta_1 = new double[sizeof(val)/sizeof(1.)-1];

  double *sigma_0 = new double[sizeof(val)/sizeof(1.)-1];
  double *sigma_1 = new double[sizeof(val)/sizeof(1.)-1];

  double eff_rhog=0;
  double eff_2pieta=0;
  
  for(int i = 0; i < sizeof(val)/sizeof(1.)-1; i++){
    
    double btv;
    stream_Nf1 >> btv >> btv >> btv;
    Nf1_pipieta +=btv;
    stream_Nf1 >> btv;
    dNf1_pipieta += btv*btv;

    stream_cross0 >> btv >> btv >> btv >> btv >> btv;
    eff_pipieta_0[i] = btv;
    stream_cross0 >> btv >> btv;
    sigma_0[i] = btv;
    stream_cross0 >> btv >> btv >> btv;

    stream_cross1 >> btv >> btv >> btv >> btv >> btv;
    eff_pipieta_1[i] = btv;
    stream_cross1 >> btv >> btv;
    sigma_1[i] = btv;
    stream_cross1  >> btv >> btv >> btv;
    
  }
  dNf1_pipieta = sqrt(dNf1_pipieta);
  for(int i = 0; i < sizeof(val)/sizeof(1.)-1; i++){
    sigma0+=sigma_0[i];
    sigma1+=sigma_1[i];
  }

  double Nf1 = paramm[0];
  double dNf1 = paramm[1];
 
  for(int i = 0; i < sizeof(val)/sizeof(1.)-1; i++){

    eff_rhog+= eff_rhog_0[i][0]*sigma_0[i] + eff_rhog_1[i][0]*sigma_1[i];
    eff_2pieta+=eff_pipieta_0[i]*sigma_0[i] + eff_pipieta_1[i]*sigma_1[i];
    
  }

  eff_rhog/=(sigma0+sigma1);
  eff_2pieta/=(sigma0+sigma1);

  cout << Nf1_pipieta << " " << eff_2pieta << " " << Nf1 << " " << eff_rhog << endl;
  double R = Nf1_pipieta/eff_2pieta/Nf1*eff_rhog;
  //R /= 0.3941;
  R*=1.5;

  double dR = pow(dNf1/Nf1,2.);// + pow(dNf1_pipieta/Nf1_pipieta,2.);
  dR = sqrt(dR);

  cout << "B(f1-> pipieta) / B(f1 -> rho gamma) = " << R << " +/- " << dR*R << endl; 
  return R;
}



double form_factor_1(double Q2){

  return 1./(1.+Q2/0.775/0.775);

}




void compare_hist(string hist = "h_mrho", string file1 = "run2.root", string hfile2 = "mc_m0.root"){

    TFile *file = TFile::Open(("" + file1).c_str());
    TH1D *hist1 = (TH1D*)file->Get(hist.c_str());

    TFile *file2 = TFile::Open(("" + hfile2).c_str());
    TH1D *hist2 = (TH1D*)file2->Get(hist.c_str());
     
    hist2->SetNormFactor(hist1->GetEntries());
    //hist2->Rebin(6);
    //hist1->Rebin(6);
    //hist1->SetXTitle(hist.c_str());
    //hist2->SetXTitle(hist.c_str());
    hist1->Draw();
    hist1->SetLineWidth(2.);
    hist2->SetLineColor(4);
    hist2->SetFillColor(4);
    hist2->SetFillStyle(3002);
    hist2->Draw("same");
    
    
}


void compare_hist_1(string hist = "h_mrho", string histb = "h_mrho_b", string file1 = "run2.root", string hfile2 = "mc_m0.root"){

    TFile *file = TFile::Open(("" + file1).c_str());
    TH1D *hist1 = (TH1D*)file->Get(hist.c_str());
    TH1D *hist1_b = (TH1D*)file->Get(histb.c_str());
 
    TFile *file2 = TFile::Open(("" + hfile2).c_str());
    TH1D *hist2 = (TH1D*)file2->Get(hist.c_str());
     
    hist2->Scale((hist1->GetEntries()-hist1_b->GetEntries())/(double)hist2->GetEntries());
    hist2->Add(hist1_b,1.);
    //hist2->Rebin(6);
    //hist1->Rebin(6);
    //hist1->SetXTitle(hist.c_str());
    //hist2->SetXTitle(hist.c_str());
    hist1->Draw("e");
    hist1->SetLineWidth(2.);
    hist2->SetLineColor(4);
    hist2->SetFillColor(4);
    hist2->SetFillStyle(3002);
    hist2->Draw("histsame");
    hist1_b->Draw("same");
    
}
