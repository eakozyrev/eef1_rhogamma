#define yield_cxx
#include "yield.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TProfile.h"
#include "TFile.h"
#include "TMatrixTSym.h"
#include "TMath.h"
#include "TVector.h"
#include "TText.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TBenchmark.h"
#include <iostream>
#include <string>


void yield::Loop(std::string histFileName)
{

   if (fChain == 0) return;
   TFile *newfile = new TFile(histFileName.c_str(),"recreate");
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   TH2D *h_P_Ptr = new TH2D("h_P_Ptr","h_P_Ptr",300,-10,10,300,-10,10);
   h_P_Ptr->SetYTitle("P_{tr}, GeV/c");
   h_P_Ptr->SetXTitle("P, GeV/c");   
   TH1D *h_Egamma = new TH1D("h_Egamma","h_Egamma",1000,0,5);
   TH1D *h_ph_ener = new TH1D("h_ph_ener","h_ph_ener",1000,0,5);
   TH1D *h_angle = new TH1D("h_angle","h_angle",100,-1.2,1.2);
   TH1D *h_angle_gamma = new TH1D("h_angle_gamma","h_angle_gamma",100,-1.2,1.2);
   TH1D *h_angle_gamma_b = new TH1D("h_angle_gamma_b","h_angle_gamma_b",100,-1.2,1.2);
   TH1D *h_el_ener = new TH1D("h_el_ener","h_el_ener",1000,0.,9.2);
   TH1D *h_mrho = new TH1D("h_mrho","h_mrho",50,0.55,0.95);
   h_mrho->SetXTitle("m_{#pi^{+}#pi^{-}}, GeV/c^{2}");
   h_mrho->SetYTitle("yields");
   TH1D *h_mrho_b = new TH1D("h_mrho_b","h_mrho_b",50,0.55,0.95);
   h_mrho_b->SetXTitle("m_{#pi^{+}#pi^{-}}, GeV/c^{2}");
   h_mrho_b->SetYTitle("yields");   
   TH1D *h_elth = new TH1D("h_elth","h_elth",200,0.,3.14);
   TH1D *h_posth = new TH1D("h_posth","h_posth",200,0.,3.14);   
   TH1D *h_mf1 = new TH1D("h_mf1","h_mf1",200,0.9,1.4);
   TH1D *h_Q2 = new TH1D("h_Q2","h_Q2",500,-50,50);
   TH1D *h_Q2_b = new TH1D("h_Q2_b","h_Q2_b",500,-50,50);
   TH2D *h_Q2_Ptr = new TH2D("h_Q2_Ptr","h_Q2_Ptr",500,-50,50,700,0,7);
   TH1D *h_Piselect = new TH1D("h_Piselect","h_Piselect",10000,-1000,70000);
   TH2D *h_rho_Z = new TH2D("h_rho_Z","h_rho_Z",4000,-1,20,400,-20,20); 
   
   TTree Tree ("Tree","Tree");
   double mf1,cosa0,Q2,Q2_0,omega,Q2gen;
   Tree.Branch("mf1",&mf1,"mf1/D");
   Tree.Branch("Q2",&Q2,"Q2/D");
   Tree.Branch("Q2gen",&Q2gen,"Q2gen/D");

   TTree Treegen ("Treegen","Treegen");
   Treegen.Branch("Q2gen",&Q2gen,"Q2gen/D");
   Treegen.Branch("Q2_0",&Q2_0,"Q2_0/D");
  
   TLorentzVector PPh[100];
   TLorentzVector Psi[100];
   TLorentzVector Prho[100];
   TLorentzVector Ppi[100];
   TLorentzVector Pf1[100];
   TLorentzVector Pel[100];
   TLorentzVector Psim[100];
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if(jentry%5000 == 0)cout << jentry << endl;
      double sqrts = sqrt(eee*eee - eepx*eepx - eepy*eepy - eepz*eepz ); 
      TLorentzVector Pbeam(0.,0.,0.,sqrts);
      TLorentzVector P_init_el(0.,0.,-sqrts/2.,sqrts/2.);
      TLorentzVector P_init_pos(0.,0.,sqrts/2.,sqrts/2.);

      
      for(int i = 0; i < ngamma; i++){
	PPh[i] = TLorentzVector(gammaenergycm[i]*sqrt(1.-gammacosthcm[i]*gammacosthcm[i])*sin(gammaphicm[i]),
				gammaenergycm[i]*sqrt(1.-gammacosthcm[i]*gammacosthcm[i])*cos(gammaphicm[i]),
				gammaenergycm[i]*gammacosthcm[i],
				gammaenergycm[i]);
      }
      for(int i = 0; i < npsi; i++){
	Psi[i] = TLorentzVector(psip3cm[i]*sqrt(1.-psicosthcm[i]*psicosthcm[i])*sin(psiphicm[i]),
				psip3cm[i]*sqrt(1.-psicosthcm[i]*psicosthcm[i])*cos(psiphicm[i]),
				psip3cm[i]*psicosthcm[i],
				psienergycm[i]);
      }
      for(int i = 0; i < nrho; i++){
	Prho[i] = TLorentzVector(rhop3cm[i]*sqrt(1.-rhocosthcm[i]*rhocosthcm[i])*sin(rhophicm[i]),
				rhop3cm[i]*sqrt(1.-rhocosthcm[i]*rhocosthcm[i])*cos(rhophicm[i]),
				rhop3cm[i]*rhocosthcm[i],
				rhoenergycm[i]);
      }
      for(int i = 0; i < npi; i++){
	Ppi[i] = TLorentzVector(pip3cm[i]*sqrt(1.-picosthcm[i]*picosthcm[i])*sin(piphicm[i]),
				 pip3cm[i]*sqrt(1.-picosthcm[i]*picosthcm[i])*cos(piphicm[i]),
				 pip3cm[i]*picosthcm[i],
				 pienergycm[i]);
      } 
      for(int i = 0; i < netap; i++){
	Pf1[i] = TLorentzVector(etapp3cm[i]*sqrt(1.-etapcosthcm[i]*etapcosthcm[i])*sin(etapphicm[i]),
				 etapp3cm[i]*sqrt(1.-etapcosthcm[i]*etapcosthcm[i])*cos(etapphicm[i]),
				 etapp3cm[i]*etapcosthcm[i],
				 etapenergycm[i]);
      }

      for(int i = 0; i < nel; i++){
	Pel[i] = TLorentzVector(elp3cm[i]*sqrt(1.-elcosthcm[i]*elcosthcm[i])*sin(elphicm[i]),
				 elp3cm[i]*sqrt(1.-elcosthcm[i]*elcosthcm[i])*cos(elphicm[i]),
				 elp3cm[i]*elcosthcm[i],
				 elenergycm[i]);
      }

      for(int i = 0; i < mclen; i++){
	Psim[i] = TLorentzVector(mcp3cm[i]*sqrt(1.-mccosthcm[i]*mccosthcm[i])*cos(mcphicm[i]),
				mcp3cm[i]*sqrt(1.-mccosthcm[i]*mccosthcm[i])*sin(mcphicm[i]),
				 mcp3cm[i]*mccosthcm[i],
				 mcenergycm[i]);
      }

      int nel_n = 3;
      if(mclund[4]==11)nel_n = 4;
      if(mclund[5]==11)nel_n = 5;
      if(mclund[6]==11)nel_n = 6;
      
      int npos = 3;
      if(mclund[4]==-11)npos = 4;
      if(mclund[5]==-11)npos = 5;
      if(mclund[6]==-11)npos = 6;
      
      Q2gen = -(Psim[0] - Psim[nel_n]).Dot(Psim[0] - Psim[nel_n]);
      Q2_0 = -(Psim[1] - Psim[npos]).Dot(Psim[1] - Psim[npos]);
      if(Q2_0 > Q2gen) {
	Q2gen = Q2_0;
	Q2_0 = -(Psim[0] - Psim[nel_n]).Dot(Psim[0] - Psim[nel_n]);
      }
      Treegen.Fill();
      if(L3outdch == 0 && L3outemc == 0)continue;
      if(Bgfmultihadron == 0)continue;
      
      for(int i = 0; i < npsi; i++){

	double diffP = (Pbeam - Psi[i]).P();
	double diffPtr = sqrt(pow((Pbeam - Psi[i]).Px(),2.) + pow((Pbeam - Psi[i]).Py(),2.));
	h_P_Ptr->Fill(diffP,diffPtr);
	if(diffPtr > 0.05)continue;
	if(diffP < 3.3)continue;
	int n_etap = psid2idx[i];
	int n_gamma = etapd2idx[n_etap];
	int n_rho =  etapd1idx[n_etap];
	int n_pip = rhod1idx[n_rho];
	int n_pim = rhod2idx[n_rho];
	int nfermion = psid1idx[i];
	mf1 = etapmass[n_etap];
	if(ntrk > 3 || mf1 < 1.15)continue;
	double photon_en = 0.;
	for(int k = 0; k < ngamma; k++){
	  if(k == n_gamma)continue;
	  photon_en = photon_en + PPh[k].E();
	}
	if(PPh[n_gamma].E() < 0.3)continue;
	
	if(photon_en > 0.5)continue;
	Ppi[n_pip].Boost(-Prho[n_rho].Px()/Prho[n_rho].E(),-Prho[n_rho].Py()/Prho[n_rho].E(),-Prho[n_rho].Pz()/Prho[n_rho].E());
	double angle = (Ppi[n_pip].E()*Prho[n_rho].E() - Ppi[n_pip].Dot(Prho[n_rho]))/Ppi[n_pip].P()/Prho[n_rho].P();
	h_angle->Fill(angle);
	//if(angle < -0.5 || angle > 0.5)continue;
	PPh[n_gamma].Boost(-Pf1[n_etap].Px()/Pf1[n_etap].E(),-Pf1[n_etap].Py()/Pf1[n_etap].E(),-Pf1[n_etap].Pz()/Pf1[n_etap].E());

	double angle_gamma = (Pf1[n_etap].E()*PPh[n_gamma].E() - Pf1[n_etap].Dot(PPh[n_gamma]))/Pf1[n_etap].P()/PPh[n_gamma].P();
	
	h_Piselect->Fill(piselectorsmap[pitrkidx[n_pip]]);
	h_Piselect->Fill(piselectorsmap[pitrkidx[n_pim]]);
	if(piselectorsmap[pitrkidx[n_pip]] < 20000 || piselectorsmap[pitrkidx[n_pim]] < 20000)continue;
	if(muselectorsmap[pitrkidx[n_pip]] > 10000)continue;
        if(muselectorsmap[pitrkidx[n_pim]] > 10000)continue;
	Q2 = pow((P_init_el - Pel[nfermion]).M(),2.);
        if(psid1lund[i] > 0.)Q2 = pow((P_init_pos - Pel[nfermion]).M(),2.);
	if(Q2 < 2)continue;

	if(mf1 > 1.22 && mf1 < 1.3){
	  h_mrho->Fill(rhomass[n_rho]);
	  h_Egamma->Fill(PPh[n_gamma].E());
	  h_ph_ener->Fill(photon_en);
	}
	if(mf1 > 1.34 && mf1 < 1.38)h_mrho_b->Fill(rhomass[n_rho]);
	if(mf1 > 1.15 && mf1 < 1.19)h_mrho_b->Fill(rhomass[n_rho]);
	//	if(rhomass[n_rho] < 0.65 || rhomass[n_rho] > 0.87)continue;


	if(mf1 > 1.22 && mf1 < 1.3)h_angle_gamma->Fill(angle_gamma);
	if(mf1 > 1.34 && mf1 < 1.38)h_angle_gamma_b->Fill(angle_gamma);
	if(mf1 > 1.15 && mf1 < 1.19)h_angle_gamma_b->Fill(angle_gamma);
	if (angle_gamma < -0.6)continue;

        //if(etapmass[n_etap] < 0.93 || etapmass[n_etap] > 0.97)continue;
        if(psid1lund[i] > 0)h_posth->Fill(Pel[nfermion].Theta());
        if(psid1lund[i] < 0)h_elth->Fill(Pel[nfermion].Theta());
//	if(psid1lund[i] > 0 && Pel[nfermion].Theta() > 0.8)continue;
//      if(psid1lund[i] < 0 && Pel[nfermion].Theta() < 2.4)continue;
	
	if(mf1 > 1.22 && mf1 < 1.3)h_Q2->Fill(Q2);
	if(mf1 > 1.34 && mf1 < 1.38)h_Q2_b->Fill(Q2);
	if(mf1 > 1.15 && mf1 < 1.19)h_Q2_b->Fill(Q2);
	h_el_ener->Fill(Pel[nfermion].E());
	h_mf1->Fill(etapmass[n_etap]);
	h_Q2_Ptr->Fill(Q2,diffP);
	h_rho_Z->Fill(Trkdocaxy_xy[pitrkidx[n_pip]],Trkdocaxy_z[pitrkidx[n_pip]]);
	h_rho_Z->Fill(Trkdocaxy_xy[pitrkidx[n_pim]],Trkdocaxy_z[pitrkidx[n_pim]]);
	Tree.Fill();
      }

      
   }
   newfile->Write();
}
