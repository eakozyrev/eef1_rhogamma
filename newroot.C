#define newroot_cxx
#include "newroot.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"

void newroot::Loop()
{
//   In a ROOT session, you can do:
//      root> .L newroot.C
//      root> newroot t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   TFile newfile("trees/data.root","recreate");
   TTree *newtree = fChain->GetTree()->CloneTree(0);
   TLorentzVector Psi[100];
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      double sqrts = sqrt(eee*eee - eepx*eepx - eepy*eepy - eepz*eepz ); 
      TLorentzVector Pbeam(0.,0.,0.,sqrts);
      
      bool writee = false;
      for(int i = 0; i < netap; i++){
	if(etapmass[i] > 1.1) writee = true;
      }

      for(int i = 0; i < npsi; i++){
	Psi[i] = TLorentzVector(psip3cm[i]*sqrt(1.-psicosthcm[i]*psicosthcm[i])*sin(psiphicm[i]),
				psip3cm[i]*sqrt(1.-psicosthcm[i]*psicosthcm[i])*cos(psiphicm[i]),
				psip3cm[i]*psicosthcm[i],
				psienergycm[i]);
      }
      bool write1 = false;
      for(int i = 0; i < npsi; i++){
	double diffP = (Pbeam - Psi[i]).P();
	if(diffP > 3.)write1 = true;
      }
      
      if(writee == true && write1 == true)newtree->Fill();

   }
   newfile.Write();
}
