//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May 16 23:48:31 2019 by ROOT version 6.15/01
// from TTree h1/myNtuple
// found on file: trees/run1.root
//////////////////////////////////////////////////////////

#ifndef yield_h
#define yield_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class yield {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           runnumber;
   Int_t           platform;
   Int_t           partition;
   Int_t           upperid;
   Int_t           lowerid;
   Int_t           majorid;
   Int_t           configkey;
   Int_t           date;
   Int_t           ddate;
   Float_t         eepx;
   Float_t         eepy;
   Float_t         eepz;
   Float_t         eee;
   UChar_t         Bgfmultihadron;
   UChar_t         Bgfneutralhadron;
   UChar_t         Bgfmumu;
   UChar_t         Bgftau;
   UChar_t         Bgftwoprong;
   UChar_t         Bgfphigamma;
   UChar_t         Bgfallneutraltwophoton;
   UChar_t         Bgfisr;
   UChar_t         Bgfradtwoprong;
   UChar_t         Bgfhighmasshadron;
   UChar_t         Bgftwophotontwotrack;
   UChar_t         Digifdchemcpreveto;
   UChar_t         Digifl1open;
   UChar_t         Digifl3open;
   UChar_t         L3outdch;
   UChar_t         L3outemc;
   Int_t           ntracks;
   Int_t           ngoodtrkloose;
   Int_t           mclen;
   Int_t           mclund[200];   //[mclen]
   Int_t           mothidx[200];   //[mclen]
   Int_t           daulen[200];   //[mclen]
   Int_t           dauidx[200];   //[mclen]
   Float_t         mcmass[200];   //[mclen]
   Float_t         mcp3cm[200];   //[mclen]
   Float_t         mccosthcm[200];   //[mclen]
   Float_t         mcphicm[200];   //[mclen]
   Float_t         mcenergycm[200];   //[mclen]
   Float_t         mcp3[200];   //[mclen]
   Float_t         mccosth[200];   //[mclen]
   Float_t         mcphi[200];   //[mclen]
   Float_t         mcenergy[200];   //[mclen]
   Int_t           npsi;
   Float_t         psimass[50];   //[npsi]
   Float_t         psimasserr[50];   //[npsi]
   Float_t         psicosth[50];   //[npsi]
   Float_t         psicosthcm[50];   //[npsi]
   Float_t         psienergy[50];   //[npsi]
   Float_t         psienergycm[50];   //[npsi]
   Float_t         psip3[50];   //[npsi]
   Float_t         psip3cm[50];   //[npsi]
   Float_t         psiphi[50];   //[npsi]
   Float_t         psiphicm[50];   //[npsi]
   Int_t           psilund[50];   //[npsi]
   Int_t           psid1lund[50];   //[npsi]
   Int_t           psid1idx[50];   //[npsi]
   Int_t           psid2lund[50];   //[npsi]
   Int_t           psid2idx[50];   //[npsi]
   Int_t           nel;
   Float_t         eldocaxy_xy[50];   //[nel]
   Float_t         eldocaxy_xyerr[50];   //[nel]
   Float_t         eldocaxy_z[50];   //[nel]
   Float_t         eldocaxy_zerr[50];   //[nel]
   Float_t         elcosth[50];   //[nel]
   Float_t         elcosthcm[50];   //[nel]
   Float_t         elenergy[50];   //[nel]
   Float_t         elenergycm[50];   //[nel]
   Float_t         elp3[50];   //[nel]
   Float_t         elp3cm[50];   //[nel]
   Float_t         elphi[50];   //[nel]
   Float_t         elphicm[50];   //[nel]
   Int_t           ellund[50];   //[nel]
   Int_t           eld1lund[50];   //[nel]
   Int_t           eld1idx[50];   //[nel]
   Int_t           eld2lund[50];   //[nel]
   Int_t           eld2idx[50];   //[nel]
   Int_t           eld3lund[50];   //[nel]
   Int_t           eld3idx[50];   //[nel]
   Int_t           eld4lund[50];   //[nel]
   Int_t           eld4idx[50];   //[nel]
   Int_t           eltrkidx[50];   //[nel]
   Int_t           netap;
   Float_t         etapmass[50];   //[netap]
   Float_t         etapmasserr[50];   //[netap]
   Float_t         etapcosth[50];   //[netap]
   Float_t         etapcosthcm[50];   //[netap]
   Float_t         etapenergy[50];   //[netap]
   Float_t         etapenergycm[50];   //[netap]
   Float_t         etapp3[50];   //[netap]
   Float_t         etapp3cm[50];   //[netap]
   Float_t         etapphi[50];   //[netap]
   Float_t         etapphicm[50];   //[netap]
   Int_t           etaplund[50];   //[netap]
   Int_t           etapd1lund[50];   //[netap]
   Int_t           etapd1idx[50];   //[netap]
   Int_t           etapd2lund[50];   //[netap]
   Int_t           etapd2idx[50];   //[netap]
   Int_t           etapd3lund[50];   //[netap]
   Int_t           etapd3idx[50];   //[netap]
   Int_t           nrho;
   Float_t         rhomass[50];   //[nrho]
   Float_t         rhomasserr[50];   //[nrho]
   Float_t         rhocosth[50];   //[nrho]
   Float_t         rhocosthcm[50];   //[nrho]
   Float_t         rhoenergy[50];   //[nrho]
   Float_t         rhoenergycm[50];   //[nrho]
   Float_t         rhop3[50];   //[nrho]
   Float_t         rhop3cm[50];   //[nrho]
   Float_t         rhophi[50];   //[nrho]
   Float_t         rhophicm[50];   //[nrho]
   Int_t           rholund[50];   //[nrho]
   Int_t           rhod1lund[50];   //[nrho]
   Int_t           rhod1idx[50];   //[nrho]
   Int_t           rhod2lund[50];   //[nrho]
   Int_t           rhod2idx[50];   //[nrho]
   Int_t           npi;
   Float_t         pidocaxy_xy[50];   //[npi]
   Float_t         pidocaxy_xyerr[50];   //[npi]
   Float_t         pidocaxy_z[50];   //[npi]
   Float_t         pidocaxy_zerr[50];   //[npi]
   Float_t         picosth[50];   //[npi]
   Float_t         picosthcm[50];   //[npi]
   Float_t         pienergy[50];   //[npi]
   Float_t         pienergycm[50];   //[npi]
   Float_t         pip3[50];   //[npi]
   Float_t         pip3cm[50];   //[npi]
   Float_t         piphi[50];   //[npi]
   Float_t         piphicm[50];   //[npi]
   Int_t           pilund[50];   //[npi]
   Int_t           pitrkidx[50];   //[npi]
   Int_t           ngamma;
   Float_t         gammacosth[50];   //[ngamma]
   Float_t         gammacosthcm[50];   //[ngamma]
   Float_t         gammaenergy[50];   //[ngamma]
   Float_t         gammaenergycm[50];   //[ngamma]
   Float_t         gammap3[50];   //[ngamma]
   Float_t         gammap3cm[50];   //[ngamma]
   Float_t         gammaphi[50];   //[ngamma]
   Float_t         gammaphicm[50];   //[ngamma]
   Int_t           gammalund[50];   //[ngamma]
   Int_t           ntrk;
   Float_t         Trkdocaxy_xy[50];   //[ntrk]
   Float_t         Trkdocaxy_xyerr[50];   //[ntrk]
   Float_t         Trkdocaxy_z[50];   //[ntrk]
   Float_t         Trkdocaxy_zerr[50];   //[ntrk]
   Float_t         Trkcosth[50];   //[ntrk]
   Float_t         Trkcosthcm[50];   //[ntrk]
   Float_t         Trkenergy[50];   //[ntrk]
   Float_t         Trkenergycm[50];   //[ntrk]
   Float_t         Trkp3[50];   //[ntrk]
   Float_t         Trkp3cm[50];   //[ntrk]
   Float_t         Trkphi[50];   //[ntrk]
   Float_t         Trkphicm[50];   //[ntrk]
   Int_t           Trklund[50];   //[ntrk]
   Int_t           piselectorsmap[50];   //[ntrk]
   Int_t           eselectorsmap[50];   //[ntrk]
   Int_t           Kselectorsmap[50];   //[ntrk]
   Int_t           muselectorsmap[50];   //[ntrk]

   // List of branches
   TBranch        *b_runnumber;   //!
   TBranch        *b_platform;   //!
   TBranch        *b_partition;   //!
   TBranch        *b_upperid;   //!
   TBranch        *b_lowerid;   //!
   TBranch        *b_majorid;   //!
   TBranch        *b_configkey;   //!
   TBranch        *b_date;   //!
   TBranch        *b_ddate;   //!
   TBranch        *b_eepx;   //!
   TBranch        *b_eepy;   //!
   TBranch        *b_eepz;   //!
   TBranch        *b_eee;   //!
   TBranch        *b_Bgfmultihadron;   //!
   TBranch        *b_Bgfneutralhadron;   //!
   TBranch        *b_Bgfmumu;   //!
   TBranch        *b_Bgftau;   //!
   TBranch        *b_Bgftwoprong;   //!
   TBranch        *b_Bgfphigamma;   //!
   TBranch        *b_Bgfallneutraltwophoton;   //!
   TBranch        *b_Bgfisr;   //!
   TBranch        *b_Bgfradtwoprong;   //!
   TBranch        *b_Bgfhighmasshadron;   //!
   TBranch        *b_Bgftwophotontwotrack;   //!
   TBranch        *b_Digifdchemcpreveto;   //!
   TBranch        *b_Digifl1open;   //!
   TBranch        *b_Digifl3open;   //!
   TBranch        *b_L3outdch;   //!
   TBranch        *b_L3outemc;   //!
   TBranch        *b_ntracks;   //!
   TBranch        *b_ngoodtrkloose;   //!
   TBranch        *b_mclen;   //!
   TBranch        *b_mclund;   //!
   TBranch        *b_mothidx;   //!
   TBranch        *b_daulen;   //!
   TBranch        *b_dauidx;   //!
   TBranch        *b_mcmass;   //!
   TBranch        *b_mcp3cm;   //!
   TBranch        *b_mccosthcm;   //!
   TBranch        *b_mcphicm;   //!
   TBranch        *b_mcenergycm;   //!
   TBranch        *b_mcp3;   //!
   TBranch        *b_mccosth;   //!
   TBranch        *b_mcphi;   //!
   TBranch        *b_mcenergy;   //!
   TBranch        *b_npsi;   //!
   TBranch        *b_psimass;   //!
   TBranch        *b_psimasserr;   //!
   TBranch        *b_psicosth;   //!
   TBranch        *b_psicosthcm;   //!
   TBranch        *b_psienergy;   //!
   TBranch        *b_psienergycm;   //!
   TBranch        *b_psip3;   //!
   TBranch        *b_psip3cm;   //!
   TBranch        *b_psiphi;   //!
   TBranch        *b_psiphicm;   //!
   TBranch        *b_psilund;   //!
   TBranch        *b_psid1lund;   //!
   TBranch        *b_psid1idx;   //!
   TBranch        *b_psid2lund;   //!
   TBranch        *b_psid2idx;   //!
   TBranch        *b_nel;   //!
   TBranch        *b_eldocaxy_xy;   //!
   TBranch        *b_eldocaxy_xyerr;   //!
   TBranch        *b_eldocaxy_z;   //!
   TBranch        *b_eldocaxy_zerr;   //!
   TBranch        *b_elcosth;   //!
   TBranch        *b_elcosthcm;   //!
   TBranch        *b_elenergy;   //!
   TBranch        *b_elenergycm;   //!
   TBranch        *b_elp3;   //!
   TBranch        *b_elp3cm;   //!
   TBranch        *b_elphi;   //!
   TBranch        *b_elphicm;   //!
   TBranch        *b_ellund;   //!
   TBranch        *b_eld1lund;   //!
   TBranch        *b_eld1idx;   //!
   TBranch        *b_eld2lund;   //!
   TBranch        *b_eld2idx;   //!
   TBranch        *b_eld3lund;   //!
   TBranch        *b_eld3idx;   //!
   TBranch        *b_eld4lund;   //!
   TBranch        *b_eld4idx;   //!
   TBranch        *b_eltrkidx;   //!
   TBranch        *b_netap;   //!
   TBranch        *b_etapmass;   //!
   TBranch        *b_etapmasserr;   //!
   TBranch        *b_etapcosth;   //!
   TBranch        *b_etapcosthcm;   //!
   TBranch        *b_etapenergy;   //!
   TBranch        *b_etapenergycm;   //!
   TBranch        *b_etapp3;   //!
   TBranch        *b_etapp3cm;   //!
   TBranch        *b_etapphi;   //!
   TBranch        *b_etapphicm;   //!
   TBranch        *b_etaplund;   //!
   TBranch        *b_etapd1lund;   //!
   TBranch        *b_etapd1idx;   //!
   TBranch        *b_etapd2lund;   //!
   TBranch        *b_etapd2idx;   //!
   TBranch        *b_etapd3lund;   //!
   TBranch        *b_etapd3idx;   //!
   TBranch        *b_nrho;   //!
   TBranch        *b_rhomass;   //!
   TBranch        *b_rhomasserr;   //!
   TBranch        *b_rhocosth;   //!
   TBranch        *b_rhocosthcm;   //!
   TBranch        *b_rhoenergy;   //!
   TBranch        *b_rhoenergycm;   //!
   TBranch        *b_rhop3;   //!
   TBranch        *b_rhop3cm;   //!
   TBranch        *b_rhophi;   //!
   TBranch        *b_rhophicm;   //!
   TBranch        *b_rholund;   //!
   TBranch        *b_rhod1lund;   //!
   TBranch        *b_rhod1idx;   //!
   TBranch        *b_rhod2lund;   //!
   TBranch        *b_rhod2idx;   //!
   TBranch        *b_npi;   //!
   TBranch        *b_pidocaxy_xy;   //!
   TBranch        *b_pidocaxy_xyerr;   //!
   TBranch        *b_pidocaxy_z;   //!
   TBranch        *b_pidocaxy_zerr;   //!
   TBranch        *b_picosth;   //!
   TBranch        *b_picosthcm;   //!
   TBranch        *b_pienergy;   //!
   TBranch        *b_pienergycm;   //!
   TBranch        *b_pip3;   //!
   TBranch        *b_pip3cm;   //!
   TBranch        *b_piphi;   //!
   TBranch        *b_piphicm;   //!
   TBranch        *b_pilund;   //!
   TBranch        *b_pitrkidx;   //!
   TBranch        *b_ngamma;   //!
   TBranch        *b_gammacosth;   //!
   TBranch        *b_gammacosthcm;   //!
   TBranch        *b_gammaenergy;   //!
   TBranch        *b_gammaenergycm;   //!
   TBranch        *b_gammap3;   //!
   TBranch        *b_gammap3cm;   //!
   TBranch        *b_gammaphi;   //!
   TBranch        *b_gammaphicm;   //!
   TBranch        *b_gammalund;   //!
   TBranch        *b_ntrk;   //!
   TBranch        *b_Trkdocaxy_xy;   //!
   TBranch        *b_Trkdocaxy_xyerr;   //!
   TBranch        *b_Trkdocaxy_z;   //!
   TBranch        *b_Trkdocaxy_zerr;   //!
   TBranch        *b_Trkcosth;   //!
   TBranch        *b_Trkcosthcm;   //!
   TBranch        *b_Trkenergy;   //!
   TBranch        *b_Trkenergycm;   //!
   TBranch        *b_Trkp3;   //!
   TBranch        *b_Trkp3cm;   //!
   TBranch        *b_Trkphi;   //!
   TBranch        *b_Trkphicm;   //!
   TBranch        *b_Trklund;   //!
   TBranch        *b_piselectorsmap;   //!
   TBranch        *b_eselectorsmap;   //!
   TBranch        *b_Kselectorsmap;   //!
   TBranch        *b_muselectorsmap;   //!

   yield(TTree *tree=0);
   virtual ~yield();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(std::string histFileName);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef yield_cxx
yield::yield(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("trees/run1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("trees/run1.root");
      }
      f->GetObject("h1",tree);

   }
   Init(tree);
}

yield::~yield()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t yield::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t yield::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void yield::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runnumber", &runnumber, &b_runnumber);
   fChain->SetBranchAddress("platform", &platform, &b_platform);
   fChain->SetBranchAddress("partition", &partition, &b_partition);
   fChain->SetBranchAddress("upperid", &upperid, &b_upperid);
   fChain->SetBranchAddress("lowerid", &lowerid, &b_lowerid);
   fChain->SetBranchAddress("majorid", &majorid, &b_majorid);
   fChain->SetBranchAddress("configkey", &configkey, &b_configkey);
   fChain->SetBranchAddress("date", &date, &b_date);
   fChain->SetBranchAddress("ddate", &ddate, &b_ddate);
   fChain->SetBranchAddress("eepx", &eepx, &b_eepx);
   fChain->SetBranchAddress("eepy", &eepy, &b_eepy);
   fChain->SetBranchAddress("eepz", &eepz, &b_eepz);
   fChain->SetBranchAddress("eee", &eee, &b_eee);
   fChain->SetBranchAddress("Bgfmultihadron", &Bgfmultihadron, &b_Bgfmultihadron);
   fChain->SetBranchAddress("Bgfneutralhadron", &Bgfneutralhadron, &b_Bgfneutralhadron);
   fChain->SetBranchAddress("Bgfmumu", &Bgfmumu, &b_Bgfmumu);
   fChain->SetBranchAddress("Bgftau", &Bgftau, &b_Bgftau);
   fChain->SetBranchAddress("Bgftwoprong", &Bgftwoprong, &b_Bgftwoprong);
   fChain->SetBranchAddress("Bgfphigamma", &Bgfphigamma, &b_Bgfphigamma);
   fChain->SetBranchAddress("Bgfallneutraltwophoton", &Bgfallneutraltwophoton, &b_Bgfallneutraltwophoton);
   fChain->SetBranchAddress("Bgfisr", &Bgfisr, &b_Bgfisr);
   fChain->SetBranchAddress("Bgfradtwoprong", &Bgfradtwoprong, &b_Bgfradtwoprong);
   fChain->SetBranchAddress("Bgfhighmasshadron", &Bgfhighmasshadron, &b_Bgfhighmasshadron);
   fChain->SetBranchAddress("Bgftwophotontwotrack", &Bgftwophotontwotrack, &b_Bgftwophotontwotrack);
   fChain->SetBranchAddress("Digifdchemcpreveto", &Digifdchemcpreveto, &b_Digifdchemcpreveto);
   fChain->SetBranchAddress("Digifl1open", &Digifl1open, &b_Digifl1open);
   fChain->SetBranchAddress("Digifl3open", &Digifl3open, &b_Digifl3open);
   fChain->SetBranchAddress("L3outdch", &L3outdch, &b_L3outdch);
   fChain->SetBranchAddress("L3outemc", &L3outemc, &b_L3outemc);
   fChain->SetBranchAddress("ntracks", &ntracks, &b_ntracks);
   fChain->SetBranchAddress("ngoodtrkloose", &ngoodtrkloose, &b_ngoodtrkloose);
   fChain->SetBranchAddress("mclen", &mclen, &b_mclen);
   fChain->SetBranchAddress("mclund", mclund, &b_mclund);
   fChain->SetBranchAddress("mothidx", mothidx, &b_mothidx);
   fChain->SetBranchAddress("daulen", daulen, &b_daulen);
   fChain->SetBranchAddress("dauidx", dauidx, &b_dauidx);
   fChain->SetBranchAddress("mcmass", mcmass, &b_mcmass);
   fChain->SetBranchAddress("mcp3cm", mcp3cm, &b_mcp3cm);
   fChain->SetBranchAddress("mccosthcm", mccosthcm, &b_mccosthcm);
   fChain->SetBranchAddress("mcphicm", mcphicm, &b_mcphicm);
   fChain->SetBranchAddress("mcenergycm", mcenergycm, &b_mcenergycm);
   fChain->SetBranchAddress("mcp3", mcp3, &b_mcp3);
   fChain->SetBranchAddress("mccosth", mccosth, &b_mccosth);
   fChain->SetBranchAddress("mcphi", mcphi, &b_mcphi);
   fChain->SetBranchAddress("mcenergy", mcenergy, &b_mcenergy);
   fChain->SetBranchAddress("npsi", &npsi, &b_npsi);
   fChain->SetBranchAddress("psimass", psimass, &b_psimass);
   fChain->SetBranchAddress("psimasserr", psimasserr, &b_psimasserr);
   fChain->SetBranchAddress("psicosth", psicosth, &b_psicosth);
   fChain->SetBranchAddress("psicosthcm", psicosthcm, &b_psicosthcm);
   fChain->SetBranchAddress("psienergy", psienergy, &b_psienergy);
   fChain->SetBranchAddress("psienergycm", psienergycm, &b_psienergycm);
   fChain->SetBranchAddress("psip3", psip3, &b_psip3);
   fChain->SetBranchAddress("psip3cm", psip3cm, &b_psip3cm);
   fChain->SetBranchAddress("psiphi", psiphi, &b_psiphi);
   fChain->SetBranchAddress("psiphicm", psiphicm, &b_psiphicm);
   fChain->SetBranchAddress("psilund", psilund, &b_psilund);
   fChain->SetBranchAddress("psid1lund", psid1lund, &b_psid1lund);
   fChain->SetBranchAddress("psid1idx", psid1idx, &b_psid1idx);
   fChain->SetBranchAddress("psid2lund", psid2lund, &b_psid2lund);
   fChain->SetBranchAddress("psid2idx", psid2idx, &b_psid2idx);
   fChain->SetBranchAddress("nel", &nel, &b_nel);
   fChain->SetBranchAddress("eldocaxy_xy", eldocaxy_xy, &b_eldocaxy_xy);
   fChain->SetBranchAddress("eldocaxy_xyerr", eldocaxy_xyerr, &b_eldocaxy_xyerr);
   fChain->SetBranchAddress("eldocaxy_z", eldocaxy_z, &b_eldocaxy_z);
   fChain->SetBranchAddress("eldocaxy_zerr", eldocaxy_zerr, &b_eldocaxy_zerr);
   fChain->SetBranchAddress("elcosth", elcosth, &b_elcosth);
   fChain->SetBranchAddress("elcosthcm", elcosthcm, &b_elcosthcm);
   fChain->SetBranchAddress("elenergy", elenergy, &b_elenergy);
   fChain->SetBranchAddress("elenergycm", elenergycm, &b_elenergycm);
   fChain->SetBranchAddress("elp3", elp3, &b_elp3);
   fChain->SetBranchAddress("elp3cm", elp3cm, &b_elp3cm);
   fChain->SetBranchAddress("elphi", elphi, &b_elphi);
   fChain->SetBranchAddress("elphicm", elphicm, &b_elphicm);
   fChain->SetBranchAddress("ellund", ellund, &b_ellund);
   fChain->SetBranchAddress("eld1lund", eld1lund, &b_eld1lund);
   fChain->SetBranchAddress("eld1idx", eld1idx, &b_eld1idx);
   fChain->SetBranchAddress("eld2lund", eld2lund, &b_eld2lund);
   fChain->SetBranchAddress("eld2idx", eld2idx, &b_eld2idx);
   fChain->SetBranchAddress("eld3lund", eld3lund, &b_eld3lund);
   fChain->SetBranchAddress("eld3idx", eld3idx, &b_eld3idx);
   fChain->SetBranchAddress("eld4lund", eld4lund, &b_eld4lund);
   fChain->SetBranchAddress("eld4idx", eld4idx, &b_eld4idx);
   fChain->SetBranchAddress("eltrkidx", eltrkidx, &b_eltrkidx);
   fChain->SetBranchAddress("netap", &netap, &b_netap);
   fChain->SetBranchAddress("etapmass", etapmass, &b_etapmass);
   fChain->SetBranchAddress("etapmasserr", etapmasserr, &b_etapmasserr);
   fChain->SetBranchAddress("etapcosth", etapcosth, &b_etapcosth);
   fChain->SetBranchAddress("etapcosthcm", etapcosthcm, &b_etapcosthcm);
   fChain->SetBranchAddress("etapenergy", etapenergy, &b_etapenergy);
   fChain->SetBranchAddress("etapenergycm", etapenergycm, &b_etapenergycm);
   fChain->SetBranchAddress("etapp3", etapp3, &b_etapp3);
   fChain->SetBranchAddress("etapp3cm", etapp3cm, &b_etapp3cm);
   fChain->SetBranchAddress("etapphi", etapphi, &b_etapphi);
   fChain->SetBranchAddress("etapphicm", etapphicm, &b_etapphicm);
   fChain->SetBranchAddress("etaplund", etaplund, &b_etaplund);
   fChain->SetBranchAddress("etapd1lund", etapd1lund, &b_etapd1lund);
   fChain->SetBranchAddress("etapd1idx", etapd1idx, &b_etapd1idx);
   fChain->SetBranchAddress("etapd2lund", etapd2lund, &b_etapd2lund);
   fChain->SetBranchAddress("etapd2idx", etapd2idx, &b_etapd2idx);
   fChain->SetBranchAddress("etapd3lund", etapd3lund, &b_etapd3lund);
   fChain->SetBranchAddress("etapd3idx", etapd3idx, &b_etapd3idx);
   fChain->SetBranchAddress("nrho", &nrho, &b_nrho);
   fChain->SetBranchAddress("rhomass", rhomass, &b_rhomass);
   fChain->SetBranchAddress("rhomasserr", rhomasserr, &b_rhomasserr);
   fChain->SetBranchAddress("rhocosth", rhocosth, &b_rhocosth);
   fChain->SetBranchAddress("rhocosthcm", rhocosthcm, &b_rhocosthcm);
   fChain->SetBranchAddress("rhoenergy", rhoenergy, &b_rhoenergy);
   fChain->SetBranchAddress("rhoenergycm", rhoenergycm, &b_rhoenergycm);
   fChain->SetBranchAddress("rhop3", rhop3, &b_rhop3);
   fChain->SetBranchAddress("rhop3cm", rhop3cm, &b_rhop3cm);
   fChain->SetBranchAddress("rhophi", rhophi, &b_rhophi);
   fChain->SetBranchAddress("rhophicm", rhophicm, &b_rhophicm);
   fChain->SetBranchAddress("rholund", rholund, &b_rholund);
   fChain->SetBranchAddress("rhod1lund", rhod1lund, &b_rhod1lund);
   fChain->SetBranchAddress("rhod1idx", rhod1idx, &b_rhod1idx);
   fChain->SetBranchAddress("rhod2lund", rhod2lund, &b_rhod2lund);
   fChain->SetBranchAddress("rhod2idx", rhod2idx, &b_rhod2idx);
   fChain->SetBranchAddress("npi", &npi, &b_npi);
   fChain->SetBranchAddress("pidocaxy_xy", pidocaxy_xy, &b_pidocaxy_xy);
   fChain->SetBranchAddress("pidocaxy_xyerr", pidocaxy_xyerr, &b_pidocaxy_xyerr);
   fChain->SetBranchAddress("pidocaxy_z", pidocaxy_z, &b_pidocaxy_z);
   fChain->SetBranchAddress("pidocaxy_zerr", pidocaxy_zerr, &b_pidocaxy_zerr);
   fChain->SetBranchAddress("picosth", picosth, &b_picosth);
   fChain->SetBranchAddress("picosthcm", picosthcm, &b_picosthcm);
   fChain->SetBranchAddress("pienergy", pienergy, &b_pienergy);
   fChain->SetBranchAddress("pienergycm", pienergycm, &b_pienergycm);
   fChain->SetBranchAddress("pip3", pip3, &b_pip3);
   fChain->SetBranchAddress("pip3cm", pip3cm, &b_pip3cm);
   fChain->SetBranchAddress("piphi", piphi, &b_piphi);
   fChain->SetBranchAddress("piphicm", piphicm, &b_piphicm);
   fChain->SetBranchAddress("pilund", pilund, &b_pilund);
   fChain->SetBranchAddress("pitrkidx", pitrkidx, &b_pitrkidx);
   fChain->SetBranchAddress("ngamma", &ngamma, &b_ngamma);
   fChain->SetBranchAddress("gammacosth", gammacosth, &b_gammacosth);
   fChain->SetBranchAddress("gammacosthcm", gammacosthcm, &b_gammacosthcm);
   fChain->SetBranchAddress("gammaenergy", gammaenergy, &b_gammaenergy);
   fChain->SetBranchAddress("gammaenergycm", gammaenergycm, &b_gammaenergycm);
   fChain->SetBranchAddress("gammap3", gammap3, &b_gammap3);
   fChain->SetBranchAddress("gammap3cm", gammap3cm, &b_gammap3cm);
   fChain->SetBranchAddress("gammaphi", gammaphi, &b_gammaphi);
   fChain->SetBranchAddress("gammaphicm", gammaphicm, &b_gammaphicm);
   fChain->SetBranchAddress("gammalund", gammalund, &b_gammalund);
   fChain->SetBranchAddress("ntrk", &ntrk, &b_ntrk);
   fChain->SetBranchAddress("Trkdocaxy_xy", Trkdocaxy_xy, &b_Trkdocaxy_xy);
   fChain->SetBranchAddress("Trkdocaxy_xyerr", Trkdocaxy_xyerr, &b_Trkdocaxy_xyerr);
   fChain->SetBranchAddress("Trkdocaxy_z", Trkdocaxy_z, &b_Trkdocaxy_z);
   fChain->SetBranchAddress("Trkdocaxy_zerr", Trkdocaxy_zerr, &b_Trkdocaxy_zerr);
   fChain->SetBranchAddress("Trkcosth", Trkcosth, &b_Trkcosth);
   fChain->SetBranchAddress("Trkcosthcm", Trkcosthcm, &b_Trkcosthcm);
   fChain->SetBranchAddress("Trkenergy", Trkenergy, &b_Trkenergy);
   fChain->SetBranchAddress("Trkenergycm", Trkenergycm, &b_Trkenergycm);
   fChain->SetBranchAddress("Trkp3", Trkp3, &b_Trkp3);
   fChain->SetBranchAddress("Trkp3cm", Trkp3cm, &b_Trkp3cm);
   fChain->SetBranchAddress("Trkphi", Trkphi, &b_Trkphi);
   fChain->SetBranchAddress("Trkphicm", Trkphicm, &b_Trkphicm);
   fChain->SetBranchAddress("Trklund", Trklund, &b_Trklund);
   fChain->SetBranchAddress("piselectorsmap", piselectorsmap, &b_piselectorsmap);
   fChain->SetBranchAddress("eselectorsmap", eselectorsmap, &b_eselectorsmap);
   fChain->SetBranchAddress("Kselectorsmap", Kselectorsmap, &b_Kselectorsmap);
   fChain->SetBranchAddress("muselectorsmap", muselectorsmap, &b_muselectorsmap);
   Notify();
}

Bool_t yield::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void yield::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t yield::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef yield_cxx
