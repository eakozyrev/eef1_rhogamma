{

 cout << "d******************************" << endl;
 cout << "*  Welcome Evgeny to ROOT v" << gROOT->GetVersion() << "  *" << endl;
 cout << "******************************" << endl;
 cout << endl;
 // gROOT->ProcessLine(".include /usr/local/include/root");
 gSystem->AddIncludePath(" -I/usr/local/include/root");
 gStyle->SetCanvasColor(kWhite);     // background is no longer mouse-dropping white
 gStyle->SetPalette(1,0);            // blue to red false color palette. Use 9 for b/w
 gStyle->SetCanvasBorderMode(0);     // turn off canvas borders
 gStyle->SetPadBorderMode(0);
 gStyle->SetPaintTextFormat("5.2f");  // What precision to put numbers if plotted with "TEXT"

 gStyle->SetHistLineWidth(1.1);
 gStyle->SetMarkerSize(1.3);
 gStyle->SetLegendBorderSize(0);

 // For publishing:
 gStyle->SetLineWidth(3.);
 gStyle->SetTextSize(1.8);
 gStyle->SetLabelSize(0.06,"xy");
 gStyle->SetTitleSize(0.07,"xy");
 gStyle->SetTitleOffset(0.9,"x");
 gStyle->SetTitleOffset(0.7,"y");
 gStyle->SetPadTopMargin(0.1);
 gStyle->SetPadRightMargin(0.1);
 gStyle->SetPadBottomMargin(0.16);
 gStyle->SetPadLeftMargin(0.12);
// gStyle->SetOptStat(0);

 

 
//gROOT->LoadMacro("/home/eakozyrev/cmd3/analysis/rootplay/4pi/selection/results/./../../mygenerator/Cmd3Generator2pi2pi0_ke_cc.so");
 //gROOT->LoadMacro("../../mygenerator/Cmd3Generator2pi2pi0_ke.cc++");
 //gROOT->LoadMacro("alpha_PWA.C+");
 //Cmd3Generator2pi2pi0_ke t;
 //gROOT->LoadMacro("../mygenerator/popov/TFindEpsilon.cc+");
 //gROOT->LoadMacro("../mygenerator/popov/TPhaseSpace.cc+");
 //gROOT->LoadMacro("../cross_sections/cross_fit4pi.C+");
 //gROOT->LoadMacro("../cross_sections/plot3pi.C+");
 //gROOT->LoadMacro("../cross_sections/cross_fit.C+");
 //gROOT->LoadMacro("alpha.C+");
}
