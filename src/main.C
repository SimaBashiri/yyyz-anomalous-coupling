#include "analysis.h"
#include "TSystem.h"

int main(){

    TString path = "/eos/cms/store/group/phys_pps/Phase2/Delphes/PU200/";
    int tp = 20;

    //zz1
    TChain* ch = new TChain("Delphes") ;
    ch->Add(path + "ZToMuMu_M-120to200_TuneCP5_14TeV-powheg-pythia8/*");
    analysis t(ch);
    t.Loop(18.72e-12, "all", "V", "ztomumu1", tp);

    TChain* ch = new TChain("Delphes") ;
    ch->Add(path + "ZToMuMu_M-120to200_TuneCP5_14TeV-powheg-pythia8/*");
    analysis t(ch);
    t.Loop(18.72e-12, "all", "H", "ztomumu1", tp);

    //zz2
    TChain* ch = new TChain("Delphes") ;
    ch->Add(path + "ZToMuMu_M-200to400_TuneCP5_14TeV-powheg-pythia8/*");
    analysis t(ch);
    t.Loop(2.682e-12, "all", "V", "ztomumu2", tp);

    TChain* ch = new TChain("Delphes") ;
    ch->Add(path + "ZToMuMu_M-200to400_TuneCP5_14TeV-powheg-pythia8/*");
    analysis t(ch);
    t.Loop(2.682e-12, "all", "H", "ztomumu2", tp);

    //zz3
    TChain* ch = new TChain("Delphes") ;
    ch->Add(path + "ZToMuMu_M-400to800_TuneCP5_14TeV-powheg-pythia8/*");
    analysis t(ch);
    t.Loop(0.2396e-12, "all", "V", "ztomumu3", tp);

    TChain* ch = new TChain("Delphes") ;
    ch->Add(path + "ZToMuMu_M-400to800_TuneCP5_14TeV-powheg-pythia8/*");
    analysis t(ch);
    t.Loop(0.2396e-12, "all", "H", "ztomumu3", tp);

    //mad
    TChain* ch = new TChain("Delphes") ;
    ch->Add(path + "MadGraph5_14TeV_pptoZA-noHadr-1.314E2pm7.0E-1_Zmumu_Delphes_PU200_v2.root");
    analysis t(ch);
    t.Loop(1.314e-10, "all", "H", "mad", tp);

    TChain* ch = new TChain("Delphes") ;
    ch->Add(path + "MadGraph5_14TeV_pptoZA-noHadr-1.314E2pm7.0E-1_Zmumu_Delphes_PU200_v2.root");
    analysis t(ch);
    t.Loop(1.314e-10, "all", "V", "mad", tp);


    //zy
    TChain* ch = new TChain("Delphes") ;
    ch->Add(path + "Zgamma_inc_SM_Madgraph5_Delphes_PU200.root");
    analysis t(ch);
    t.Loop(1, "all", "H", "zy", tp);

    TChain* ch = new TChain("Delphes") ;
    ch->Add(path + "Zgamma_inc_SM_Madgraph5_Delphes_PU200.root");
    analysis t(ch);
    t.Loop(1, "all", "V", "zy", tp);


    //signal
    TChain* ch = new TChain("Delphes") ;
    ch->Add(path + "/FPMC_bSM_14tev_AAAZeft_A1A_0E0_A2A_1E-13_pt50_horXing-noHADR_2.4390029090069734E-003_Zmumu_Delphes_PU200.root");
    analysis t(ch);
    t.Loop(2.439E-15, "all", "H", "signal1", tp);

    TChain* ch = new TChain("Delphes") ;
    ch->Add(path + "/FPMC_bSM_14tev_AAAZeft_A1A_0E0_A2A_1E-13_pt50_horXing-noHADR_2.4390029090069734E-003_Zmumu_Delphes_PU200.root");
    analysis t(ch);
    t.Loop(2.439E-15, "nopu", "H", "signal1", tp);

    TChain* ch = new TChain("Delphes") ;
    ch->Add(path + "/FPMC_bSM_14tev_AAAZeft_A1A_0E0_A2A_1E-13_pt50-noHADR_3.5563044069199587E-004_Zmumu.root");
    analysis t(ch);
    t.Loop(3.55E-16, "nopu", "V", "signal1", tp);


    TChain* ch = new TChain("Delphes") ;
    ch->Add(path + "/FPMC_bSM_14tev_AAAZeft_A1A_0E0_A2A_1E-13_pt50-noHADR_3.5563044069199587E-004_Zmumu.root");
    analysis t(ch);
    t.Loop(3.55E-16, "all", "V", "signal1", tp);

    return 0;
}
