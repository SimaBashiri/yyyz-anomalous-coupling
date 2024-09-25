#include "analysis.h"
#include "TSystem.h"

int main(){

    TString path = "/eos/cms/store/group/phys_pps/Phase2/Delphes/PU200/";
    int tp = 20;

    TChain* ch = new TChain("Delphes");
    analysis t(ch);

    // zz1
    ch->Add(path + "ZToMuMu_M-120to200_TuneCP5_14TeV-powheg-pythia8/*");
    t.Loop(18.72e-12, "all", "V", "ztomumu1", tp);

    ch->Reset();  // Reset the TChain to clear the previous files
    ch->Add(path + "ZToMuMu_M-120to200_TuneCP5_14TeV-powheg-pythia8/*");
    t.Loop(18.72e-12, "all", "H", "ztomumu1", tp);

    // zz2
    ch->Reset();
    ch->Add(path + "ZToMuMu_M-200to400_TuneCP5_14TeV-powheg-pythia8/*");
    t.Loop(2.682e-12, "all", "V", "ztomumu2", tp);

    ch->Reset();
    ch->Add(path + "ZToMuMu_M-200to400_TuneCP5_14TeV-powheg-pythia8/*");
    t.Loop(2.682e-12, "all", "H", "ztomumu2", tp);

    // zz3
    ch->Reset();
    ch->Add(path + "ZToMuMu_M-400to800_TuneCP5_14TeV-powheg-pythia8/*");
    t.Loop(0.2396e-12, "all", "V", "ztomumu3", tp);

    ch->Reset();
    ch->Add(path + "ZToMuMu_M-400to800_TuneCP5_14TeV-powheg-pythia8/*");
    t.Loop(0.2396e-12, "all", "H", "ztomumu3", tp);

    // // mad
    // ch->Reset();
    // ch->Add(path + "MadGraph5_14TeV_pptoZA-noHadr-1.314E2pm7.0E-1_Zmumu_Delphes_PU200_v2.root");
    // t.Loop(1.314e-10, "all", "H", "mad", tp);

    // ch->Reset();
    // ch->Add(path + "MadGraph5_14TeV_pptoZA-noHadr-1.314E2pm7.0E-1_Zmumu_Delphes_PU200_v2.root");
    // t.Loop(1.314e-10, "all", "V", "mad", tp);

    // zy
    ch->Reset();
    ch->Add(path + "Zgamma_inc_SM_Madgraph5_Delphes_PU200.root");
    t.Loop(1, "all", "H", "zy", tp);

    ch->Reset();
    ch->Add(path + "Zgamma_inc_SM_Madgraph5_Delphes_PU200.root");
    t.Loop(1, "all", "V", "zy", tp);

    // signal1
    ch->Reset();
    ch->Add(path + "/FPMC_bSM_14tev_AAAZeft_A1A_0E0_A2A_1E-13_pt50_horXing-noHADR_2.4390029090069734E-003_Zmumu_Delphes_PU200.root");
    t.Loop(2.439E-15, "all", "H", "signal1", tp);

    ch->Reset();
    ch->Add(path + "/FPMC_bSM_14tev_AAAZeft_A1A_0E0_A2A_1E-13_pt50_horXing-noHADR_2.4390029090069734E-003_Zmumu_Delphes_PU200.root");
    t.Loop(2.439E-15, "nopu", "H", "signal1", tp);

    ch->Reset();
    ch->Add(path + "/FPMC_bSM_14tev_AAAZeft_A1A_0E0_A2A_1E-13_pt50-noHADR_3.5563044069199587E-004_Zmumu.root");
    t.Loop(3.55E-16, "nopu", "V", "signal1", tp);

    ch->Reset();
    ch->Add(path + "/FPMC_bSM_14tev_AAAZeft_A1A_0E0_A2A_1E-13_pt50-noHADR_3.5563044069199587E-004_Zmumu.root");
    t.Loop(3.55E-16, "all", "V", "signal1", tp);

    // zjet
    ch->Reset();
    ch->Add(path + "ZJets_inc_SM_Madgraph5_JetPT200GeV_Delphes_PU200/*.root");
    t.Loop(60.517398e-12, "all", "V", "zjet", tp);

    ch->Reset();
    ch->Add(path + "ZJets_inc_SM_Madgraph5_JetPT200GeV_Delphes_PU200/*.root");
    t.Loop(60.517398e-12, "all", "H", "zjet", tp);

    return 0;
}
