#include "analysis.h"
#include "TSystem.h"

int main() {
    TString path = "/eos/cms/store/group/phys_pps/Phase2/Delphes/PU200/";
    int tp = 20;

    TChain* ch = new TChain("Delphes");
    analysis t(ch);

    // // zy1
    // ch->Add(path + "ZToMuMu_M-120to200_TuneCP5_14TeV-powheg-pythia8/*");
    // t.Loop(18.72e-12, "all", "V", "ztomumu1", tp);
    // t.Loop(18.72e-12, "all", "H", "ztomumu1", tp);

    // // zy2
    // ch->SetEntries(0);
    // ch->Add(path + "ZToMuMu_M-200to400_TuneCP5_14TeV-powheg-pythia8/*");
    // t.Loop(2.682e-12, "all", "V", "ztomumu2", tp);
    // t.Loop(2.682e-12, "all", "H", "ztomumu2", tp);

    // // zy3
    // ch->SetEntries(0);
    // ch->Add(path + "ZToMuMu_M-400to800_TuneCP5_14TeV-powheg-pythia8/*");
    // t.Loop(0.2396e-12, "all", "V", "ztomumu3", tp);
    // t.Loop(0.2396e-12, "all", "H", "ztomumu3", tp);

    // // mad
    // ch->SetEntries(0);
    // ch->Add(path + "MadGraph5_14TeV_pptoZA-noHadr-1.314E2pm7.0E-1_Zmumu_Delphes_PU200_v2.root");
    // t.Loop(1.314e-10, "all", "H", "mad", tp);
    
    // ch->SetEntries(0);
    // ch->Add(path + "MadGraph5_14TeV_pptoZA-noHadr-1.314E2pm7.0E-1_Zmumu_Delphes_PU200_v2.root");
    // t.Loop(1.314e-10, "all", "V", "mad", tp);

    // zy SM
    ch->SetEntries(0);
    ch->Add(path + "/Zgamma_inc_SM_Madgraph5_PhotonPT200GeV_Delphes_PU200/*.root");
    t.Loop(0.152e-12, "all", "H", "zy", tp);
    t.Loop(0.152e-12, "all", "V", "zy", tp);

    // signal1
    ch->SetEntries(0);
    ch->Add(path + "../horXing_PU200/FPMC_bSM_14tev_AAAZeft_A1A_0E0_A2A_1E-13_pt50_horXing-noHADR_2.4390029090069734E-003_Zmumu_Delphes_PU200.root");
    t.Loop(2.439E-15, "all", "H", "signal1", tp);

    ch->SetEntries(0);
    ch->Add(path + "../horXing_PU200/FPMC_bSM_14tev_AAAZeft_A1A_0E0_A2A_1E-13_pt50_horXing-noHADR_2.4390029090069734E-003_Zmumu_Delphes_PU200.root");
    t.Loop(2.439E-15, "nopu", "H", "signal1", tp);

    ch->SetEntries(0);
    ch->Add(path + "/FPMC_bSM_14tev_AAAZeft_A1A_0E0_A2A_1E-13_pt50-noHADR_3.5563044069199587E-004_Zmumu.root");
    t.Loop(3.55E-16, "nopu", "V", "signal1", tp);

    ch->SetEntries(0);
    ch->Add(path + "/FPMC_bSM_14tev_AAAZeft_A1A_0E0_A2A_1E-13_pt50-noHADR_3.5563044069199587E-004_Zmumu.root");
    t.Loop(3.55E-16, "all", "V", "signal1", tp);

    // zjet
    ch->SetEntries(0);
    ch->Add(path + "ZJets_inc_SM_Madgraph5_JetPT200GeV_Delphes_PU200/*.root");
    t.Loop(60.517398e-12, "all", "V", "zjet", tp);

    ch->SetEntries(0);
    ch->Add(path + "ZJets_inc_SM_Madgraph5_JetPT200GeV_Delphes_PU200/*.root");
    t.Loop(60.517398e-12, "all", "H", "zjet", tp);

    delete ch;
    return 0;
}
