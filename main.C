#include "analysis.h"
#include "TSystem.h"

int main(){
    TString path = "/home/bashiri/Downloads/";
    int tp = 20;
    // FPMC
    gROOT->ProcessLine(".L analysis.cc");
    TChain* ch = new TChain("Delphes") ;
    ch->Add(path + "/FPMC_bSM_14tev_AAAZeft_A1A_0E0_A2A_1E-13_pt50-noHADR_3.5563044069199587E-004_Zmumu.root");
    analysis t(ch);
    t.Loop(3.55E-16, "all", "V", "signal1", tp);
    return 0;
}
