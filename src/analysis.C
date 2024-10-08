#define analysis_cxx
#include "analysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TRandom3.h>
#include <TLatex.h>
#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <TSystem.h>
#include <utility>



using namespace std;
void displayProgress(long current, long max){
    using std::cerr;
    if (max<2500) return;
    if (current%(max/2500)!=0 && current<max-1) return;

    int width = 52; // Hope the terminal is at least that wide.
    int barWidth = width - 2;
    cerr << "\x1B[2K"; // Clear line
    cerr << "\x1B[2000D"; // Cursor left
    cerr << '[';
    for(int i=0 ; i<barWidth ; ++i){ if(i<barWidth*current/max){ cerr << '=' ; }  else{ cerr << ' ' ; } }
    cerr << ']';
    cerr << " " << Form("%8d/%8d (%5.2f%%)", (int)current, (int)max, 100.0*current/max) ;
    cerr.flush();
}

void cleanupVectors(std::vector<lepton_candidate*>* leptons, std::vector<photon_candidate*>* photons, std::vector<proton_candidate*>* protons) {
    for (auto& lepton : *leptons) {
        delete lepton;
    }
    for (auto& photon : *photons) {
        delete photon;
    }
    for (auto& proton : *protons) {
        delete proton;
    }
    leptons->clear();
    photons->clear();
    protons->clear();
    delete leptons;
    delete photons;
    delete protons;
}

bool ComparePtLep(lepton_candidate *a, lepton_candidate *b) { return a->pt_ > b->pt_; }
bool ComparePtPhoton(photon_candidate *a, photon_candidate *b) { return a->pt_ > b->pt_; }
bool ComparePzProton(proton_candidate *a, proton_candidate *b) { return abs(a->p4_.Pz()) < abs(b->p4_.Pz()); }
bool CompareChargeProton(proton_candidate *a, proton_candidate *b) { return TMath::Sign(1, a->p4_.Pz()) > TMath::Sign(1, b->p4_.Pz()); }

Double_t deltaPhi(Double_t phi1, Double_t phi2) {
    Double_t dPhi = phi1 - phi2;
    if (dPhi > TMath::Pi()) dPhi -= 2.*TMath::Pi();
    if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();
    return dPhi;
}


Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2) {
    Double_t dEta, dPhi ;
    dEta = eta1 - eta2;
    dPhi = deltaPhi(phi1, phi2);
    return sqrt(dEta*dEta+dPhi*dPhi);
}

Double_t Rapidity(Double_t EE, Double_t ZZ){
    return 0.5*TMath::Log((EE + ZZ)/(EE - ZZ));
}


std::pair<float, float> zt_Vertex_(std::vector<proton_candidate*> *selectedProtons, 
                  Float_t *GenProton_T, 
                  Float_t *GenProton_Z, 
                  int selp, int selm) {
    float z_Vertex = 0;  
    float t_Vertex = 0;

    if (selectedProtons->size() > 1) { 
        float tp1, tp2;
        float C = 30;  // cm/nsec

        if ((*selectedProtons)[selp]->p4_.Pz() > 0) {
          if( ((*selectedProtons)[selp]->p4_).Pz() > 0 )
              tp1 = (1e9 * GenProton_T[(*selectedProtons)[selp]->indice_] + (23400 - 0.1 * GenProton_Z[(*selectedProtons)[selp]->indice_]) / 30) + (*selectedProtons)[selp]->smearfac_;
          else
              tp1 = (1e9 * GenProton_T[(*selectedProtons)[selp]->indice_] + (23400 + 0.1 * GenProton_Z[(*selectedProtons)[selp]->indice_]) / 30) + (*selectedProtons)[selp]->smearfac_;

          if( ((*selectedProtons)[selm]->p4_).Pz() > 0 )
              tp2 = (1e9 * GenProton_T[(*selectedProtons)[selm]->indice_] + (23400 - 0.1 * GenProton_Z[(*selectedProtons)[selm]->indice_]) / 30) + (*selectedProtons)[selm]->smearfac_;
          else
              tp2 = (1e9 * GenProton_T[(*selectedProtons)[selm]->indice_] + (23400 + 0.1 * GenProton_Z[(*selectedProtons)[selm]->indice_]) / 30) + (*selectedProtons)[selm]->smearfac_;

          z_Vertex = (tp2 - tp1) * C / 2;
          t_Vertex = ((tp1 + tp2) / 2 - 23400 / C);
        }
    }
    return std::make_pair(z_Vertex, t_Vertex);
}


void histogram(TH1F *Hists[2][8][32+2+2], TH2F *Hists2[2][8][7], int ch, int cn, 
    std::vector<lepton_candidate*> *selectedLeptons, 
    std::vector<photon_candidate*> *selectedPhotons, 
    std::vector<proton_candidate*> *selectedProtons, 
    float weight, Float_t *GenProton_Rapidity, 
    Float_t *GenProton_T, Float_t *GenProton_Z, 
    float *Vertex_T, float *Vertex_Z, 
    Int_t *GenProton_IsPU, float timepc, 
    int selp, int selm) {


    double YZ_Rapidity;
    float xi_cms1;
    float xi_cms2;
    float MX;
    float YX;
    float sum_p;
    double sum_pz;
    float pp_rapidity;
    float t_Vertex;
    float z_Vertex;
    float diff1;
    float diff2;

    if((*selectedPhotons).size() > 0 && (selectedLeptons->size() > 1)){

      xi_cms1 = ( ((*selectedLeptons)[0]->p4_).E() + ((*selectedLeptons)[1]->p4_).E() + ((*selectedPhotons)[0]->p4_).E() + ((*selectedLeptons)[0]->p4_).Pz() + ((*selectedLeptons)[1]->p4_).Pz() + ((*selectedPhotons)[0]->p4_).Pz() )/14000;
      xi_cms2 = ( ((*selectedLeptons)[0]->p4_).E() + ((*selectedLeptons)[1]->p4_).E() + ((*selectedPhotons)[0]->p4_).E() - (((*selectedLeptons)[0]->p4_).Pz() + ((*selectedLeptons)[1]->p4_).Pz() + ((*selectedPhotons)[0]->p4_).Pz()) )/14000;
   
      MX = 14000 * TMath::Sqrt(xi_cms1 * xi_cms2);
      YX = 0.5 * TMath::Log(xi_cms1 / xi_cms2);
      YZ_Rapidity = Rapidity( ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_ + (*selectedPhotons)[0]->p4_).E(),  ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_ + (*selectedPhotons)[0]->p4_).Pz() );
      sum_p = TMath::Sqrt( ((*selectedLeptons)[0]->p4_).Px()*((*selectedLeptons)[0]->p4_).Px() + ((*selectedLeptons)[0]->p4_).Py()*((*selectedLeptons)[0]->p4_).Py() +
      ((*selectedLeptons)[0]->p4_).Pz()*((*selectedLeptons)[0]->p4_).Pz() ) + TMath::Sqrt( ((*selectedLeptons)[1]->p4_).Px()*((*selectedLeptons)[1]->p4_).Px() + ((*selectedLeptons)[1]->p4_).Py()*((*selectedLeptons)[1]->p4_).Py() +
      ((*selectedLeptons)[1]->p4_).Pz()*((*selectedLeptons)[1]->p4_).Pz() ) + TMath::Sqrt( ((*selectedPhotons)[0]->p4_).Px()*((*selectedPhotons)[0]->p4_).Px() + ((*selectedPhotons)[0]->p4_).Py()*((*selectedPhotons)[0]->p4_).Py() +
      ((*selectedPhotons)[0]->p4_).Pz()*((*selectedPhotons)[0]->p4_).Pz() );

      sum_pz = ((*selectedLeptons)[0]->p4_).Pz() + ((*selectedLeptons)[1]->p4_).Pz() + ((*selectedPhotons)[0]->p4_).Pz();
   }

    if(selectedLeptons->size() > 1){
        Hists[ch][cn][0]->Fill((*selectedLeptons)[0]->pt_,weight);
        Hists[ch][cn][1]->Fill((*selectedLeptons)[0]->eta_,weight);
        Hists[ch][cn][2]->Fill((*selectedLeptons)[0]->phi_,weight);

        Hists[ch][cn][3]->Fill((*selectedLeptons)[1]->pt_,weight);
        Hists[ch][cn][4]->Fill((*selectedLeptons)[1]->eta_,weight);
        Hists[ch][cn][5]->Fill((*selectedLeptons)[1]->phi_,weight);
    }
    if(selectedPhotons->size() > 0){
      Hists[ch][cn][6]->Fill((*selectedPhotons)[0]->pt_,weight);
      Hists[ch][cn][7]->Fill((*selectedPhotons)[0]->eta_,weight);
      Hists[ch][cn][8]->Fill((*selectedPhotons)[0]->phi_,weight);
    }
    if(selectedLeptons->size() > 1){ 
        Hists[ch][cn][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M() ,weight);
        Hists[ch][cn][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt() ,weight);
        Hists[ch][cn][11]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_), weight);
        Hists[ch][cn][12]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_), weight);
    }
    if(selectedProtons->size() > 1){
      pp_rapidity = 0.5*TMath::Log((((*selectedProtons)[0]->p4_ + (*selectedProtons)[1]->p4_).E() + ((*selectedProtons)[0]->p4_ + (*selectedProtons)[1]->p4_).Pz())/(((*selectedProtons)[0]->p4_ + (*selectedProtons)[1]->p4_).E() - ((*selectedProtons)[0]->p4_ + (*selectedProtons)[1]->p4_).Pz()));

      float C = 30;  //cm/nsec
      float tp1, tp2;
      float tr = 1E9 * timepc * 1E-12; //nsec


      z_Vertex = zt_Vertex_(selectedProtons, GenProton_T, GenProton_Z, selp, selm).first;
      t_Vertex = zt_Vertex_(selectedProtons, GenProton_T, GenProton_Z, selp, selm).second;
      // cout << "zv1: " << z_Vertex << "  " << "zv2: " << Vertex_Z[0]*0.1  << "   minzV: " << (Vertex_Z[0]*0.1) - z_Vertex << endl;
//         float pp_eta = ((*selectedProtons)[0]->p4_ + (*selectedProtons)[1]->p4_).Eta();
      Hists[ch][cn][13]->Fill((*selectedProtons)[0]->pt_,weight);
      Hists[ch][cn][14]->Fill(GenProton_Rapidity[(*selectedProtons)[0]->indice_],weight);
      Hists[ch][cn][15]->Fill((*selectedProtons)[0]->phi_,weight);
      Hists[ch][cn][16]->Fill((*selectedProtons)[1]->pt_,weight);
      Hists[ch][cn][17]->Fill(GenProton_Rapidity[(*selectedProtons)[1]->indice_],weight);
      Hists[ch][cn][18]->Fill((*selectedProtons)[1]->phi_,weight);
    }
    if((*selectedProtons).size() > 0)   diff1 = (*selectedProtons)[0]->xi_-((sum_p + sum_pz)/14000);
    if((*selectedProtons).size() > 1)   diff2 = (*selectedProtons)[1]->xi_-((sum_p - sum_pz)/14000);
          
      
    if(selectedPhotons->size() > 0 && selectedLeptons->size() > 1){
      Hists[ch][cn][19]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_ + (*selectedPhotons)[0]->p4_).M(),weight);
      Hists[ch][cn][20]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_ + (*selectedPhotons)[0]->p4_).Pt(),weight);
    }
    if((*selectedProtons).size() > 1 && (selectedPhotons->size() > 0))  Hists[ch][cn][21]->Fill(YZ_Rapidity,weight);
    if((*selectedProtons).size() > 1 && (selectedPhotons->size() > 0))  Hists[ch][cn][22]->Fill( YX - YZ_Rapidity ,weight);
    if(selectedPhotons->size() > 0)  Hists[ch][cn][23]->Fill(YX ,weight);
    if((*selectedProtons).size() > 0 && (selectedPhotons->size() > 0))  Hists[ch][cn][24]->Fill((*selectedProtons)[0]->xi_ - xi_cms1 ,weight);
    if((*selectedProtons).size() > 1 && (selectedPhotons->size() > 0))  Hists[ch][cn][24]->Fill((*selectedProtons)[1]->xi_ - xi_cms2 ,weight);
    if((*selectedProtons).size() > 0 && (selectedPhotons->size() > 0))  Hists[ch][cn][25]->Fill(diff1,weight);
    if((*selectedProtons).size() > 1 && (selectedPhotons->size() > 0))  Hists[ch][cn][26]->Fill(diff2,weight);
    if((*selectedProtons).size() > 1)  Hists[ch][cn][27]->Fill( Vertex_T[0]*1e9 - t_Vertex ,weight);
    if((*selectedProtons).size() > 1)  Hists[ch][cn][28]->Fill( (Vertex_Z[0]*0.1) - z_Vertex ,weight);
    if((*selectedPhotons).size() > 0 && selectedLeptons->size() > 1)  Hists2[ch][cn][0]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_ + (*selectedPhotons)[0]->p4_).M(), MX, weight);
    if((*selectedProtons).size() > 1 && selectedLeptons->size() > 1)  Hists2[ch][cn][1]->Fill( YZ_Rapidity , YX, weight);
  //     Hists2[ch][3][0][2]->Fill(YZ_Rapidity , pp_rapidity);
  //     Hists2[ch][3][0][3]->Fill( pp_rapidity , YX);
    if((*selectedPhotons).size() > 0 && (selectedLeptons->size() > 1))  Hists2[ch][cn][2]->Fill(xi_cms1 , xi_cms2, weight);
    if((*selectedProtons).size() > 1)  Hists2[ch][cn][3]->Fill(Vertex_T[0]*1e9 , t_Vertex, weight);
    if( (*selectedPhotons).size() > 0 && (selectedLeptons->size() > 1))  Hists2[ch][cn][4]->Fill( xi_cms1, YZ_Rapidity, weight);
    if((*selectedProtons).size() > 0 && (*selectedPhotons).size() > 0 && (selectedLeptons->size() > 1))  Hists2[ch][cn][5]->Fill( xi_cms1, (*selectedProtons)[0]->xi_, weight);
    if((*selectedProtons).size() > 1 && (*selectedPhotons).size() > 0 && (selectedLeptons->size() > 1))  Hists2[ch][cn][5]->Fill( xi_cms2, (*selectedProtons)[1]->xi_ , weight);
    if((*selectedProtons).size() > 1)  Hists2[ch][cn][6]->Fill(  Vertex_Z[0]*0.1,  z_Vertex, weight);
    if((*selectedProtons).size() > 1)  Hists[ch][cn][29]->Fill(GenProton_IsPU[(*selectedProtons)[0]->indice_] , weight);
    if((*selectedProtons).size() > 1)  Hists[ch][cn][30]->Fill(GenProton_IsPU[(*selectedProtons)[1]->indice_] , weight);
    if((*selectedPhotons).size() > 0 && (selectedLeptons->size() > 1)){
        float ZY_dPhi = deltaPhi(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Phi(), (*selectedPhotons)[0]->phi_);
        Hists[ch][cn][31]->Fill(ZY_dPhi , weight);

    }


}


void analysis::Loop(double cross_section, TString puflag, TString xiflag, TString ssig, int timepc){
    int ch = -1;
    double weight = 1;
    Long64_t nAccept=0;
    float smear = 0;
    float smear_p = 0;
    float smear_p2 = 0;
    double lumi = 3000e15;
    int protonSizeCut = 0;
    int ele_protonSizeCut = 0;
    int mu_protonSizeCut = 0;
    int channelCut = 0;
    int mu_channelCut = 0;
    int photonSizeCut = 0;
    int ele_photonSizeCut = 0;
    int mu_photonSizeCut = 0;
    int leptonSizeCut = 0;
    int ele_leptonSizeCut = 0;
    int mu_leptonSizeCut = 0;
    int cms_xiCut = 0;
    int ele_r1cms_xiCut = 0;
    int ele_r2cms_xiCut = 0;
    int mu_r1cms_xiCut = 0;
    int mu_r2cms_xiCut = 0;
    int ele_xiCut = 0;
    int mu_xiCut = 0;
    int ZVertexCut = 0;
    int ele_ZVertexCut = 0;
    int mu_ZVertexCut = 0;
    int timingCut = 0;
    int ele_timingCut = 0;
    int mu_timingCut = 0;
    int MassCut = 0;
    int ele_MassCut = 0;
    int mu_MassCut = 0;
    int XiResulutionCut = 0;
    int ele_r1XiResulutionCut = 0;
    int mu_r2XiResulutionCut = 0;
    int ele_r2XiResulutionCut = 0;
    int mu_r1XiResulutionCut = 0;
    int RapidityCut = 0;
    int ele_RapidityCut = 0;
    int mu_RapidityCut = 0;
    int Mzwindow = 0;
    int ele_Mzwindow = 0;
    int mu_Mzwindow = 0;
    int xicondition = 0;
    int ele_xicondition = 0;
    int mu_xicondition = 0;
    int PPSXiCut = 0;
    int ele_PPSXiCut = 0;
    int mu_PPSXiCut = 0;
    int ele_zyDeltaPhi = 0;
    int mu_zyDeltaPhi = 0;
    int zyDeltaPhiCut = 0;
    TString spu;
    TString sxi;
    float xi_min;
    float xi_max;
    int region = 0;
    int ZPt_cut = 0;
    int SampleCount = 0;
    float xi_cms1;
    float xi_cms2;
    float z_Vertex;
    float t_Vertex;
    int selp = -1;
    int selm = -1;
    float smear_p1;
    float C = 30;  //cm/nsec
    float tp1, tp2;
    TRandom3 *r1 = new TRandom3();
    r1->SetSeed(1000000);
    float tr = 1E9* timepc * 1E-12; //nsec
    float pz_min = (1-xi_max)*7000;
    float pz_max = (1-xi_min)*7000;
    float xi_smear;
    int indP_ = -1;
    int indM_ = -1;
    float xiPPS;
    float z_vertex_result;
    float minZVertex;

    if(puflag == "nopu")  spu = "noPU";
    else     spu = "PU";

    if(xiflag == "V"){
        //Vertical
      xi_min = 0.0147; // using the 234m station
      xi_max = 0.196; // using the 196m station
      sxi = "Vertical";
    }
    if(xiflag == "H"){
      //Horizental
        xi_min = 0.0472; // using the 234m station
        xi_max = 0.287; // using the 196m station
        sxi = "Horizontal" ;
    }
    TString path = "./";
    std::vector<TString> regions{"ZPt", "ny", "Mzwindow", "zyDeltaPhi", "nProton", "XiResolutionCut", "ZVertexCut", "timingCut"};
    std::vector<TString> channels{"ee", "mumu"};
    std::vector<TString> vars   {"lep1Pt","lep1Eta","lep1Phi","lep2Pt","lep2Eta","lep2Phi","photonPt","photonEta","photonPhi","Mz","Ptz","Drz","Dphiz","proton1Pt","proton1Rapidity","proton1Phi","proton2Pt","proton2Rapidity","proton2Phi", "zgammaM", "zgammaPt", "zgammaRapidity", "YXyzgamma", "YX", "xi_Resolution(p1&p2)", "diff1", "diff2", "time_Resolution", "ZVertex_resolution",  "isPU1", "isPU2", "ZYdPhi"}; //, "Tmp1", "Tmp2", "sin2", "total"};
    std::vector<TString> HTitles{
        "p_{T}(\\text{leading lepton}) \\ [\\text{GeV}]",
        "\\eta(\\text{leading lepton})",
        "\\Phi(\\text{leading lepton}) \\ [\\text{Rad}]",
        "p_{T}(\\text{sub-leading lepton}) \\ [\\text{GeV}]",
        "\\eta(\\text{sub-leading lepton})",
        "\\Phi(\\text{sub-leading lepton}) \\ [\\text{Rad}]",
        "p_{T}(\\gamma) \\ [\\text{GeV}]",
        "\\eta(\\gamma)",
        "\\Phi(\\gamma) \\ [\\text{Rad}]",
        "M_{\mu\mu} \ [\text{GeV}]",
        "p_{T}(\mu\mu) \\ [\\text{GeV}]",
        "\\Delta R(\mu\mu) \\ [\\text{Rad}]",
        "\\Delta \\Phi(\mu\mu) \\ [\\text{Rad}]",
        "p_{T}(\\text{proton1}) \\ [\\text{GeV}]",
        "\\text{rapidity}(\\text{proton1})",
        "\\Phi(\\text{proton1}) \\ [\\text{Rad}]",
        "p_{T}(\\text{proton2}) \\ [\\text{GeV}]",
        "\\text{rapidity}(\\text{proton2})",
        "\\Phi(\\text{proton2}) \\ [\\text{Rad}]",
        "M_{\\gamma z} \\ [\\text{GeV}]",
        "p_{T}(\\gamma z) \\ [\\text{GeV}]",
        "\\text{rapidity}(\\gamma z)",
        "Y_{X} - Y_{\\gamma z}",
        "Y_{X}",
        "\\text{Resolution}_{\\xi}",
        "|\\xi_{\\text{pps1}} - \\xi_{\\text{cms1}}|",
        "|\\xi_{\\text{pps2}} - \\xi_{\\text{cms2}}|",
        "|Vertex.T_{\\text{pps}} - Vertex.T_{\\text{cms}}| \\ [\\text{ns}]",
        "|Vertex.Z_{\\text{pps}} - Vertex.Z_{\\text{cms}}| \\ [\\text{cm}]",

        "\\text{isPU1}",
        "\\text{isPU2}",
        "\\Delta \\Phi(Z\\gamma)",
        "\\text{Tmp1}",
        "\\text{Tmp2}",
        "\\sin(2)",
        "\\text{total}"
    };

    std::vector<int>    nbins   {30      ,20       ,25       ,20      ,20       ,25       ,30      ,20       ,25       ,30   ,35    ,25    ,15    ,30      ,20       ,25       ,30      ,20       ,25       ,30    ,30   ,20,  40,  40 ,   80,   50,   50,   100,  200, 10, 10, 20, 10, 10, 10, 100};
    std::vector<float> lowEdge  {0       ,-3       ,-4       ,0       ,-3       ,-4       ,0       ,-3       ,-4       ,0    ,0     ,0     ,0       ,0       ,-15       ,-3.2       ,0       ,-15       ,-3.2       ,0       ,0     ,-4      ,  -2e-7,  -2,   -2.,   -0.015*3,   -0.015*3,  -1/*-1.3*/,   -7/*-5e-3*/, -5, -5, -6 , 0, 0, 0, 0};
    std::vector<float> highEdge {2000     ,3        ,4        ,500     ,3        ,4        ,2000     ,3        ,4        ,150   ,3500   ,7     ,4     ,6.5     ,15        ,3.2        ,6.5     ,15        ,3.2       ,8000    ,100    ,4     ,2e-7,   2,  2. ,   0.015*6,    0.015*6,   1/*1.3*/,   7/*5e-3*/, 5, 5,  6, 360000000, 360000000, 1, 1};

    TH1F *Hists[2][8][32+2+2] ;
    std::stringstream name;
    TH1F *h_test;
    for (int i=0;i<channels.size();++i){
      for (int k=0;k<regions.size();++k){
        for (int l=0;l<vars.size();++l){
          name<<channels[i]<<"_"<<regions[k]<<"_"<<vars[l];
          h_test = new TH1F((name.str()).c_str(),(name.str()).c_str(),nbins[l],lowEdge[l],highEdge[l]);
          h_test->StatOverflows(kTRUE);
          h_test->Sumw2(kTRUE);
          Hists[i][k][l] = h_test;
          name.str("");
        }
      }
    }

    TH2F *h_test2;
    std::vector<int>    nbins2X {100, 20, 50, 20, 50, 50, 20};
    std::vector<int>    nbins2Y {100, 20, 50, 20, 50, 50, 20};
    std::vector<float> lowEdge2X {0, -3, 0, -1.3, 0, 0, -20};
    std::vector<float> highEdge2X {2000, 3, 0.23, 1.3, 0.23, 0.32, 20};
    std::vector<float> lowEdge2Y {0, -3, 0, -1.3, -3, 0, -20};
    std::vector<float> highEdge2Y {2000, 3, 0.23, 1.3, 3, 0.32, 20};
    std::vector<TString> vars2 {"MgammazMX", "rapidtygammazYX", "xi cms", "timing", "rapidtygammazXi", "xi_cmsVsxi_pps", "zvertex"};
    std::vector<TString> xax{
    "M(\\gamma z)",
    "\\text{rapidity}(\\gamma z)",
    "\\xi_{\\text{cms1}}",
    "Vertex.T_{\\text{cms}} \\ [\\text{ns}]",
    "\\xi_{\\text{cms1}}",
    "\\xi(\\text{cms})",
    "Vertex.Z_{\\text{cms}} \\ [\\text{cm}]"
    };
    std::vector<TString> yax{
        "M_X",
        "Y_X",
        "\\xi_{\\text{cms2}}",
        "\\frac{(t_{p1} + t_{p2})}{2} - \\frac{Z_{\\text{pps}}}{C} \\ [\\text{ns}]",
        "\\text{rapidity}(\\gamma z)",
        "\\xi_{\\text{pps}}",
        "\\frac{(t_{p2} - t_{p1}) C}{2} \\ [\\text{cm}]"
    };


    TH2F *Hists2[2][8][7] ;

    for (int i=0;i<channels.size();++i){
        for (int k=0;k<regions.size();++k){
            for (int l=0;l<vars2.size();++l){
                name<<channels[i]<<"_"<<regions[k]<<"_"<<vars2[l];
                h_test2 = new TH2F((name.str()).c_str(),(name.str()).c_str(),nbins2X[l],lowEdge2X[l],highEdge2X[l],nbins2Y[l],lowEdge2Y[l],highEdge2Y[l]);
                h_test2->StatOverflows(kTRUE);
                h_test2->Sumw2(kTRUE);
                Hists2[i][k][l] = h_test2;
                name.str("");
            }
        }
    }

    ofstream fSampleCount(path +  "Events_" + ssig + "_" + sxi + "_" + spu + "_" + timepc + "_" + "mumu" + ".txt", std::ofstream::out);  
    
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntries();
    // cout << "nentries: " << nentries << endl;
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      displayProgress(jentry, nentries) ;
      // if (Cut(ientry) < 0) continue;
      // cout << "entry: "<< jentry <<  endl;

    ch = -1;
    region = 0;
    selp = -1;
    selm = -1;
    indP_ = -1;
    indM_ = -1;
    xiPPS = 10000000;
    z_vertex_result = 100000;
    minZVertex = 1000000;
    weight = lumi * cross_section / nentries;

      selectedLeptons = new std::vector<lepton_candidate*>();
      selectedPhotons = new std::vector<photon_candidate*>();
      selectedProtons = new std::vector<proton_candidate*>();

      selectedLeptons->clear();
      selectedPhotons->clear();
      selectedProtons->clear();

      for(int l=0; l<ElectronLoose_size; l++){
  //            if(ElectronLoose_PT[l] < 20)   continue;
  //            if(abs((ElectronLoose_Eta)[l]) > 2.4)    continue;
  //            cout << "Ele Pt: " << ElectronLoose_PT[l] << endl;
          selectedLeptons->push_back(new lepton_candidate(ElectronLoose_PT[l],ElectronLoose_Eta[l],ElectronLoose_Phi[l],ElectronLoose_Charge[l],l,1));
      }

      for(int l=0; l<MuonLoose_size; l++){
//           if(MuonLoose_PT[l] <20 || abs(MuonLoose_Eta[l]) > 2.4)    continue;
//           if(abs(MuonLoose_Eta[l]) > 2.4)    continue;
//           cout << "MuonLoose Pt: " << MuonLoose_PT[l] << endl;
          selectedLeptons->push_back(new lepton_candidate(MuonLoose_PT[l],MuonLoose_Eta[l],MuonLoose_Phi[l],MuonLoose_Charge[l],l,10));
      }

      for(int l=0; l<PhotonTight_size; l++){
          if(PhotonTight_SumPtCharged[l] > 5 || PhotonTight_SumPtCharged[l] < 0 )    continue;
        //   if(PhotonTight_PT[l] < 100 )    continue;
        //   if(abs(PhotonTight_Eta[l]) > 2.4 )    continue;
          //selectedPhotons->push_back(new photon_candidate(PhotonTight_PT[l],PhotonTight_Eta[l],PhotonTight_Phi[l], PhotonTight_E[l] ,0,l,22));
          selectedPhotons->push_back(new photon_candidate(PhotonTight_PT[l],PhotonTight_Eta[l],PhotonTight_Phi[l], PhotonTight_E[l] ,0,l,22));
      }

      sort(selectedLeptons->begin(), selectedLeptons->end(), ComparePtLep);
      sort(selectedPhotons->begin(), selectedPhotons->end(), ComparePtPhoton);

      for(int l=0; l<GenProton_size; l++){
          if(puflag == "nopu")
              if(GenProton_IsPU[l]==1) continue;
          smear = r1->Gaus(0,0.02);
         if(selectedLeptons->size() > 1 && selectedPhotons->size() > 0){
              xi_cms1 = ( ((*selectedLeptons)[0]->p4_).E() + ((*selectedLeptons)[1]->p4_).E() + ((*selectedPhotons)[0]->p4_).E() + ((*selectedLeptons)[0]->p4_).Pz() + ((*selectedLeptons)[1]->p4_).Pz() + ((*selectedPhotons)[0]->p4_).Pz() )/14000;
              xi_cms2 = ( ((*selectedLeptons)[0]->p4_).E() + ((*selectedLeptons)[1]->p4_).E() + ((*selectedPhotons)[0]->p4_).E() - (((*selectedLeptons)[0]->p4_).Pz() + ((*selectedLeptons)[1]->p4_).Pz() + ((*selectedPhotons)[0]->p4_).Pz()) )/14000;
          }
            smear_p = r1->Gaus(0, tr);
            smear = r1->Gaus(0,0.02);
          xiPPS = (1-abs(GenProton_Pz[l])/7000) * (1+smear);
          selectedProtons->push_back(new proton_candidate(GenProton_PT[l],GenProton_Eta[l],GenProton_Phi[l], GenProton_E[l] ,GenProton_Charge[l],l, xiPPS, smear_p));
      }

        // Sort protons based on Pz in ascending order (smallest Pz first)
        sort(selectedProtons->begin(), selectedProtons->end(), ComparePzProton);
        //select first and second protons
        selp = -1;
        selm = -1;
        for (auto it = selectedProtons->begin(); it != selectedProtons->end(); ++it) {
            if (selp == -1 && (*it)->p4_.Pz() > 0) {
                selp = std::distance(selectedProtons->begin(), it);
            }
            if (selm == -1 && (*it)->p4_.Pz() < 0) {
                selm = std::distance(selectedProtons->begin(), it);
            }
            if (selp != -1 && selm != -1) {
                break;
            }
        }

      // lepton size cut
      if(selectedLeptons->size()<2 ||  ((*selectedLeptons)[0]->charge_ * (*selectedLeptons)[1]->charge_ == 1)   ) {
        cleanupVectors(selectedLeptons, selectedPhotons, selectedProtons);
        continue;
      }
      leptonSizeCut++;
      
      // choose channel
      if( (*selectedLeptons)[0]->lep_ == 1 && (*selectedLeptons)[1]->lep_ == 1)  ch = 0;
      if( (*selectedLeptons)[0]->lep_ == 10 && (*selectedLeptons)[1]->lep_ == 10)  ch = 1;
      if(!(ch ==0 || ch ==1) )   {
        cleanupVectors(selectedLeptons, selectedPhotons, selectedProtons);
        continue;
      }

     // Z-boson Pt cut
      if( ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt() < 100 )  {
        cleanupVectors(selectedLeptons, selectedPhotons, selectedProtons);
        continue;
      }
      ZPt_cut++;
      histogram(Hists, Hists2, ch, region, selectedLeptons, selectedPhotons, selectedProtons, weight, GenProton_Rapidity, GenProton_T, GenProton_Z, Vertex_T, Vertex_Z, GenProton_IsPU, timepc, selp, selm);
      region++;

      // photon Pt,eta cut
      sort(selectedPhotons->begin(), selectedPhotons->end(), ComparePtPhoton);
      if(selectedPhotons->size() < 1 || (*selectedPhotons)[0]->pt_ < 200 || abs((*selectedPhotons)[0]->eta_ ) > 2.4 ) {
        cleanupVectors(selectedLeptons, selectedPhotons, selectedProtons);
        continue;
      }
      photonSizeCut++;
      histogram(Hists, Hists2, ch, region, selectedLeptons, selectedPhotons, selectedProtons, weight, GenProton_Rapidity, GenProton_T, GenProton_Z, Vertex_T, Vertex_Z, GenProton_IsPU, timepc, selp, selm);
      region++;
      if(ch==0)    ele_photonSizeCut++;
      if(ch==1)    mu_photonSizeCut++;

      //Mz cut
      float Mz = ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M();
      if (abs(Mz - 90) > 15) {
        cleanupVectors(selectedLeptons, selectedPhotons, selectedProtons);
        continue;
      }
      Mzwindow++;
      histogram(Hists, Hists2, ch, region, selectedLeptons, selectedPhotons, selectedProtons, weight, GenProton_Rapidity, GenProton_T, GenProton_Z, Vertex_T, Vertex_Z, GenProton_IsPU, timepc, selp, selm);
      region++;
      if(ch==0)    ele_Mzwindow++;
      if(ch==1)    mu_Mzwindow++;

      //DeltaPhi(z,y) Cut
      float ZY_dPhi = deltaPhi(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Phi(), (*selectedPhotons)[0]->phi_);
      if (abs(abs(ZY_dPhi) - TMath::Pi()) > 0.1){
        cleanupVectors(selectedLeptons, selectedPhotons, selectedProtons);
        continue;
      }
      zyDeltaPhiCut++;
      histogram(Hists, Hists2, ch, region, selectedLeptons, selectedPhotons, selectedProtons, weight, GenProton_Rapidity, GenProton_T, GenProton_Z, Vertex_T, Vertex_Z, GenProton_IsPU, timepc, selp, selm);
      region++;
      if(ch==0)    ele_zyDeltaPhi++;
      if(ch==1)    mu_zyDeltaPhi++;

      /// PSS Xi cut and delete those of protons outside of range
      for (auto it = selectedProtons->begin(); it != selectedProtons->end(); ) {
          proton_candidate* proton = *it;
          if (proton->xi_ > xi_max || proton->xi_ < xi_min) {
              delete proton;  
              it = selectedProtons->erase(it);  
          } else {
              ++it;
          }
      }


      // Sort protons based on Pz in ascending order (smallest Pz first)
      selp = -1;
      selm = -1;
      sort(selectedProtons->begin(), selectedProtons->end(), ComparePzProton);
      //select first and second protons
      for (auto it = selectedProtons->begin(); it != selectedProtons->end(); ++it) {
          if (selp == -1 && (*it)->p4_.Pz() > 0) {
              selp = std::distance(selectedProtons->begin(), it); 
          }
          if (selm == -1 && (*it)->p4_.Pz() < 0) {
              selm = std::distance(selectedProtons->begin(), it); 
          }
          if (selp != -1 && selm != -1) {
              break;  
          }
      }

      // find two protons with Minimum Z_Vertex
      minZVertex = 100000;

      selp = -1;
      selm = -1;
      indP_ = -1;
      indM_ = -1;
      for (size_t i = 0; i < selectedProtons->size(); ++i) {
          for (size_t j = i + 1; j < selectedProtons->size(); ++j) {
              if (TMath::Sign(1, (*selectedProtons)[i]->p4_.Pz()) == TMath::Sign(1, (*selectedProtons)[j]->p4_.Pz())) {
                  continue;
              }

              if (TMath::Sign(1, (*selectedProtons)[i]->p4_.Pz()) > 0) {
                  selp = i;
                  selm = j;
              } else {
                  selp = j;
                  selm = i;
              }
              z_vertex_result = abs(zt_Vertex_(selectedProtons, GenProton_T, GenProton_Z, selp, selm).first - (Vertex_Z[0] * 0.1));
              if (selp != -1 && selm != -1 && z_vertex_result < minZVertex) {
                  indP_ = selp;
                  indM_ = selm;
                  minZVertex = z_vertex_result;
              }
          }
      }


      /// keep only two protons with Minimum Z_Vertex
      float xip,xim;
      if(indP_ > -1)    xip = (*selectedProtons)[indP_]->xi_;
      if(indM_ > -1)    xim = (*selectedProtons)[indM_]->xi_;

      if(selectedProtons->size()<2) {
        cleanupVectors(selectedLeptons, selectedPhotons, selectedProtons);
        continue;
      }

      for (int index = selectedProtons->size() - 1; index >= 0; --index) {
          if (index != indP_ && index != indM_) {
              proton_candidate* proton = selectedProtons->at(index);
              delete proton;
              selectedProtons->erase(selectedProtons->begin() + index);
          }
      }

      sort(selectedProtons->begin(), selectedProtons->end(), ComparePzProton);
      selp = -1;   selm = -1;
      for (auto it = selectedProtons->begin(); it != selectedProtons->end(); ++it) {
          if (selp == -1 && (*it)->p4_.Pz() > 0) {
              selp = std::distance(selectedProtons->begin(), it);
          }
          if (selm == -1 && (*it)->p4_.Pz() < 0) {
              selm = std::distance(selectedProtons->begin(), it);
          }
          if (selp != -1 && selm != -1) {
              break;
          }
      }

       //proton xi(size) cut
      if(selectedProtons->size()<2  || ((*selectedProtons)[0]->p4_).Pz()/abs(((*selectedProtons)[0]->p4_).Pz()) * ((*selectedProtons)[1]->p4_).Pz()/abs(((*selectedProtons)[1]->p4_).Pz()) ==1 ) {
        cleanupVectors(selectedLeptons, selectedPhotons, selectedProtons);
        continue;
      }
      protonSizeCut++;
      histogram(Hists, Hists2, ch, region, selectedLeptons, selectedPhotons, selectedProtons, weight, GenProton_Rapidity, GenProton_T, GenProton_Z, Vertex_T, Vertex_Z, GenProton_IsPU, timepc, selp, selm);
      region++;
      if(ch==0)    ele_protonSizeCut++;
      if(ch==1)    mu_protonSizeCut++;

      //if ((TMath::Power((xi_smear[(*selectedProtons)[0]->indice_]-((sum_p + sum_pz)/14000)), 2) + TMath::Power((xi_smear[(*selectedProtons)[1]->indice_]-((sum_p - sum_pz)/14000)), 2)) > (TMath::Power((1-((sum_p + sum_pz)/14000)), 2) + TMath::Power((1-((sum_p - sum_pz)/14000)),2)))  continue;
      // xicondition++;
      // if(ch==0)   ele_xicondition++;
      // if(ch==1)   mu_xicondition++;

      xi_cms1 = ( ((*selectedLeptons)[0]->p4_).E() + ((*selectedLeptons)[1]->p4_).E() + ((*selectedPhotons)[0]->p4_).E() + ((*selectedLeptons)[0]->p4_).Pz() + ((*selectedLeptons)[1]->p4_).Pz() + ((*selectedPhotons)[0]->p4_).Pz() )/14000;
      xi_cms2 = ( ((*selectedLeptons)[0]->p4_).E() + ((*selectedLeptons)[1]->p4_).E() + ((*selectedPhotons)[0]->p4_).E() - (((*selectedLeptons)[0]->p4_).Pz() + ((*selectedLeptons)[1]->p4_).Pz() + ((*selectedPhotons)[0]->p4_).Pz()) )/14000;
      xi_cms1 = xi_cms1*(1+smear);
      xi_cms2 = xi_cms2*(1+smear);

      // |xi_pps - xi_cms| cut
      if(abs((*selectedProtons)[0]->xi_ - xi_cms1) > 0.2 && abs((*selectedProtons)[0]->xi_ - xi_cms2) > 0.2 )   {
        cleanupVectors(selectedLeptons, selectedPhotons, selectedProtons);
        continue;
      }
      if(ch==0)   ele_r1XiResulutionCut++;
      if(ch==1)   mu_r1XiResulutionCut++;

      // |xi_pps - xi_cms| cut
      if(abs((*selectedProtons)[1]->xi_ - xi_cms2) > 0.2 && abs((*selectedProtons)[1]->xi_ - xi_cms1) > 0.2 )   {
        cleanupVectors(selectedLeptons, selectedPhotons, selectedProtons);
        continue;
      }
      XiResulutionCut++;
      if(ch==0)   ele_r2XiResulutionCut++;
      if(ch==1)   mu_r2XiResulutionCut++;
      histogram(Hists, Hists2, ch, region, selectedLeptons, selectedPhotons, selectedProtons, weight, GenProton_Rapidity, GenProton_T, GenProton_Z, Vertex_T, Vertex_Z, GenProton_IsPU, timepc, selp, selm);
      region++;
      
      if(ch==1){
        if(SampleCount < 10) {
          //open file for writing
          if (fSampleCount.is_open()){
          // cout << ((*selectedProtons)[0]->p4_).Pz() << endl;
              fSampleCount << "pz1: " << ((*selectedProtons)[0]->p4_).Pz() << "   \t" << "pz2: " <<  ((*selectedProtons)[1]->p4_).Pz() << "   \t" << "r1CMSXi: " << xi_cms1 << "   \t" << "r2CMSXi: " << xi_cms2 << "   \t" << "r1XiResolution: " << std::min(abs(1-(*selectedProtons)[0]->xi_/xi_cms1), abs(1-(*selectedProtons)[0]->xi_/xi_cms2)) << "   \t" << "r2XiResolution: " << min(abs(1-(*selectedProtons)[1]->xi_/xi_cms1), abs(1-(*selectedProtons)[1]->xi_/xi_cms2)) << "\t" << "\n";
          }
          }
          else  fSampleCount.close();
          SampleCount++;
      }


      // Vertex.Z cut
      z_Vertex =  zt_Vertex_(selectedProtons, GenProton_T, GenProton_Z, selp, selm).first;
      float zv_cut;
      if(xiflag == "H")
        zv_cut = 0.65;
      if(xiflag == "V")
        zv_cut = 1.;      
      if( abs((Vertex_Z[0]*0.1) - z_Vertex) > zv_cut )  {
        cleanupVectors(selectedLeptons, selectedPhotons, selectedProtons);
        continue;
      }
      ZVertexCut++;
      if(ch==0)   ele_ZVertexCut++;
      if(ch==1)   mu_ZVertexCut++;
      histogram(Hists, Hists2, ch, region, selectedLeptons, selectedPhotons, selectedProtons, weight, GenProton_Rapidity, GenProton_T, GenProton_Z, Vertex_T, Vertex_Z, GenProton_IsPU, timepc, selp, selm);
      region++;


      // Vertex.T cut
      t_Vertex = zt_Vertex_(selectedProtons, GenProton_T, GenProton_Z, selp, selm).second;
      if(abs(Vertex_T[0]*1e9 - t_Vertex) > 0.2 ){
        cleanupVectors(selectedLeptons, selectedPhotons, selectedProtons);
        continue;
      }
      timingCut++;
      if(ch==0)   ele_timingCut++;
      if(ch==1)   mu_timingCut++;
      histogram(Hists, Hists2, ch, region, selectedLeptons, selectedPhotons, selectedProtons, weight, GenProton_Rapidity, GenProton_T, GenProton_Z, Vertex_T, Vertex_Z, GenProton_IsPU, timepc, selp, selm);
      region++;
    //   if( abs(YX - YZ_Rapidity) > 0.05e-6 )    continue;
    //   RapidityCut++;
    //   histogram(Hists, Hists2, ch, region, selectedLeptons, selectedPhotons, selectedProtons, weight, GenProton_Rapidity, GenProton_T, GenProton_Z, Vertex_T, Vertex_Z, GenProton_IsPU, timepc);
    //   region++;
    //   if(ch==0)   ele_RapidityCut++;
    //   if(ch==1)   mu_RapidityCut++;

      // cout << "xi_cms: " <<  xi_cms1 << "  " << xi_cms2 << "  xi_smear: " << xi_smear[(*selectedProtons)[0]->indice_] << "  " << xi_smear[(*selectedProtons)[1]->indice_] << endl;

      // for (int l=0;l<selectedProtonsplus->size();l++){
      //   delete (*selectedProtonsplus)[l];
      // }

      cleanupVectors(selectedLeptons, selectedPhotons, selectedProtons);
      continue;

      nAccept++;
    }
    cout<<"from "<<nentries<<" evnets, "<<nAccept*weight<<" events are accepted"<<endl;
    TFile file_out (ssig + "_" + sxi + "_" + spu + "_" + timepc + ".root", "RECREATE");
    gROOT->SetBatch(kTRUE);
    int counter = 0;
    for (int i=0;i<channels.size();++i){
        if(ssig != "zy" and i<1)	continue;
        for (int k=0;k<regions.size();++k){
          //if(k<4)	continue;
            for (int l=0;l<vars.size();++l){
              counter++;
              gSystem->MakeDirectory(path + "/" + ssig);
              gSystem->MakeDirectory(path + "/" + ssig + "/" + sxi);
              gSystem->MakeDirectory(path + "/" + ssig + "/" + sxi + "/" + regions[k]);
              Hists[i][k][l]->Write("",TObject::kOverwrite);//
              TCanvas *c = new TCanvas("c2", "", 800, 600);
              TString hname = Hists[i][k][l]->GetName();
              Hists[i][k][l]->SetFillColor(2);
              TString s = ssig + ", " + sxi + ", " + spu;   //std::to_string(counter);
              Hists[i][k][l]->SetTitle(s);
              Hists[i][k][l]->GetXaxis()->SetTitle(HTitles[l]);
              Hists[i][k][l]->GetXaxis()->SetTitleSize(0.05);
              Hists[i][k][l]->GetYaxis()->SetTitleSize(0.05);
              Hists[i][k][l]->GetXaxis()->SetTitleOffset(0.92);
              Hists[i][k][l]->GetYaxis()->SetTitleOffset(0.92);
              Hists[i][k][l]->GetYaxis()->SetTitle("NEvents");
              Hists[i][k][l]->Draw("hist");
              float y_max=1.6*Hists[i][k][l]->GetMaximum();
              Hists[i][k][l]->GetYaxis()->SetRangeUser(Hists[i][k][l]->GetMinimum(), y_max);
              TLatex   *Label_channel = new TLatex(0.2,0.8,channels[i]);
              Label_channel->SetNDC();
              Label_channel->SetTextFont(42);
              Label_channel->Draw("same");
              cout << path + ssig + "/" + sxi + "/" + hname + "_" + spu + "_" + timepc << endl;
              c->Print(path + ssig + "/" + sxi + "/" + regions[k] + "/" + hname + "_" + spu + "_" + timepc + ".png");
              c->Close();
              c->Delete();
              delete c;
            }
        }
    }
//   counter = 0;
    for (int i=0;i<channels.size();++i){
        if(ssig != "zy" and i<1)	continue;
        for (int k=0;k<regions.size();++k){
          //if(k<4)	continue;
            for (int l=0;l<vars2.size();++l){
                counter++;
                gSystem->MakeDirectory(path + "/" + ssig);
                gSystem->MakeDirectory(path + "/" + ssig + "/" + sxi);
                gSystem->MakeDirectory(path + "/" + ssig + "/" + sxi + "/" + regions[k]);
                Hists2[i][k][l]  ->Write("",TObject::kOverwrite);//
                TCanvas *c = new TCanvas("c2", "", 800, 600);
                TString hname = Hists2[i][k][l]->GetName();
                Hists2[i][k][l]->GetXaxis()->SetTitle(xax[l]);
                Hists2[i][k][l]->GetYaxis()->SetTitle(yax[l]);
                Hists2[i][k][l]->GetXaxis()->SetTitleSize(0.05);
                Hists2[i][k][l]->GetYaxis()->SetTitleSize(0.05);
                Hists2[i][k][l]->GetXaxis()->SetTitleOffset(0.92);
                Hists2[i][k][l]->GetYaxis()->SetTitleOffset(0.92);
                TString s =  ssig + ", " + sxi + ", " + spu;    //std::to_string(counter);
                Hists2[i][k][l]->SetTitle(s);
                float y_max=1.6*Hists2[i][k][l]->GetMaximum();
                //Hists2[i][k][l]->GetYaxis()->SetRangeUser(Hists2[i][k][l]->GetMinimum() ,y_max);
                TLatex   *Label_channel = new TLatex(0.2,0.8,channels[i]);
                Label_channel->SetNDC();
                Label_channel->SetTextFont(42);
                Label_channel->Draw("same");
                Hists2[i][k][l]->Draw("col z");
                c->Print(path + ssig + "/" + sxi + "/" + regions[k] + "/" + hname  + "_" + spu + "_" + timepc + ".png");
                delete c;
            }
        }
    }



    file_out.Close();

    //open file for writing
    ofstream fw(path +  ssig + "_" + sxi + "_" + spu + "_" + timepc + "_" + "ee" + ".txt", std::ofstream::out);
    if (fw.is_open()){
        fw << "AllEvents" << "\t" << nentries*weight << "\n";
        fw << "leptonSizeCut" << "\t" << leptonSizeCut*weight << "\n";
        fw << "ZPtCut " << "\t" << ZPt_cut*weight << "\n";
        fw << "photonSizeCut" << "\t" << ele_photonSizeCut*weight << "\n";
        fw << "MzCut" << "\t" << ele_Mzwindow*weight << "\n";
        fw << "zyDeltaPhi" << "\t" << ele_zyDeltaPhi*weight << "\n";
        fw << "protonSizeCut" << "\t" << ele_protonSizeCut*weight << "\n";
        fw << "XiResolutionCut1" << "\t" << ele_r1XiResulutionCut*weight << "\n";
	      fw << "XiResolutionCut2" << "\t" << ele_r2XiResulutionCut*weight << "\n";
        fw << "ZVertexCut" << "\t" << ele_ZVertexCut*weight << "\n";
        fw << "timingCut" << "\t" << ele_timingCut*weight << "\n";
        fw.close();
    }

    //open file for writing
    ofstream fwm(path +  ssig + "_" + sxi + "_" + spu + "_" + timepc + "_" + "mumu" + ".txt", std::ofstream::out);
    if (fwm.is_open()){
        fwm << "AllEvents" << "\t" << nentries*weight << "\n";
        fwm << "leptonSizeCut" << "\t" << leptonSizeCut*weight << "\n";
        fwm << "ZPtCut " << "\t" << ZPt_cut*weight << "\n";
        fwm << "photonSizeCut" << "\t" << mu_photonSizeCut*weight << "\n";
        fwm << "MzCut" << "\t" << mu_Mzwindow*weight << "\n";
        fwm << "zyDeltaPhi" << "\t" << mu_zyDeltaPhi*weight << "\n";
        fwm << "protonSizeCut" << "\t" << mu_protonSizeCut*weight << "\n";
        fwm << "XiResolutionCut1" << "\t" << mu_r1XiResulutionCut*weight << "\n";
        fwm << "XiResolutionCut2" << "\t" << mu_r2XiResulutionCut*weight << "\n";
        fwm << "ZVertexCut" << "\t" << mu_ZVertexCut*weight << "\n";        
        fwm << "timingCut" << "\t" << mu_timingCut*weight << "\n";

        fwm.close();

    }
    else cout << "Problem with opening file";


}



