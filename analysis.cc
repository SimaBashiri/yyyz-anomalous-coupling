#define analysis_cxx
#include "analysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "lepton_candidate.h"
#include "photon_candidate.h"
#include "proton_candidate.h"
#include "lepton_candidate.cc"
#include "photon_candidate.cc"
#include "proton_candidate.cc"
#include <TRandom3.h>
#include <TLatex.h>
#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <TSystem.h>



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

bool ComparePtLep(lepton_candidate *a, lepton_candidate *b) { return a->pt_ > b->pt_; }
bool ComparePtPhoton(photon_candidate *a, photon_candidate *b) { return a->pt_ > b->pt_; }
bool ComparePzProton(proton_candidate *a, proton_candidate *b) { return abs(a->p4_.Pz()) < abs(b->p4_.Pz()); }
//bool ComparePzProton(proton_candidate *a, proton_candidate *b) { return a->energy_ < b->energy_; }

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

void histogram(TH1F *Hists[2][10][32], TH2F *Hists2[2][10][7], int ch, int cn, std::vector<lepton_candidate*> *selectedLeptons, std::vector<photon_candidate*> *selectedPhotons, std::vector<proton_candidate*> *selectedProtons, float weight, Float_t *GenProton_Rapidity, float YZ_Rapidity, float YX, float diff1, float diff2, float *Vertex_T, float *Vertex_Z, float xi_cms1, float xi_cms2, float z_V, Int_t *GenProton_IsPU, float tVertex, float MX ){

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
    Hists[ch][cn][9]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M(),weight);
    Hists[ch][cn][10]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt(),weight);
    Hists[ch][cn][11]->Fill(deltaR((*selectedLeptons)[0]->eta_,(*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->eta_,(*selectedLeptons)[1]->phi_),weight);
    Hists[ch][cn][12]->Fill(deltaPhi((*selectedLeptons)[0]->phi_,(*selectedLeptons)[1]->phi_),weight);
    if(selectedProtons->size() > 1){
//         float pp_eta = ((*selectedProtons)[0]->p4_ + (*selectedProtons)[1]->p4_).Eta();
      Hists[ch][cn][13]->Fill((*selectedProtons)[0]->pt_,weight);
      Hists[ch][cn][14]->Fill(GenProton_Rapidity[(*selectedProtons)[0]->indice_],weight);
      Hists[ch][cn][15]->Fill((*selectedProtons)[0]->phi_,weight);
      Hists[ch][cn][16]->Fill((*selectedProtons)[1]->pt_,weight);
      Hists[ch][cn][17]->Fill(GenProton_Rapidity[(*selectedProtons)[1]->indice_],weight);
      Hists[ch][cn][18]->Fill((*selectedProtons)[1]->phi_,weight);
    }
    if(selectedPhotons->size() > 0){
      Hists[ch][cn][19]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_ + (*selectedPhotons)[0]->p4_).M(),weight);
      Hists[ch][cn][20]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_ + (*selectedPhotons)[0]->p4_).Pt(),weight);
    }
    if((*selectedProtons).size() > 1 && (selectedPhotons->size() > 0))  Hists[ch][cn][21]->Fill(YZ_Rapidity,weight);
    if((*selectedProtons).size() > 1 && (selectedPhotons->size() > 0))  Hists[ch][cn][22]->Fill( YX - YZ_Rapidity ,weight);
    if(selectedPhotons->size() > 0)  Hists[ch][cn][23]->Fill(YX ,weight);
    if((*selectedProtons).size() > 1 && (selectedPhotons->size() > 0))  Hists[ch][cn][24]->Fill(1-(*selectedProtons)[0]->xi_/xi_cms1 ,weight);
    if((*selectedProtons).size() > 1 && (selectedPhotons->size() > 0))  Hists[ch][cn][24]->Fill(1-(*selectedProtons)[1]->xi_/xi_cms2 ,weight);
    if((*selectedProtons).size() > 1 && (selectedPhotons->size() > 0))  Hists[ch][cn][25]->Fill(diff1,weight);
    if((*selectedProtons).size() > 1 && (selectedPhotons->size() > 0))  Hists[ch][cn][26]->Fill(diff2,weight);
    if((*selectedProtons).size() > 1)  Hists[ch][cn][27]->Fill( 1 - Vertex_T[0]*1e9/tVertex ,weight);
    if((*selectedProtons).size() > 1)  Hists[ch][cn][28]->Fill( 1 - abs(Vertex_Z[0])*100/abs(z_V) ,weight);
    //
    if((*selectedPhotons).size() > 0)  Hists2[ch][cn][0]->Fill(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_ + (*selectedPhotons)[0]->p4_).M(), MX, weight);
    if((*selectedProtons).size() > 1 && selectedLeptons->size() > 1)  Hists2[ch][cn][1]->Fill(0.3, 0.5); //( YZ_Rapidity , YX, weight);
  //     Hists2[ch][3][0][2]->Fill(YZ_Rapidity , pp_rapidity);
  //     Hists2[ch][3][0][3]->Fill( pp_rapidity , YX);
    if(selectedPhotons->size() > 0)  Hists2[ch][cn][2]->Fill(xi_cms1 , xi_cms2, weight);
    if((*selectedProtons).size() > 1)  Hists2[ch][cn][3]->Fill(Vertex_T[0]*1e9 , tVertex, weight);
    if((*selectedProtons).size() > 1 && selectedPhotons->size() > 0 && selectedPhotons->size() > 0)  Hists2[ch][cn][4]->Fill( xi_cms1, YZ_Rapidity, weight);
    if((*selectedProtons).size() > 1 && (selectedPhotons->size() > 0))  Hists2[ch][cn][5]->Fill( xi_cms1, (*selectedProtons)[0]->xi_, weight);
    if((*selectedProtons).size() > 1 && (selectedPhotons->size() > 0))  Hists2[ch][cn][5]->Fill( xi_cms2, (*selectedProtons)[1]->xi_ , weight);
    if((*selectedProtons).size() > 1)  Hists2[ch][cn][6]->Fill(  Vertex_Z[0]*100,  -z_V, weight);
    if((*selectedProtons).size() > 1)  Hists[ch][cn][29]->Fill(GenProton_IsPU[(*selectedProtons)[0]->indice_] , weight);
    if((*selectedProtons).size() > 1)  Hists[ch][cn][30]->Fill(GenProton_IsPU[(*selectedProtons)[1]->indice_] , weight);
    if((*selectedPhotons).size() > 0 && (selectedLeptons->size() > 0)){
        float ZY_dPhi = deltaPhi(((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Phi(), (*selectedPhotons)[0]->phi_);
        Hists[ch][cn][31]->Fill(ZY_dPhi , weight);
    }
}


void analysis::Loop(double cross_section, TString puflag, TString xiflag, TString ssig, int timepc){

    int ch = -1;
    double weight = 1;
    Long64_t nAccept=0;
    float smear = 0;
    float smear_p1 = 0;
    float smear_p2 = 0;
    int pass_e = 0;
    int pass_mu = 0;
    int pass_photon = 0;
    int pass_proton = 0;
    double lumi = 3000e15;
    //double cross_section =  0.2735e-12;

    int leptonCut_pt20 = 0;
    int leptonCut_Eta = 0;
    int ProtonsCut_Pz = 0;
    int ele_ProtonsCut_Pz = 0;
    int mu_ProtonsCut_Pz = 0;
    int protonSizeCut = 0;
    int ele_protonSizeCut = 0;
    int mu_protonSizeCut = 0;
    int channelCut = 0;
    int ele_channelCut = 0;
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
    int ele_XiResulutionCut = 0;
    int mu_XiResulutionCut = 0;
    int RapidityCut = 0;
    int ele_RapidityCut = 0;
    int mu_RapidityCut = 0;
    int Mzwindow = 0;
    int ele_Mzwindow = 0;
    int mu_Mzwindow = 0;
    int xicondition = 0;
    int ele_xicondition = 0;
    int mu_xicondition = 0;
    TString spu;
    TString sxi;
    float xi_min;
    float xi_max;
    int region = 0;
    int ZPt_cut = 0;


    if(puflag == "nop")  spu = "noPU";
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
        sxi = "Horizental" ;
    }

    TString path = "/home/sima/zz/";

           
    std::vector<TString> regions{"channelCut", "photonSizeCut", "ProtonsCut_Pz", "protonSizeCut", "Mzwindow", "CMSXiCut", "ZVertexCut", "timingCut", "XiResulutionCut", "RapidityCut"};
    std::vector<TString> channels{"ee", "mumu"};
    std::vector<TString> vars   {"lep1Pt","lep1Eta","lep1Phi","lep2Pt","lep2Eta","lep2Phi","photonPt","photonEta","photonPhi","Mz","Ptz","Drz","Dphiz","proton1Pt","proton1Rapidity","proton1Phi","proton2Pt","proton2Rapidity","proton2Phi", "zgammaM", "zgammaPt", "zgammaRapidity", "YXyzgamma", "YX", "xi_Resolution(p1&p2)", "diff1", "diff2", "time_Resolution", "ZVertex_resolution",  "isPU", "isnotPU", "ZYdPhi"};

    std::vector<TString> HTitles{"p_{T}(leading lepton) [GeV]","#eta(leading lepton)","#Phi(leading lepton) [Rad]","p_{T}(sub-leading lepton) [GeV]","#eta(sub-leading lepton)","#Phi(sub-leading lepton) [Rad]","p_{T}(#gamma) [GeV]","#eta(#gamma)","#Phi(#gamma) [Rad]","M_{ll} [GeV]","p_{T}(ll) [GeV]","#Delta R(ll) [Rad]","#Delta #Phi(ll) [Rad]","p_{T}(proton1) [GeV]","rapidity(proton1)","#Phi(proton1) [Rad]","p_{T}(proton2) [GeV]","rapidity(proton2)","#Phi(proton2) [Rad]", "M_{#gammaz} [GeV]", "p_{T}(#gammaz) [GeV]", "rapidity(#gammaz)", "Y_{X}-Y_{#gammaz}", "Y_{X}", "Resolution_{#xi}", "|#xi_pps1-#xi_cms1|", "|#xi_pps2-#xi_cms2|", "Resolution_{protons-timing} [ns]", "Resolution_{ZVertex} [cm]", "isPU", "isnotPU", "Z#gammad#Delta#Phi"};

    std::vector<int>    nbins   {30      ,20       ,25       ,20      ,20       ,25       ,30      ,20       ,25       ,30   ,15    ,25    ,15    ,30      ,20       ,25       ,30      ,20       ,25       ,30    ,30   ,20,  40,  40 ,   80,   50,   50,   30,  50, 10, 10, 20};
    std::vector<float> lowEdge  {0       ,-3       ,-4       ,0       ,-3       ,-4       ,0       ,-3       ,-4       ,0    ,0     ,0     ,0       ,0       ,-15       ,-3.2       ,0       ,-15       ,-3.2       ,0       ,0     ,-4      ,  -2e-7,  -2,   -2.,   -0.015*3,   -0.015*3,  -1.3,   -5e-3, -5, -5, -6 };
    std::vector<float> highEdge {1000     ,3        ,4        ,500     ,3        ,4        ,1500     ,3        ,4        ,150   ,2000   ,7     ,4     ,6.5     ,15        ,3.2        ,6.5     ,15        ,3.2       ,4000    ,50    ,4     ,2e-7,   2,  2. ,   0.015*6,    0.015*6,   1.3,   5e-3, 5, 5,  6};


    TH1F *Hists[2][10][32] ;
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
    std::vector<float> lowEdge2X {0, -3, 0, -1.3, 0, 0, 0.5};
    std::vector<float> highEdge2X {2000, 3, 0.23, 1.3, 0.23, 0.32, 0.5};
    std::vector<float> lowEdge2Y {0, -3, 0, -1.3, -3, 0, -2};
    std::vector<float> highEdge2Y {2000, 3, 0.23, 1.3, 3, 0.32, 2};
    std::vector<TString> vars2 {"MgammazMX", "rapidtygammazYX", "xi cms", "timing", "rapidtygammazXi", "xi_cmsVsxi_pps", "zvertex"};
    std::vector<TString> xax{"M(#gammaz)", "rapidity(#gammaz)", "#xi_cms1", "vertex_T[ns]", "#xi cms1", "#xi(cms)", "Vertex_Z [cm]"};
    std::vector<TString> yax{"MX", "YX", "#xi cms2", "(tp1+tp2)/2-Z_{pps}/C [ns]", "rapidity(#gammaz)", "#xi(pps)", "(t1-t2)*c/2 [cm]"};

  
    TH2F *Hists2[2][10][7] ;

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

  
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntries();
    cout << "nentries: " << nentries << endl;
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        displayProgress(jentry, nentries) ;
        // if (Cut(ientry) < 0) continue;
        //cout << "entry: "<< jentry <<  endl;
        //cout << "Tree Num: " << fChain->GetTreeNumber() << endl;

        ch = -1;
        region = 0;

        weight = lumi * cross_section / nentries;
        if (ssig=="zy")
            //zy SM
            weight = lumi * cross_section * Event_Weight[0] * 0.06729 * 1e-12 / nentries;

      std::vector<lepton_candidate*> *selectedLeptons;
      std::vector<photon_candidate*> *selectedPhotons;
      std::vector<proton_candidate*> *selectedProtons;
      std::vector<proton_candidate*> *selectedProtonsplus;
      std::vector<proton_candidate*> *selectedProtonsplusSmear;
      std::vector<proton_candidate*> *selectedProtonsminus;
      std::vector<proton_candidate*> *selectedProtonsminusSmear;
      selectedLeptons = new std::vector<lepton_candidate*>();
      selectedPhotons = new std::vector<photon_candidate*>();
      selectedProtons = new std::vector<proton_candidate*>();
      selectedProtonsplus = new std::vector<proton_candidate*>();
      selectedProtonsminus = new std::vector<proton_candidate*>();
      selectedProtonsplusSmear = new std::vector<proton_candidate*>();
      selectedProtonsminusSmear = new std::vector<proton_candidate*>();

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
//       leptonCut_pt20++;
//       leptonCut_Eta++;
      for(int l=0; l<PhotonTight_size; l++){
          if(PhotonTight_PT[l] < 100 )    continue;
          if(abs(PhotonTight_Eta[l]) > 2.4 )    continue;
          selectedPhotons->push_back(new photon_candidate(PhotonTight_PT[l],PhotonTight_Eta[l],PhotonTight_Phi[l], PhotonTight_E[l] ,0,l,22));
      }

      float pz_min = (1-xi_max)*7000;
      float pz_max = (1-xi_min)*7000;
      std::vector<float> xi;
      std::vector<float> xi_smear;
//       TRandom3 *r = new TRandom3();
//       r->SetSeed(19680801);
      int n_p = 0;
      int n_n = 0;

      sort(selectedLeptons->begin(), selectedLeptons->end(), ComparePtLep);
      sort(selectedPhotons->begin(), selectedPhotons->end(), ComparePtPhoton);

      if(selectedLeptons->size()>2){
          for(int l=0; l<selectedLeptons->size(); l++){
              if( (*selectedLeptons)[0]->lep_ != (*selectedLeptons)[1]->lep_ ){

                  if( (*selectedLeptons)[0]->lep_ == (*selectedLeptons)[2]->lep_ )  selectedLeptons->erase(selectedLeptons->begin()+1);
                  else if( (*selectedLeptons)[1]->lep_ == (*selectedLeptons)[2]->lep_ )  selectedLeptons->erase(selectedLeptons->begin());
              }
          }
      }

      if(selectedLeptons->size()<2 ||  ((*selectedLeptons)[0]->charge_ * (*selectedLeptons)[1]->charge_ == 1)   ) {
          for (unsigned int l=0;l<selectedLeptons->size();l++){
              delete (*selectedLeptons)[l];
          }
          selectedLeptons->clear();
          selectedLeptons->shrink_to_fit();
          delete selectedLeptons;
          continue;
      }
      leptonSizeCut++;


      if( ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).Pt() < 100 )  continue;
      ZPt_cut++;

      if( (*selectedLeptons)[0]->lep_ == 1 && (*selectedLeptons)[1]->lep_ == 1)  ch = 0;
      if( (*selectedLeptons)[0]->lep_ == 10 && (*selectedLeptons)[1]->lep_ == 10)  ch = 1;
      if(!(ch ==0 || ch ==1) )   continue;
      histogram(Hists, Hists2, ch, region, selectedLeptons, selectedPhotons, selectedProtons, weight, GenProton_Rapidity, 0, 0, 0, 0, Vertex_T, Vertex_Z, 0, 0, 0, GenProton_IsPU, 0, 0);
      region++;
      channelCut++;

      if(ch==0)    ele_channelCut++;
      if(ch==1)    mu_channelCut++;

      if((*selectedLeptons)[0]->lep_ == 1  && (*selectedLeptons)[1]->lep_ == 1)   pass_e++;
      if((*selectedLeptons)[0]->lep_ == 10  && (*selectedLeptons)[1]->lep_ == 10)   pass_mu++;
      if(selectedPhotons->size()<1 ) {
          for (unsigned int l=0;l<selectedPhotons->size();l++){
              delete (*selectedPhotons)[l];
          }
          selectedPhotons->clear();
          selectedPhotons->shrink_to_fit();
          delete selectedPhotons;
          continue;
      }
      photonSizeCut++;

      float xi_cms1 = ( ((*selectedLeptons)[0]->p4_).E() + ((*selectedLeptons)[1]->p4_).E() + ((*selectedPhotons)[0]->p4_).E() + ((*selectedLeptons)[0]->p4_).Pz() + ((*selectedLeptons)[1]->p4_).Pz() + ((*selectedPhotons)[0]->p4_).Pz() )/14000;
      float xi_cms2 = ( ((*selectedLeptons)[0]->p4_).E() + ((*selectedLeptons)[1]->p4_).E() + ((*selectedPhotons)[0]->p4_).E() - (((*selectedLeptons)[0]->p4_).Pz() + ((*selectedLeptons)[1]->p4_).Pz() + ((*selectedPhotons)[0]->p4_).Pz()) )/14000;

      float MX = 14000 * TMath::Sqrt(xi_cms1 * xi_cms2);
      float YX = 0.5 * TMath::Log(xi_cms1 / xi_cms2);
      double YZ_Rapidity = Rapidity( ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_ + (*selectedPhotons)[0]->p4_).E(),  ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_ + (*selectedPhotons)[0]->p4_).Pz() );


      histogram(Hists, Hists2, ch, region, selectedLeptons, selectedPhotons, selectedProtons, weight, GenProton_Rapidity, YZ_Rapidity, YX, 0, 0, Vertex_T, Vertex_Z, xi_cms1, xi_cms2, 0, GenProton_IsPU, 0, MX);
      region++;
      if(ch==0)    ele_photonSizeCut++;
      if(ch==1)    mu_photonSizeCut++;

      for(int l=0; l<GenProton_size; l++){
          //if(GenProton_Mass[l]<0) continue;
          //if(GenProton_IsPU[l]==1) continue;
          if(puflag == "nop")
              if(GenProton_IsPU[l]==1) continue;
        //   if( -9 < GenProton_Rapidity[l]  && GenProton_Rapidity[l] < 9)	continue;

          float PE = TMath::Sqrt(GenProton_Px[l]*GenProton_Px[l] + GenProton_Py[l]*GenProton_Py[l] + GenProton_Pz[l]*GenProton_Pz[l] + GenProton_Mass[l]*GenProton_Mass[l]);

          xi.push_back(1-abs(GenProton_Pz[l])/7000);
          smear = r->Gaus(0,0.02);
          xi_smear.push_back( xi[l]*(1+smear) );

  //           if (!(abs(GenProton_Pz[l]) > pz_min and (abs(GenProton_Pz[l])) < pz_max))   continue;
          if (! (xi[l] > xi_min and xi[l] < xi_max) )   continue;

          if(GenProton_Pz[l] > 0){
              // if( abs(xi_cms1 - (1-abs(GenProton_Pz[l])/7000)*(1+smear)) > 0.02  )  continue;
              selectedProtonsplus->push_back(new proton_candidate(GenProton_PT[l],GenProton_Eta[l],GenProton_Phi[l], GenProton_E[l] ,GenProton_Charge[l],l, xi_smear[l]));
              TLorentzVector pp ;
              pp.SetPxPyPzE(((*selectedProtonsplus)[n_p]->p4_).Px(), ((*selectedProtonsplus)[n_p]->p4_).Py(), 7000*(1-xi_smear[l]), ((*selectedProtonsplus)[n_p]->p4_).E());
              selectedProtonsplusSmear->push_back(new proton_candidate(pp.Pt(),pp.Eta(),pp.Phi(),pp.E() ,GenProton_Charge[l],l, xi_smear[l]));
              n_p++;
          }

          if(GenProton_Pz[l] < 0){
              //   if( abs(xi_cms2 - (1-abs(GenProton_Pz[l])/7000)*(1+smear)) > 0.02 )  continue;
              selectedProtonsminus->push_back(new proton_candidate(GenProton_PT[l],GenProton_Eta[l],GenProton_Phi[l], GenProton_E[l] ,GenProton_Charge[l],l, xi_smear[l]));
              TLorentzVector pn ;
              pn.SetPxPyPzE(((*selectedProtonsminus)[n_n]->p4_).Px(), ((*selectedProtonsminus)[n_n]->p4_).Py(), 7000*(1-xi_smear[l]), ((*selectedProtonsminus)[n_n]->p4_).E());
              selectedProtonsminusSmear->push_back(new proton_candidate(pn.Pt(),pn.Eta(),pn.Phi(),pn.E() ,GenProton_Charge[l],l, xi_smear[l]));
                n_n++;
          }
      }
      ProtonsCut_Pz++;
      histogram(Hists, Hists2, ch, region, selectedLeptons, selectedPhotons, selectedProtons, weight, GenProton_Rapidity, YZ_Rapidity, YX, 0, 0, Vertex_T, Vertex_Z, xi_cms1, xi_cms2, 0, GenProton_IsPU, 0, MX);
      region++;
      if(ch==0)    ele_ProtonsCut_Pz++;
      if(ch==1)    mu_ProtonsCut_Pz++;

      sort(selectedProtonsplus->begin(), selectedProtonsplus->end(), ComparePzProton);
      sort(selectedProtonsminus->begin(), selectedProtonsminus->end(), ComparePzProton);
      sort(selectedProtonsplusSmear->begin(), selectedProtonsplusSmear->end(), ComparePzProton);
      sort(selectedProtonsminusSmear->begin(), selectedProtonsminusSmear->end(), ComparePzProton);

      if(selectedProtonsplus->size() > 0)  selectedProtons->push_back( (*selectedProtonsplus)[0] );
      if(selectedProtonsminus->size() > 0)  selectedProtons->push_back( (*selectedProtonsminus)[0] );

      if(selectedProtons->size()<2  || ((*selectedProtons)[0]->p4_).Pz()/abs(((*selectedProtons)[0]->p4_).Pz()) * ((*selectedProtons)[1]->p4_).Pz()/abs(((*selectedProtons)[1]->p4_).Pz()) ==1 ) {
          for (unsigned int l=0;l<selectedProtons->size();l++){
              delete (*selectedProtons)[l];
          }
          selectedProtons->clear();
          selectedProtons->shrink_to_fit();
          delete selectedProtons;
          continue;
      }
      float pp_rapidity = 0.5*TMath::Log((((*selectedProtons)[0]->p4_ + (*selectedProtons)[1]->p4_).E() + ((*selectedProtons)[0]->p4_ + (*selectedProtons)[1]->p4_).Pz())/(((*selectedProtons)[0]->p4_ + (*selectedProtons)[1]->p4_).E() - ((*selectedProtons)[0]->p4_ + (*selectedProtons)[1]->p4_).Pz()));
      float sum_p = TMath::Sqrt( ((*selectedLeptons)[0]->p4_).Px()*((*selectedLeptons)[0]->p4_).Px() + ((*selectedLeptons)[0]->p4_).Py()*((*selectedLeptons)[0]->p4_).Py() +
      ((*selectedLeptons)[0]->p4_).Pz()*((*selectedLeptons)[0]->p4_).Pz() ) + TMath::Sqrt( ((*selectedLeptons)[1]->p4_).Px()*((*selectedLeptons)[1]->p4_).Px() + ((*selectedLeptons)[1]->p4_).Py()*((*selectedLeptons)[1]->p4_).Py() +
      ((*selectedLeptons)[1]->p4_).Pz()*((*selectedLeptons)[1]->p4_).Pz() ) + TMath::Sqrt( ((*selectedPhotons)[0]->p4_).Px()*((*selectedPhotons)[0]->p4_).Px() + ((*selectedPhotons)[0]->p4_).Py()*((*selectedPhotons)[0]->p4_).Py() +
      ((*selectedPhotons)[0]->p4_).Pz()*((*selectedPhotons)[0]->p4_).Pz() );

      double sum_pz = ((*selectedLeptons)[0]->p4_).Pz() + ((*selectedLeptons)[1]->p4_).Pz() + ((*selectedPhotons)[0]->p4_).Pz();
      float diff1 = xi_smear[(*selectedProtons)[0]->indice_]-((sum_p + sum_pz)/14000);
      float diff2 = xi_smear[(*selectedProtons)[1]->indice_]-((sum_p - sum_pz)/14000);

      float C = 30;
      float tp1, tp2;
      float tr = (pow(10,9))* timepc * 1E-12;
      float tr2 = (pow(10,9))* timepc * 1E-12;

      //TRandom3 *r1 = new TRandom3();
      //r1->SetSeed(19680801);
      smear_p1 = r1->Gaus(0, tr);
      smear_p2 = r1->Gaus(0, tr2);

      if( ((*selectedProtons)[0]->p4_).Pz()>0 )
          tp1 = (1e9*GenProton_T[(*selectedProtons)[0]->indice_]+(23400-100*GenProton_Z[(*selectedProtons)[0]->indice_])/30) + smear_p1 ;
      else
          tp1 = (1e9*GenProton_T[(*selectedProtons)[0]->indice_]+(23400+100*GenProton_Z[(*selectedProtons)[0]->indice_])/30) + smear_p1;

      if( ((*selectedProtons)[1]->p4_).Pz()>0 )
          tp2 = (1e9*GenProton_T[(*selectedProtons)[1]->indice_]+(23400-100*GenProton_Z[(*selectedProtons)[1]->indice_])/30) + smear_p2;
      else
          tp2 = (1e9*GenProton_T[(*selectedProtons)[1]->indice_]+(23400+100*GenProton_Z[(*selectedProtons)[1]->indice_])/30) + smear_p2;

      float tVertex = ((tp1 + tp2)/2 - 23400/C);
      float z_V = (tp1 - tp2)*C/2;
      protonSizeCut++;
      histogram(Hists, Hists2, ch, region, selectedLeptons, selectedPhotons, selectedProtons, weight, GenProton_Rapidity, YZ_Rapidity, YX, diff1, diff2, Vertex_T, Vertex_Z, xi_cms1, xi_cms2, z_V, GenProton_IsPU, tVertex, MX);
      region++;
      if(ch==0)    ele_protonSizeCut++;
      if(ch==1)    mu_protonSizeCut++;

      float Mz = ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M();
      if( 90 - 15 > Mz  or  Mz > 90 + 15 )	continue;
      Mzwindow++;
      histogram(Hists, Hists2, ch, region, selectedLeptons, selectedPhotons, selectedProtons, weight, GenProton_Rapidity, YZ_Rapidity, YX, diff1, diff2, Vertex_T, Vertex_Z, xi_cms1, xi_cms2, z_V, GenProton_IsPU, tVertex, MX);
      region++;
      if(ch==0)    ele_Mzwindow++;
      if(ch==1)    mu_Mzwindow++;

      //if ((TMath::Power((xi_smear[(*selectedProtons)[0]->indice_]-((sum_p + sum_pz)/14000)), 2) + TMath::Power((xi_smear[(*selectedProtons)[1]->indice_]-((sum_p - sum_pz)/14000)), 2)) > (TMath::Power((1-((sum_p + sum_pz)/14000)), 2) + TMath::Power((1-((sum_p - sum_pz)/14000)),2)))  continue;
      xicondition++;
      if(ch==0)   ele_xicondition++;
      if(ch==1)   mu_xicondition++;

      cout << "xi " << xi_cms1 << "  " << xi_cms2  << "  " << xi_smear[(*selectedProtons)[0]->indice_]<< "  " << (*selectedProtons)[1]->xi_ <<   endl;

      xi_cms1 = xi_cms1*(1+smear);
      xi_cms2 = xi_cms2*(1+smear);


      if(xi_cms1 < xi_min || xi_cms1 > xi_max) continue;
      if(ch==0)   ele_r1cms_xiCut++;
      if(ch==1)   mu_r1cms_xiCut++;
      if(xi_cms2 < xi_min || xi_cms2 > xi_max) continue;
//       cms_xiCut++;
      if(ch==0)   ele_r2cms_xiCut++;
      if(ch==1)   mu_r2cms_xiCut++;
      histogram(Hists, Hists2, ch, region, selectedLeptons, selectedPhotons, selectedProtons, weight, GenProton_Rapidity, YZ_Rapidity, YX, diff1, diff2, Vertex_T, Vertex_Z, xi_cms1, xi_cms2, z_V, GenProton_IsPU, tVertex, MX);
      region++;


      if(1 - abs(Vertex_Z[0])*100/abs(z_V) > 0.002 )  continue;
      ZVertexCut++;
      if(ch==0)   ele_ZVertexCut++;
      if(ch==1)   mu_ZVertexCut++;

      histogram(Hists, Hists2, ch, region, selectedLeptons, selectedPhotons, selectedProtons, weight, GenProton_Rapidity, YZ_Rapidity, YX, diff1, diff2, Vertex_T, Vertex_Z, xi_cms1, xi_cms2, z_V, GenProton_IsPU, tVertex, MX);
      region++;

      if(( 1 - Vertex_T[0]*1e9/tVertex) < 0.8 ){  //continue;
          timingCut++;
          if(ch==0)   ele_timingCut++;
          if(ch==1)   mu_timingCut++;
          histogram(Hists, Hists2, ch, region, selectedLeptons, selectedPhotons, selectedProtons, weight, GenProton_Rapidity, YZ_Rapidity, YX, diff1, diff2, Vertex_T, Vertex_Z, xi_cms1, xi_cms2, z_V, GenProton_IsPU, tVertex, MX);
          region++;
      }

      if(abs(1-(*selectedProtons)[0]->xi_/xi_cms1) > 0.15 )   continue;
      if(abs(1-(*selectedProtons)[1]->xi_/xi_cms2) > 0.15 )   continue;
      XiResulutionCut++;
      if(ch==0)   ele_XiResulutionCut++;
      if(ch==1)   mu_XiResulutionCut++;

      histogram(Hists, Hists2, ch, region, selectedLeptons, selectedPhotons, selectedProtons, weight, GenProton_Rapidity, YZ_Rapidity, YX, diff1, diff2, Vertex_T, Vertex_Z, xi_cms1, xi_cms2, z_V, GenProton_IsPU, tVertex, MX);
      region++;


      if( abs(YX - YZ_Rapidity) > 0.05e-6 )    continue;
      RapidityCut++;
      histogram(Hists, Hists2, ch, region, selectedLeptons, selectedPhotons, selectedProtons, weight, GenProton_Rapidity, YZ_Rapidity, YX, diff1, diff2, Vertex_T, Vertex_Z, xi_cms1, xi_cms2, z_V, GenProton_IsPU, tVertex, MX);
      region++;
      if(ch==0)   ele_RapidityCut++;
      if(ch==1)   mu_RapidityCut++;

      cout << "protons mass :" << ((*selectedProtons)[0]->p4_ + (*selectedProtons)[1]->p4_).M() << endl;
      cout << "xi_cms: " <<  xi_cms1 << "  " << xi_cms2 << "  xi_smear: " << xi_smear[(*selectedProtons)[0]->indice_] << "  " << xi_smear[(*selectedProtons)[1]->indice_] << endl;



      for (int l=0;l<selectedLeptons->size();l++){
        delete (*selectedLeptons)[l];
      }
      for (int l=0;l<selectedPhotons->size();l++){
        delete (*selectedPhotons)[l];
      }
      for (int l=0;l<selectedProtons->size();l++){
        delete (*selectedProtons)[l];
      }
      selectedLeptons->clear();
      selectedLeptons->shrink_to_fit();
      delete selectedLeptons;
      selectedPhotons->clear();
      selectedPhotons->shrink_to_fit();
      delete selectedPhotons;
      selectedProtons->clear();
      selectedProtons->shrink_to_fit();
      delete selectedProtons;
      nAccept++;

    }

    cout<<"from "<<nentries<<" evnets, "<<nAccept*weight<<" events are accepted"<<endl;
    TFile file_out (ssig + "_" + sxi + "_" + spu + "_" + timepc + ".root", "RECREATE");
    gROOT->SetBatch(kTRUE);
    int counter = 0;
    for (int i=0;i<channels.size();++i){
        if(i<1)	continue;
        for (int k=0;k<regions.size();++k){
          //if(k<4)	continue;
            for (int l=0;l<vars.size();++l){
              counter++;
              gSystem->MakeDirectory(path + "/" + ssig + "/" + sxi + "/" + regions[k]);
              Hists[i][k][l]->Write("",TObject::kOverwrite);//
              TCanvas *c = new TCanvas("c2", "", 800, 600);
              //TCanvas c = TCanvas("c2", "", 800, 600);

              //TPad *p = new TPad("pad1", "pad1", 0, 0.315, 1, 0.99 , 0);
              //p->cd();
              //p->Draw();
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
        if(i<1)	continue;
        for (int k=0;k<regions.size();++k){
          //if(k<4)	continue;
            for (int l=0;l<vars2.size();++l){
                counter++;
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
  
//   h_llyInvMass->Write("",TObject::kOverwrite);//
//   h_llyrapidty->Write("",TObject::kOverwrite);//
//   h_ppInvMass->Write("",TObject::kOverwrite);//
//   TCanvas *c = new TCanvas();
//   h_llyInvMass->Draw("col z");
//   c
  

    file_out.Close();
    cout << "Finish" << endl;
    cout << "pass_e  " << pass_e*weight << endl;
    cout << "pass_mu  " << pass_mu*weight << endl;
    cout << "pass_photon  " << pass_photon*weight << endl;
    cout << "pass_proton  " << pass_proton*weight << endl;
    cout << "integral: " << Hists[0][0][0]->Integral() << endl;
    cout << "integral: " << Hists[1][0][0]->Integral() << endl;
    cout << "weight: " << weight << endl;
    cout << "Accept: " << nAccept << endl;

    cout << "leptonCut_pt20: " << leptonCut_pt20*weight << "  leptonCut_Eta: " << leptonCut_Eta*weight << "ZPt_cut: " << ZPt_cut*weight << "  ProtonsCut_Pz: " << ProtonsCut_Pz*weight << "  channelCut: " << channelCut*weight << "  phtonSizeCut: " << photonSizeCut*weight << " leptonSizeCut: " << leptonSizeCut*weight  << " protonSizeCut: " << protonSizeCut*weight << " Mzwindow: " << Mzwindow*weight << "  cms_xiCut: " << cms_xiCut*weight << " xiCondition: " << xicondition*weight << "  ZVertexCut: " << ZVertexCut*weight << "  timingCut: " << timingCut*weight << "  XiResulutionCut: " << XiResulutionCut*weight << "  RapidityCut: " << RapidityCut*weight << endl;

    //open file for writing
    ofstream fw(path +  ssig + "_" + sxi + "_" + spu + "_" + timepc + "_" + "ee" + ".txt", std::ofstream::out);
    if (fw.is_open()){
        fw << "AllEvents" << "\t" << nentries*weight << "\n";
    //     fw << "lepCutPtEta30" << "\t" << leptonCut_pt20*weight << "\n";
        fw << "leptonSizeCut" << "\t" << leptonSizeCut*weight << "\n";
        fw << "ZPtCut " << "\t" << ZPt_cut*weight << "\n";
        fw << "channelCut" << "\t" << ele_channelCut*weight << "\n";
        fw << "photonSizeCut" << "\t" << ele_photonSizeCut*weight << "\n";
        fw << "PPSXiCut" << "\t" << ele_ProtonsCut_Pz*weight << "\n";
        fw << "protonSizeCut" << "\t" << ele_protonSizeCut*weight << "\n";
        fw << "MzCut" << "\t" << ele_Mzwindow*weight << "\n";
        fw << "r1CMSXiCut" << "\t" << ele_r1cms_xiCut*weight << "\n";
        fw << "r2CMSXiCut" << "\t" << ele_r2cms_xiCut*weight << "\n";
        fw << "ZVertexCut" << "\t" << ele_ZVertexCut*weight << "\n";
        //fw << "timingCut" << "\t" << timingCut*weight << "\n";
        fw << "XiResolutionCut" << "\t" << ele_XiResulutionCut*weight << "\n";

        fw.close();
    }
    
    //open file for writing
    ofstream fwm(path +  ssig + "_" + sxi + "_" + spu + "_" + timepc + "_" + "mumu" + ".txt", std::ofstream::out);
    if (fwm.is_open()){
        fwm << "AllEvents" << "\t" << nentries*weight << "\n";
    //     fwm << "lepCutPtEta20" << "\t" << leptonCut_pt20*weight << "\n";
        fwm << "leptonSizeCut" << "\t" << leptonSizeCut*weight << "\n";
        fwm << "ZPtCut " << "\t" << ZPt_cut*weight << "\n";
        fwm << "channelCut" << "\t" << mu_channelCut*weight << "\n";
        fwm << "photonSizeCut" << "\t" << mu_photonSizeCut*weight << "\n";
        fwm << "PPSXiCut" << "\t" << mu_ProtonsCut_Pz*weight << "\n";
        fwm << "protonSizeCut" << "\t" << mu_protonSizeCut*weight << "\n";
        fwm << "MzCut" << "\t" << mu_Mzwindow*weight << "\n";
//         fwm << "XiCut" << "\t" << mu_xiCut*weight << "\n";
        fwm << "r1CMSXiCut" << "\t" << mu_r1cms_xiCut*weight << "\n";
        fwm << "r2CMSXiCut" << "\t" << mu_r2cms_xiCut*weight << "\n";
        fwm << "ZVertexCut" << "\t" << mu_ZVertexCut*weight << "\n";
        //fw << "timingCut" << "\t" << timingCut*weight << "\n";
        fwm << "XiResolutionCut" << "\t" << mu_XiResulutionCut*weight << "\n";

        fwm.close();

    }
    else cout << "Problem with opening file";


}


