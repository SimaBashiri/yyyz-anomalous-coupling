#include "lepton_candidate.h"

lepton_candidate::lepton_candidate(float pt_in, float eta_in, float phi_in, int charge_in, int ind_in, int lep_in ){
  pt_ = pt_in;
  eta_ = eta_in;
  phi_ = phi_in;
  //energy_ = energy_in;
  charge_ = charge_in;
  lep_ = lep_in;
  if(lep_in == 1)  p4_.SetPtEtaPhiM(pt_, eta_, phi_, 0.000511) ;
  if(lep_in == 10)  p4_.SetPtEtaPhiM(pt_, eta_, phi_, 0.10566) ;
  indice_ = ind_in;
}


lepton_candidate::~lepton_candidate(){}


