#include "photon_candidate.h"

photon_candidate::photon_candidate(float pt_in, float eta_in, float phi_in, float energy_in,  int charge_in, int ind_in, int lep_in ){
  pt_ = pt_in;
  eta_ = eta_in;
  phi_ = phi_in;
  energy_ = energy_in;
  charge_ = charge_in;
  lep_ = lep_in;
  if(lep_in == 22)  p4_.SetPtEtaPhiE(pt_, eta_, phi_, energy_) ;
  indice_ = ind_in;
}


photon_candidate::~photon_candidate(){}


