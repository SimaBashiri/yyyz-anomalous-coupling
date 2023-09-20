#include "proton_candidate.h"

proton_candidate::proton_candidate(float pt_in, float eta_in, float phi_in, float energy_in, int charge_in, int ind_in, float xi_in ){
  pt_ = pt_in;
  eta_ = eta_in;
  phi_ = phi_in;
  energy_ = energy_in;
  charge_ = charge_in;
  p4_.SetPtEtaPhiE(pt_, eta_, phi_, energy_in) ;
  //p4_.SetPxPyPzE(pt_, eta_, phi_, energy_in) ;
  indice_ = ind_in;
  xi_ = xi_in;
}


proton_candidate::~proton_candidate(){}


