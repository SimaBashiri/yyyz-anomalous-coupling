#ifndef MY_proton_candidate
#define MY_proton_candidate

#include<cmath>
#include<string>
#include<iostream>
#include<vector>
#include<complex>
#include <TLorentzVector.h>

using namespace std;
//using namespace math;
class proton_candidate {
  
public:
  proton_candidate(float, float, float, float, int, int, float, float );
  ~proton_candidate();
  float pt_;
  float eta_;
  float phi_;
  float energy_;
  int charge_;
  int indice_;
  TLorentzVector p4_;
  float xi_;
  float smearfac_;


private:
  
};

#endif

