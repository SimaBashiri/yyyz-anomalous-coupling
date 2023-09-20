#ifndef MY_photon_candidate
#define MY_photon_candidate

#include<cmath>
#include<string>
#include<iostream>
#include<vector>
#include<complex>
#include <TLorentzVector.h>

using namespace std;
//using namespace math;
class photon_candidate {
  
public:
  photon_candidate(float, float, float,float, int, int, int );
  ~photon_candidate();
  float pt_;
  float eta_;
  float phi_;
  float energy_;
  int charge_;
  int indice_;
  int lep_;
  TLorentzVector p4_;


private:
  
};

#endif

