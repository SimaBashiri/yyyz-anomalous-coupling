#ifndef MY_lepton_candidate
#define MY_lepton_candidate

#include<cmath>
#include<string>
#include<iostream>
#include<vector>
#include<complex>
#include <TLorentzVector.h>

using namespace std;
//using namespace math;
class lepton_candidate {
  
public:
  lepton_candidate(float, float, float, int, int, int );
  ~lepton_candidate();
  float pt_;
  float eta_;
  float phi_;
  //float energy_;
  int charge_;
  int indice_;
  int lep_;
  TLorentzVector p4_;


private:
  
};

#endif

