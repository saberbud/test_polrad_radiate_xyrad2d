#include "g1g2_nlo_cf.h"

extern "C" int readg1g2dataf_(){
  return readg1g2data();
}

extern "C" int findg1g2f_(double *Q2,double *x){
  double Q2in=*Q2;
  double xin=*x;
  return findg1g2(Q2in,xin);
}

extern "C" double reg1pf_(){
  return reg1p();
}

extern "C" double reg2pf_(){
  return reg2p();
}

extern "C" double reg1nf_(){
  return reg1n();
}

extern "C" double reg2nf_(){
  return reg2n();
}
