#include "g1g2_grsv_cf.h"

extern "C" int readg1g2grsvf_(){
  return readg1g2grsv();
}

extern "C" int findg1g2grsvf_(double *Q2,double *x){
  double Q2in=*Q2;
  double xin=*x;
  return findg1g2grsv(Q2in,xin);
}

extern "C" double reg1pgrsvf_(){
  return reg1pgrsv();
}

extern "C" double reg2pgrsvf_(){
  return reg2pgrsv();
}

extern "C" double reg1ngrsvf_(){
  return reg1ngrsv();
}

extern "C" double reg2ngrsvf_(){
  return reg2ngrsv();
}
