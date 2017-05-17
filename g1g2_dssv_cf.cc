#include "g1g2_dssv_cf.h"

extern "C" int readg1g2dssvf_(){
  return readg1g2dssv();
}

extern "C" int findg1g2dssvf_(double *Q2,double *x){
  double Q2in=*Q2;
  double xin=*x;
  return findg1g2dssv(Q2in,xin);
}

extern "C" double reg1pdssvf_(){
  return reg1pdssv();
}

extern "C" double reg2pdssvf_(){
  return reg2pdssv();
}

extern "C" double reg1ndssvf_(){
  return reg1ndssv();
}

extern "C" double reg2ndssvf_(){
  return reg2ndssv();
}
