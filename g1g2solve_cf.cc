#include "xy_g1g2_dxs_para_perp_calc.h"

extern "C" double g1_perp_para_sol__(double *e0_in,double *ep_in,double *theta_in,double *xsperp,double *xspara){
  double e0=*e0_in;double ep=*ep_in;double theta=*theta_in;
  double xs[2];
  xs[0]=*xsperp;xs[1]=*xspara;
  xs[0]=xs[0]*1000.;xs[1]=xs[1]*1000.;
  double g12[2];
  int ig12=g1g2_dxs_solve(e0,ep,theta,xs,g12);

  double g1=g12[0];

  return g1;
}

extern "C" double g2_perp_para_sol__(double *e0_in,double *ep_in,double *theta_in,double *xsperp,double *xspara){
  double e0=*e0_in;double ep=*ep_in;double theta=*theta_in;
  double xs[2];
  xs[0]=*xsperp;xs[1]=*xspara;
  xs[0]=xs[0]*1000.;xs[1]=xs[1]*1000.;
  double g12[2];
  int ig12=g1g2_dxs_solve(e0,ep,theta,xs,g12);

  double g2=g12[1];

  return g2;
}

extern "C" double g2_perp_g1_sol__(double *e0_in,double *ep_in,double *theta_in,double *xsperp,double *g1_in){
  double e0=*e0_in;double ep=*ep_in;double theta=*theta_in;
  double xs[2];
  xs[0]=*xsperp;
  xs[0]=xs[0]*1000.;
  double g1=*g1_in;

  double g2=g2_dxs_solve(e0,ep,theta,xs[0],g1);

  return g2;
}

