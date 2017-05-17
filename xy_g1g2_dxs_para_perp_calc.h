#include <iostream>
#include <cstring>
#include <cmath>

int g1g2_dxs_solve(double e0,double ep,double theta,double *xs,double *g12){
  double M=0.93827;
  double ALPH=1./137.035999;
  double Q2=4.*e0*ep*pow(sin(theta/2.),2);
  double nu=e0-ep;
  double xbj=Q2/2.0/M/nu;

  double xsperp=xs[0];double xspara=xs[1];

  double c1=1000.*389.379*4.0*ALPH*ALPH*ep;
  double c2=1./(M*Q2*e0*nu);
  double c3=ep*sin(theta);

  double a11,a12,a21,a22;
  a11=(e0+ep*cos(theta))*c1*c2;a12=-2.0*M*xbj*c1*c2;
  a21=c1*c2*c3;a22=(2.0*e0/nu)*c1*c2*c3;

  double g1=1.;double g2=2.;
  g1=(a12*xsperp-a22*xspara)/(a12*a21-a11*a22);
  g2=(a21*xspara-a11*xsperp)/(a12*a21-a11*a22);

  if(g1!=g1||g2!=g2){
    g12[0]=0.;g12[1]=0.;
    return 0;
  }

  g12[0]=g1;g12[1]=g2;

  return 1;
}

double g2_dxs_solve(double e0,double ep,double theta,double xsperp,double g1){
  double M=0.93827;
  double ALPH=1./137.035999;
  double Q2=4.*e0*ep*pow(sin(theta/2.),2);
  double nu=e0-ep;
  double xbj=Q2/2.0/M/nu;

  double c1=1000.*389.379*4.0*ALPH*ALPH*ep;
  double c2=1./(M*Q2*e0*nu);
  double c3=ep*sin(theta);

  double a11,a12,a21,a22;
  a11=(e0+ep*cos(theta))*c1*c2;a12=-2.0*M*xbj*c1*c2;
  a21=c1*c2*c3;a22=(2.0*e0/nu)*c1*c2*c3;

  double g2=(xsperp-a21*g1)/a22;

  return g2;
}

int g1g2_dxs_calc(double e0,double ep,double theta,double *xs,double *g12){
  double M=0.93827;
  double ALPH=1./137.035999;
  double Q2=4.*e0*ep*pow(sin(theta/2.),2);
  double nu=e0-ep;
  double xbj=Q2/2.0/M/nu;

  double c1=1000.*389.379*4.0*ALPH*ALPH*ep;
  double c2=1./(M*Q2*e0*nu);
  double c3=ep*sin(theta);

  double a11,a12,a21,a22;
  a11=(e0+ep*cos(theta))*c1*c2;a12=-2.0*M*xbj*c1*c2;
  a21=c1*c2*c3;a22=(2.0*e0/nu)*c1*c2*c3;

  xs[1]=a11*g12[0]+a12*g12[1];
  xs[0]=a21*g12[0]+a22*g12[1];

  return 1;
}
