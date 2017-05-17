#include "Lstructure.h"

extern "C" double g1cc_(double *x, double *Q2, double *np, double* nn, double* type){
        double var[4];
        var[0]=*x;var[1]=*Q2;var[2]=0.;var[3]=0.;
        double AZc[4]={1.0,1.0,1.0,1.0};
        AZc[0]=*np;AZc[1]=*nn;

        // std::cout << var[0] << " " << var[1] << " " << AZc[0] << " " << AZc[1] << std::endl;

        double g1temp[6],result;

        if(*type<1){
          Lstructure::g1N(AZc,var,g1temp,"dssv");
        }else{
          Lstructure::g1N(AZc,var,g1temp,"grsv");
        }

        // std::cout << " " << g1temp[0] << std::endl;

        // return g1temp[0];

        result=0.5*(4.0/9.0)*g1temp[0]+0.5*(1.0/9.0)*g1temp[1]+0.5*(1.0/9.0)*g1temp[2]+0.5*(4.0/9.0)*g1temp[3]+0.5*(1.0/9.0)*g1temp[4]+0.5*(1.0/9.0)*g1temp[5];

        return result;

}

extern "C" double g2cc_(double *x, double *Q2, double *np, double* nn, double* type, double* nin){
        double var[4];
        var[0]=*x;var[1]=*Q2;var[2]=0.;var[3]=0.;
        double AZc[4]={1.0,1.0,1.0,1.0};
        AZc[0]=*np;AZc[1]=*nn;
        double nint=*nin;
        int n=nint;

        double g1temp[6],g2temp[6];

        if(*type<1){
          Lstructure::g1N(AZc,var,g1temp,"dssv");
        }else{
          Lstructure::g1N(AZc,var,g1temp,"grsv");
        }

        for(int i=0;i<6;i++){
          g2temp[i]=-g1temp[i];
        }

        // std::cout << "cc:x,g2temp[0] before intg: " << var[0] << " " << g2temp[0] << std:: endl;

        double x0=var[0];
        double dx=(1.0-x0)/nint;

        for(int j=0;j<n;j++){
          var[0]=x0+dx*(double)j;
          if(*type<1){
            Lstructure::g1N(AZc,var,g1temp,"dssv");
          }else{
            Lstructure::g1N(AZc,var,g1temp,"grsv");
          }
          for(int i=0;i<6;i++){
            g2temp[i]=g2temp[i]+dx*g1temp[i]/(x0+dx*(double)j)/2.0;
          }

          var[0]=x0+dx*(double)j+dx;
          if(*type<1){
            Lstructure::g1N(AZc,var,g1temp,"dssv");
          }else{
            Lstructure::g1N(AZc,var,g1temp,"grsv");
          }
          for(int i=0;i<6;i++){
            g2temp[i]=g2temp[i]+dx*g1temp[i]/(x0+dx*(double)j)/2.0;
          }
        }

        // std::cout << "cc: x0,dx,nint,var[0]=" << x0 << " " << dx << " " << nint << " " << var[0] << std::endl;

        // return g2temp[0];

        double result;
        result=0.5*(4.0/9.0)*g2temp[0]+0.5*(1.0/9.0)*g2temp[1]+0.5*(1.0/9.0)*g2temp[2]+0.5*(4.0/9.0)*g2temp[3]+0.5*(1.0/9.0)*g2temp[4]+0.5*(1.0/9.0)*g2temp[5];

        return result;

}
