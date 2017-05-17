#include "g1g2_nlo.h"

int readg1g2data(){
  TString filename="/var/phy/project/mepg/xy33/TMD/G1G2NLO/G1G2NLO.dat";
  cout << "g1g2 file: " << filename << endl;
  ifstream infile(filename);
  if(!infile.is_open()){
    cout << "File not exist: " << endl;
    return 1;
  }
  int i=0;
  while(infile >> Q2[i]>>X[i]>>g1p[i]>>g1p_error[i]>>g2p[i]>>g2p_error[i]>>g1n[i]>>g1n_error[i]>>g2n[i]>>g2n_error[i]){
    i=i+1;
  }
  cout << "# of lines= " << i << endl;
  infile.close();

  int nline=i;

  double dQ2[10000],dQ2t,prod;
  int idQ2[200];
  int ndQ2=0;
  for (int i=0;i<(nline-1);i++){
    dQ2t=Q2[i+1]-Q2[i];
    prod=1.;
    if(ndQ2>0){
      if(dQ2[ndQ2-1]!=dQ2t){
	dQ2[ndQ2]=dQ2t;
	idQ2[ndQ2]=i;
	ndQ2=ndQ2+1;
      }
    }else{
      dQ2[ndQ2]=dQ2t;
      idQ2[ndQ2]=i;
      ndQ2=ndQ2+1;
    }
  }
  // cout << ndQ2 << endl;

  int diddQ2[30];
  nddQ2=0;
  for(int i=0;i<ndQ2;i++){
    if(nddQ2>0){
      if((ddQ2[nddQ2-1]<dQ2[i])&&dQ2[i]>0.){
	// cout << ddQ2[nddQ2-1] << " " << dQ2[i] << endl;
	ddQ2[nddQ2]=dQ2[i];
	iddQ2[nddQ2]=idQ2[i];
	diddQ2[nddQ2]=idQ2[i]-idQ2[i-1];
	nddQ2=nddQ2+1;
      }
    }else{
      ddQ2[nddQ2]=dQ2[i];
      iddQ2[nddQ2]=idQ2[i];
      diddQ2[nddQ2]=idQ2[i];
      nddQ2=nddQ2+1;
    }
  }
  iddQ2[nddQ2]=idQ2[ndQ2-1];
  // cout << nddQ2 << endl;

  return 0;
}

int findg1g2(double Q2in,double xin){
  // cout << "Q2in= " << Q2in << " xin= " << xin << endl;
  g1pr=0.;
  g1p_errorr=0.;
  g2pr=0.;
  g2p_errorr=0.;

  g1nr=0.;
  g1n_errorr=0.;
  g2nr=0.;
  g2n_errorr=0.;

  //Find Q2in
  int nb=0;
  if(Q2in<Q2[0]||Q2in>Q2[9599]){
    cout << "Q2in out of range; return" << endl;
    cout << "Q2in= " << Q2in << " range: " << Q2[0] << " to " << Q2[9599] << endl;
    return 1; 
  }
  for(int i=0;i<(nddQ2+1);i++){
    if(Q2in>=Q2[iddQ2[i]]&&nb<i){
      nb=i;
    }
  }
  // cout << "Q2in= " << Q2in << " ; " << nb << " " << Q2[iddQ2[nb]] << " " << Q2[iddQ2[nb+1]] << endl;

  int nb1=(Q2in-Q2[iddQ2[nb]])/ddQ2[nb];
  nb1=nb1*128+iddQ2[nb];
  // cout << "nb1= " << nb1 << " Q2in= " << Q2in << " " << ddQ2[nb] << " " << Q2[iddQ2[nb]] << " " << nb1 << " ; " << Q2[nb1] << " below " << Q2[nb1-1] << " above " << Q2[nb1+128] << endl;
  //end find Q2in

  //Find X
  if(xin<X[0]||xin>X[127]){
    cout << "xin out of range: return" << endl;
    cout << "xin= " << xin << " range: " << X[0] << " to " << X[127] << endl;
    return 1;
  }

  int nx=0;
  for(int i=0;i<128;i++){
    if(xin<X[i]&&nx<1){
      nx=i;
    }
  }
  // cout << "nx= " << nx << " ; " << xin << " " << X[nx] << " " << X[nx-1] << " " << X[nx+1] << endl;
  //end find X

  nxy[0]=nb1+nx;nxy[1]=nxy[0]+1;nxy[2]=nxy[0]-128;nxy[3]=nxy[1]-128;

  double Q2c[5],Xc[5];
  for(int i=0;i<4;i++){
    Q2c[i]=Q2[nxy[i]];
    Xc[i]=X[nxy[i]];

    // cout << "Q2in= " << Q2in << " Q2 bin= " << Q2c[i] << " ; xin= " << xin << " x bin= " << Xc[i] << endl;
  }

  double wQ2[5],wX[5];
  wQ2[0]=1.0-(Q2c[0]-Q2in)/(Q2c[0]-Q2c[2]);
  wQ2[1]=1.0-(Q2c[1]-Q2in)/(Q2c[1]-Q2c[3]);
  wQ2[2]=1.0-(Q2in-Q2c[2])/(Q2c[0]-Q2c[2]);
  wQ2[3]=1.0-(Q2in-Q2c[3])/(Q2c[1]-Q2c[3]);

  wX[0]=1.0-(xin-Xc[0])/(Xc[1]-Xc[0]);
  wX[1]=1.0-(Xc[1]-xin)/(Xc[1]-Xc[0]);
  wX[2]=1.0-(xin-Xc[2])/(Xc[3]-Xc[2]);
  wX[3]=1.0-(Xc[3]-xin)/(Xc[3]-Xc[2]);

  // double wprod=0.;
  // for(int i=0;i<4;i++){
  //   wprod=wprod+wQ2[i]*wX[i];
  //   cout << i << " wQ2= " << wQ2[i] << " wX= " << wX[i] << endl;
  // }
  // cout << "wprod= " << wprod << endl;

  for(int i=0;i<4;i++){
    g1pr=g1pr+g1p[nxy[i]]*wQ2[i]*wX[i];
    g1p_errorr=g1p_errorr+g1p_error[nxy[i]]*wQ2[i]*wX[i];
    g2pr=g2pr+g2p[nxy[i]]*wQ2[i]*wX[i];
    g2p_errorr=g2p_errorr+g2p_error[nxy[i]]*wQ2[i]*wX[i];

    g1nr=g1nr+g1n[nxy[i]]*wQ2[i]*wX[i];
    g1n_errorr=g1n_errorr+g1n_error[nxy[i]]*wQ2[i]*wX[i];
    g2nr=g2nr+g2n[nxy[i]]*wQ2[i]*wX[i];
    g2n_errorr=g2n_errorr+g2n_error[nxy[i]]*wQ2[i]*wX[i];

    // cout << i << " g1p= " << g1p[nxy[i]] << " wQ2= " << wQ2[i] << " wX= " << wX[i] << endl;
  }
  // cout << "g1pr= " << g1pr << endl;
  // cout << "g2pr= " << g2pr << endl;
  // cout << "g1nr= " << g1nr << endl;
  // cout << "g2nr= " << g2nr << endl;

  return 0;
}

double reg1p(){
  return g1pr;
}

double reg2p(){
  return g2pr;
}

double reg1n(){
  return g1nr;
}

double reg2n(){
  return g2nr;
}

void g1g2_nlo(){
  int rd=readg1g2data();
  if(rd>0){
    cout << "Reading data problem, stop." << endl;
    return;
  }else{
    cout << "Data reading fine." << endl;
  }

  // for(int i=0;i<2;i++){
  //   cout << i << setprecision(9) << " " << Q2[i] << " " << X[i] << " " << g1p[i] << " " << g1p_error[i] << " " << g2p[i] << " " << g2p_error[i] << " " << g1n[i] << " " << g1n_error[i] << " " << g2n[i] << " " << g2n_error[i] << endl;
  // }

  // for(int i=0;i<nddQ2;i++){
  //   cout << i << " " << iddQ2[i] << " " << setprecision(9) << ddQ2[i] << " ; " << Q2[iddQ2[i]] << " " << (Q2[iddQ2[i]+1]-Q2[iddQ2[i]]) << endl;
  // }

  double Q2in=2.01;
  double xin=0.155;
  int f=findg1g2(Q2in,xin);

  cout << "Q2in= " << Q2in << " xin= " << xin << endl;

  for(int i=0;i<4;i++){
    cout << i << " Q2c= " << Q2[nxy[i]] << " Xc= " << X[nxy[i]] << " g1p= " << g1p[nxy[i]] << " g1p_error= " << g1p_error[nxy[i]] << " g2p= " << g2p[nxy[i]] << " g2p_error= " << g2p_error[nxy[i]] << endl;
  }
  cout << "g1pr= " << g1pr << " g1p_errorr= " << g1p_errorr << " g2pr= " << g2pr << " g2p_errorr= " << g2p_errorr << endl;
  cout << "g1pr= " << reg1p() << " g2pr= " << reg2p() << endl;

  cout << endl;
  cout << "Q2in= " << Q2in << " xin= " << xin << endl;

  for(int i=0;i<4;i++){
    cout << i << " Q2c= " << Q2[nxy[i]] << " Xc= " << X[nxy[i]] << " g1n= " << g1n[nxy[i]] << " g1n_error= " << g1n_error[nxy[i]] << " g2n= " << g2n[nxy[i]] << " g2n_error= " << g2n_error[nxy[i]] << endl;
  }
  cout << "g1nr= " << g1nr << " g1n_errorr= " << g1n_errorr << " g2nr= " << g2nr << " g2n_errorr= " << g2n_errorr << endl;
  cout << "g1nr= " << reg1n() << " g2nr= " << reg2n() << endl;

}
