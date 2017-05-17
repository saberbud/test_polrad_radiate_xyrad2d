#include "g1g2_grsv_ww.h"

int readg1g2grsv(){
  TString filename="/var/phy/project/mepg/xy33/TMD/xy_SIMC/fff_polrad_radiate_test/test_g1g2comb/g1g2table_grsv_ww.txt";
  cout << "g1g2 grsv file: " << filename << endl;
  ifstream infile(filename);
  if(!infile.is_open()){
    cout << "File not exist: " << endl;
    return 1;
  }
  int i=0;
  while(infile >> Q2_grsv[i]>>X_grsv[i]>>g1p_grsv[i]>>g2p_grsv[i]>>g1n_grsv[i]>>g2n_grsv[i]){
    i=i+1;
  }
  cout << "# of lines= " << i << endl;
  infile.close();

  int nline=i;

  double dQ2_grsv[10000],dQ2t_grsv,prod_grsv;
  int idQ2_grsv[200];
  int ndQ2_grsv=0;
  for (int i=0;i<(nline-1);i++){
    dQ2t_grsv=Q2_grsv[i+1]-Q2_grsv[i];
    prod_grsv=1.;
    if(ndQ2_grsv>0){
      if(dQ2_grsv[ndQ2_grsv-1]!=dQ2t_grsv){
	dQ2_grsv[ndQ2_grsv]=dQ2t_grsv;
	idQ2_grsv[ndQ2_grsv]=i;
	ndQ2_grsv=ndQ2_grsv+1;
      }
    }else{
      dQ2_grsv[ndQ2_grsv]=dQ2t_grsv;
      idQ2_grsv[ndQ2_grsv]=i;
      ndQ2_grsv=ndQ2_grsv+1;
    }
  }
  // cout << ndQ2 << endl;

  int diddQ2_grsv[30];
  nddQ2_grsv=0;
  for(int i=0;i<ndQ2_grsv;i++){
    if(nddQ2_grsv>0){
      if((ddQ2_grsv[nddQ2_grsv-1]<dQ2_grsv[i])&&dQ2_grsv[i]>0.){
	// cout << ddQ2[nddQ2-1] << " " << dQ2[i] << endl;
	ddQ2_grsv[nddQ2_grsv]=dQ2_grsv[i];
	iddQ2_grsv[nddQ2_grsv]=idQ2_grsv[i];
	diddQ2_grsv[nddQ2_grsv]=idQ2_grsv[i]-idQ2_grsv[i-1];
	nddQ2_grsv=nddQ2_grsv+1;
      }
    }else{
      ddQ2_grsv[nddQ2_grsv]=dQ2_grsv[i];
      iddQ2_grsv[nddQ2_grsv]=idQ2_grsv[i];
      diddQ2_grsv[nddQ2_grsv]=idQ2_grsv[i];
      nddQ2_grsv=nddQ2_grsv+1;
    }
  }
  iddQ2_grsv[nddQ2_grsv]=idQ2_grsv[ndQ2_grsv-1];
  // cout << nddQ2 << endl;

  return 0;
}

int findg1g2grsv(double Q2in,double xin){
  // cout << "Q2in= " << Q2in << " xin= " << xin << endl;
  g1pr_grsv=0.;
  g2pr_grsv=0.;

  g1nr_grsv=0.;
  g2nr_grsv=0.;

  //Find Q2in
  int nb=0;
  if(Q2in<Q2_grsv[0]||Q2in>Q2_grsv[9599]){
    cout << "Q2in out of range; return" << endl;
    cout << "Q2in= " << Q2in << " range: " << Q2_grsv[0] << " to " << Q2_grsv[9599] << endl;
    return 1; 
  }
  for(int i=0;i<(nddQ2_grsv+1);i++){
    if(Q2in>=Q2_grsv[iddQ2_grsv[i]]&&nb<i){
      nb=i;
    }
  }
  // cout << "Q2in= " << Q2in << " ; " << nb << " " << Q2[iddQ2[nb]] << " " << Q2[iddQ2[nb+1]] << endl;

  int nb1=(Q2in-Q2_grsv[iddQ2_grsv[nb]])/ddQ2_grsv[nb];
  nb1=nb1*128+iddQ2_grsv[nb];
  // cout << "nb1= " << nb1 << " Q2in= " << Q2in << " " << ddQ2[nb] << " " << Q2[iddQ2[nb]] << " " << nb1 << " ; " << Q2[nb1] << " below " << Q2[nb1-1] << " above " << Q2[nb1+128] << endl;
  //end find Q2in

  //Find X
  if(xin<X_grsv[0]||xin>X_grsv[127]){
    cout << "xin out of range: return" << endl;
    cout << "xin= " << xin << " range: " << X_grsv[0] << " to " << X_grsv[127] << endl;
    return 1;
  }

  int nx=0;
  for(int i=0;i<128;i++){
    if(xin<X_grsv[i]&&nx<1){
      nx=i;
    }
  }
  // cout << "nx= " << nx << " ; " << xin << " " << X[nx] << " " << X[nx-1] << " " << X[nx+1] << endl;
  //end find X

  nxy_grsv[0]=nb1+nx;nxy_grsv[1]=nxy_grsv[0]+1;nxy_grsv[2]=nxy_grsv[0]-128;nxy_grsv[3]=nxy_grsv[1]-128;

  double Q2c[5],Xc[5];
  for(int i=0;i<4;i++){
    Q2c[i]=Q2_grsv[nxy_grsv[i]];
    Xc[i]=X_grsv[nxy_grsv[i]];

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
    g1pr_grsv=g1pr_grsv+g1p_grsv[nxy_grsv[i]]*wQ2[i]*wX[i];
    g2pr_grsv=g2pr_grsv+g2p_grsv[nxy_grsv[i]]*wQ2[i]*wX[i];

    g1nr_grsv=g1nr_grsv+g1n_grsv[nxy_grsv[i]]*wQ2[i]*wX[i];
    g2nr_grsv=g2nr_grsv+g2n_grsv[nxy_grsv[i]]*wQ2[i]*wX[i];

    // cout << i << " g1p= " << g1p[nxy[i]] << " wQ2= " << wQ2[i] << " wX= " << wX[i] << endl;
  }
  // cout << "g1pr= " << g1pr << endl;
  // cout << "g2pr= " << g2pr << endl;
  // cout << "g1nr= " << g1nr << endl;
  // cout << "g2nr= " << g2nr << endl;

  return 0;
}

double reg1pgrsv(){
  return g1pr_grsv;
}

double reg2pgrsv(){
  return g2pr_grsv;
}

double reg1ngrsv(){
  return g1nr_grsv;
}

double reg2ngrsv(){
  return g2nr_grsv;
}

