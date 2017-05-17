#include "g1g2_dssv_ww.h"

int readg1g2dssv(){
  TString filename="/var/phy/project/mepg/xy33/TMD/xy_SIMC/fff_polrad_radiate_test/test_g1g2comb/g1g2table_dssv_ww.txt";
  cout << "g1g2 dssv file: " << filename << endl;
  ifstream infile(filename);
  if(!infile.is_open()){
    cout << "File not exist: " << endl;
    return 1;
  }
  int i=0;
  while(infile >> Q2_dssv[i]>>X_dssv[i]>>g1p_dssv[i]>>g2p_dssv[i]>>g1n_dssv[i]>>g2n_dssv[i]){
    i=i+1;
  }
  cout << "# of lines= " << i << endl;
  infile.close();

  int nline=i;

  double dQ2_dssv[10000],dQ2t_dssv,prod_dssv;
  int idQ2_dssv[200];
  int ndQ2_dssv=0;
  for (int i=0;i<(nline-1);i++){
    dQ2t_dssv=Q2_dssv[i+1]-Q2_dssv[i];
    prod_dssv=1.;
    if(ndQ2_dssv>0){
      if(dQ2_dssv[ndQ2_dssv-1]!=dQ2t_dssv){
	dQ2_dssv[ndQ2_dssv]=dQ2t_dssv;
	idQ2_dssv[ndQ2_dssv]=i;
	ndQ2_dssv=ndQ2_dssv+1;
      }
    }else{
      dQ2_dssv[ndQ2_dssv]=dQ2t_dssv;
      idQ2_dssv[ndQ2_dssv]=i;
      ndQ2_dssv=ndQ2_dssv+1;
    }
  }
  // cout << ndQ2 << endl;

  int diddQ2_dssv[30];
  nddQ2_dssv=0;
  for(int i=0;i<ndQ2_dssv;i++){
    if(nddQ2_dssv>0){
      if((ddQ2_dssv[nddQ2_dssv-1]<dQ2_dssv[i])&&dQ2_dssv[i]>0.){
	// cout << ddQ2[nddQ2-1] << " " << dQ2[i] << endl;
	ddQ2_dssv[nddQ2_dssv]=dQ2_dssv[i];
	iddQ2_dssv[nddQ2_dssv]=idQ2_dssv[i];
	diddQ2_dssv[nddQ2_dssv]=idQ2_dssv[i]-idQ2_dssv[i-1];
	nddQ2_dssv=nddQ2_dssv+1;
      }
    }else{
      ddQ2_dssv[nddQ2_dssv]=dQ2_dssv[i];
      iddQ2_dssv[nddQ2_dssv]=idQ2_dssv[i];
      diddQ2_dssv[nddQ2_dssv]=idQ2_dssv[i];
      nddQ2_dssv=nddQ2_dssv+1;
    }
  }
  iddQ2_dssv[nddQ2_dssv]=idQ2_dssv[ndQ2_dssv-1];
  // cout << nddQ2 << endl;

  return 0;
}

int findg1g2dssv(double Q2in,double xin){
  // cout << "Q2in= " << Q2in << " xin= " << xin << endl;
  g1pr_dssv=0.;
  g2pr_dssv=0.;

  g1nr_dssv=0.;
  g2nr_dssv=0.;

  //Find Q2in
  int nb=0;
  if(Q2in<Q2_dssv[0]||Q2in>Q2_dssv[9599]){
    cout << "Q2in out of range; return" << endl;
    cout << "Q2in= " << Q2in << " range: " << Q2_dssv[0] << " to " << Q2_dssv[9599] << endl;
    return 1; 
  }
  for(int i=0;i<(nddQ2_dssv+1);i++){
    if(Q2in>=Q2_dssv[iddQ2_dssv[i]]&&nb<i){
      nb=i;
    }
  }
  // cout << "Q2in= " << Q2in << " ; " << nb << " " << Q2[iddQ2[nb]] << " " << Q2[iddQ2[nb+1]] << endl;

  int nb1=(Q2in-Q2_dssv[iddQ2_dssv[nb]])/ddQ2_dssv[nb];
  nb1=nb1*128+iddQ2_dssv[nb];
  // cout << "nb1= " << nb1 << " Q2in= " << Q2in << " " << ddQ2[nb] << " " << Q2[iddQ2[nb]] << " " << nb1 << " ; " << Q2[nb1] << " below " << Q2[nb1-1] << " above " << Q2[nb1+128] << endl;
  //end find Q2in

  //Find X
  if(xin<X_dssv[0]||xin>X_dssv[127]){
    cout << "xin out of range: return" << endl;
    cout << "xin= " << xin << " range: " << X_dssv[0] << " to " << X_dssv[127] << endl;
    return 1;
  }

  int nx=0;
  for(int i=0;i<128;i++){
    if(xin<X_dssv[i]&&nx<1){
      nx=i;
    }
  }
  // cout << "nx= " << nx << " ; " << xin << " " << X[nx] << " " << X[nx-1] << " " << X[nx+1] << endl;
  //end find X

  nxy_dssv[0]=nb1+nx;nxy_dssv[1]=nxy_dssv[0]+1;nxy_dssv[2]=nxy_dssv[0]-128;nxy_dssv[3]=nxy_dssv[1]-128;

  double Q2c[5],Xc[5];
  for(int i=0;i<4;i++){
    Q2c[i]=Q2_dssv[nxy_dssv[i]];
    Xc[i]=X_dssv[nxy_dssv[i]];

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
    g1pr_dssv=g1pr_dssv+g1p_dssv[nxy_dssv[i]]*wQ2[i]*wX[i];
    g2pr_dssv=g2pr_dssv+g2p_dssv[nxy_dssv[i]]*wQ2[i]*wX[i];

    g1nr_dssv=g1nr_dssv+g1n_dssv[nxy_dssv[i]]*wQ2[i]*wX[i];
    g2nr_dssv=g2nr_dssv+g2n_dssv[nxy_dssv[i]]*wQ2[i]*wX[i];

    // cout << i << " g1p= " << g1p[nxy[i]] << " wQ2= " << wQ2[i] << " wX= " << wX[i] << endl;
  }
  // cout << "g1pr= " << g1pr << endl;
  // cout << "g2pr= " << g2pr << endl;
  // cout << "g1nr= " << g1nr << endl;
  // cout << "g2nr= " << g2nr << endl;

  return 0;
}

double reg1pdssv(){
  return g1pr_dssv;
}

double reg2pdssv(){
  return g2pr_dssv;
}

double reg1ndssv(){
  return g1nr_dssv;
}

double reg2ndssv(){
  return g2nr_dssv;
}

