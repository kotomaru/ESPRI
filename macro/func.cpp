#include <stdio.h>
#include <cmath>
#include "./include/func.h"

void track(double x[2][4],double &th,double &ph,double &x0,double &y0){//BDC1 XAYB, BDC2 XAYB

  const double geom[2][3] = {{4.125,4.125,-142.72},{4.125,4.125,-42.39}}; 
  const double r2d = 180./M_PI;
  double xx[2][4]; //[BDC1 or 2][x,A,y,B]

  for(int i=0;i<2;i++){
    xx[i][0] = x[i][0] - geom[0][0];
    xx[i][1] = x[i][1];
    xx[i][2] = x[i][2] - geom[0][1];
    xx[i][3] = x[i][3];
  }

  th = atan((xx[1][0]-xx[0][0])/(geom[1][2]-geom[0][2]));
  ph = atan((xx[1][2]-xx[0][2])/(geom[1][2]-geom[0][2]));
  th = th*r2d; ph = ph*r2d;
  x0 = (xx[0][0]*geom[1][2]-xx[1][0]*geom[0][2])/(geom[1][2]-geom[0][2]);
  y0 = (xx[0][2]*geom[1][2]-xx[1][2]*geom[0][2])/(geom[1][2]-geom[0][2]);

}

double EdEfunc1(double x){
  const double a=1.5E-4;
  const double b=-0.866;  
  const double c=1697;  
  double y;
  y = a*x*x+b*x+c;
  return y;
}
double EdEfunc2(double x){
  const double a=1.0E-4;
  const double b=-0.866;  
  const double c=1697;  
  double y;
  y = a*x*x+b*x+c;
  return y;
}
double EdEfunc3(double x){
  const double a=9E-5;
  const double b=-0.866;  
  const double c=1897;  
  double y;
  y = a*x*x+b*x+c;
  return y;
}

