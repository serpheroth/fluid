#ifndef __SPH_H_
#define __SPH_H_

#include<cmath>
#include <thread>
#include"water.h"
#include"grid.h"

/* Constants*/

// for particles
const float pdist=0.08;
const int pxnum=10,pynum=10,pznum=10;
const int ptotalNum=pxnum*pynum*pznum;
const float pstartx=-pdist*pxnum/2,pstarty=-pdist*pynum/2,pstartz=-pdist*pznum/2;
const float pradius=20;
const float pmass=1.25e-5;

// for threads
const int threadNum=4;
const int workload=ptotalNum/threadNum;

// global constants
const float gravity=10;
const float tstep=1.0/24;
const float PI=3.1415926;

// for SPH
const float h=0.14;
const float hsqr=h*h;
const float sigma=0.1;
const float wk=315/(64*PI*pow(h,9));
const float dwk=15/(PI*pow(h,6));
const float scorrk=wk*pow(0.99*hsqr,3);
const float prho0=1000;


// for bubbles
const float minta=10;
const float maxta=45;
const float mink=0.1;
const float maxk=0.3;
const int kta=50;
const float maxlifetime=3;
const int SPRAY=0;
const int BUBBLE=1;
const int FOAM=2;
const int LIMIT1=5;
const int LIMIT2=50;
const float kb=0.8;
const float kd=0.1;
const float mass=0.1;

/* functions */
void start();
void update();

#endif
