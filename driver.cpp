/*
SOLUTION OF NAVIER STOKES EQUATION IN POLAR COORDINATES (R-TH)
CONSERVATIVE FINITE DIFFERENCE DISCRETIZATION TECHNIQUE IN SEMI-COLOCATED GRID
USING BALANCED FORCE PROJECTION ALGORITHM
THE ADVECTION SCHEME IS QUICK + FIRST ORDER UPWIND
FDQGMRES OR BI-CGSTAB CAN BE USED TO SOLVE THE LINEAR SYSTEMS
*/
#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<csignal>
#define EPS 1e-6	//DO NOT CHANGE THESE TWO CONSTANTS
#define SMALL 1e-8
#define TOL 1e-6	//convergence criteria for FDQGMRES solver
using namespace std;
int LC=1;	//loop condition
void control(int n)	//signal control function
{
	cout<<"PROGRAM INTERUPTED!"<<endl;
	cout<<"Enter 0 to stop: ";
	cin>>LC;
	if(LC==0) cout<<"STOP SIGNAL RECIEVED. PROGRAM WILL STOP AFTER SIMULATING CURRENT TIME-STEP."<<endl;
	else {LC=1; cout<<"SIMULATION CONTINUES."<<endl;}
}
const int I=256,J=256;	//no of grids in i and j directions
const double R_in=1.0,R_out=2.0;	//domain size in radial direction
const double TH_s=-0.5,TH_e=0.5;	//domain size in azimuthal direction (in radian)
#include "cfd_solvers.h"
#include "mbase.cpp"
#include "GMG2.cpp"
#include "MG_FDQGMRES.cpp"
//#include "MG_BICGSTAB.cpp"
#include "NS_g.cpp"
int main()
{
	signal(SIGINT,control);	//define signal and its function (FOR MANUALLY CONTROLLING THE PROGRAM DURING RUNTIME)
	int CNT=0;
	NS ms;
	ms.MBASE::ini(CNT,350.0,1.0,8e-4);	//count,rho,mu,dt
	ms.grid_write();
	for(int COUNT=CNT;(LC && (COUNT<100000));COUNT++)	//manually controlled loop
	{
		ms.NS::solve();
		if((((COUNT+1)%100)==0)||(COUNT==0))
		{
			ms.max_CFL();
			ms.write_bin(COUNT+1);
			//ms.write(COUNT+1);
			//ms.write_s(COUNT+1);
			//ms.write_adv(COUNT+1);
		}
	}
	return 0;
}
