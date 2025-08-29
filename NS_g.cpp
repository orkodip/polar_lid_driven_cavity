/*
NAVIER-STOKES SOLVER IN POLAR COORDINATES
CONVENTION USED ARE AS FOLLOWS
(R,TH)->(X,Y)->(i,j)->(I,J)
V=(V_r,V_th)->(u,v)
BOUNDARY CONDITIONS ARE NO SLIP AND NO PENETRATION

R VELOCITY BC -> DIRICHLET

TH VELOCITY BC -> DIRICHLET
LEFT BOUNDARY -> MOVING LID

PRESSURE BC -> WALL BC (PURE NEUMANN)
 */
class NS:public MBASE,public MG_FDQGMRES
//class NS:public MBASE,public MG_BICGSTAB
{
	protected:
	double **us,**vs;	//cell centered intermediate velocity field
	double **u_adv_EW,**v_adv_EW,**u_adv_NS,**v_adv_NS;	//interpolated velocities at the cell faces using advection scheme
	double **D_x,**D_y;	//cell centered diffusion term for R and TH momentum equations

	void vel_bc();	//velocity boundary condition(ghost nodes are handled here)
	void P_bc();	//pressure boundary condition(ghost nodes are handled here)
	void advection();	//explicit advection term
	void diff();	//explicit viscous term
	void mom();	//solve the inviscid momentum equations
	void calc_fc();	//calculate face centered velocities
	void Press();	//solve the pressure Poisson equation
	void update();	//velocity and pressure update
	void continuity();	//calculate the continuity
	public:
			NS(); ~NS();
			void solve();
			void max_CFL();
			void write_s(int t);
			void write_bin(int count);	//export intermediate file in binary format
			void read_bin(string fname);	//import intermediate data
};
NS::NS():MG_FDQGMRES(8,I*J,100)	//FDQGMRES parameterized constructor
//NS::NS():MG_BICGSTAB(I*J,100)	//BICGSTAB parameterized constructor
{
	us=new double*[J+2];
	vs=new double*[J+2];
	u_adv_EW=new double*[J+1];
	u_adv_NS=new double*[J+1];
	v_adv_EW=new double*[J+1];
	v_adv_NS=new double*[J+1];
	D_x=new double*[J+1];
	D_y=new double*[J+1];
	for(int j=0;j<J+2;j++)
	{
		us[j]=new double[I+2];
		vs[j]=new double[I+2];
		if(j<J+1)
		{
			u_adv_EW[j]=new double[I+1];
			u_adv_NS[j]=new double[I+1];
			v_adv_EW[j]=new double[I+1];
			v_adv_NS[j]=new double[I+1];
			D_x[j]=new double[I+1];
			D_y[j]=new double[I+1];
		}
	}
	cout<<"NS: MEMORY ALLOCATED"<<endl;
}
NS::~NS()
{
	for(int j=0;j<J+2;j++)
	{
		delete[] us[j];
		delete[] vs[j];
		if(j<J+1)
		{
			delete[] u_adv_EW[j];
			delete[] u_adv_NS[j];
			delete[] v_adv_EW[j];
			delete[] v_adv_NS[j];
			delete[] D_x[j];
			delete[] D_y[j];
		}
	}
	delete[] us;
	delete[] vs;
	delete[] u_adv_EW;
	delete[] u_adv_NS;
	delete[] v_adv_EW;
	delete[] v_adv_NS;
	delete[] D_x;
	delete[] D_y;
	cout<<"NS: MEMORY RELEASED"<<endl;
}
void NS::vel_bc()
{
	for(int i=1;i<=I;i++)	//no slip and no penetration
	{
		v_NS[0][i]=0.0;	//bottom boundary
		u[0][i]=-u[1][i];
		v[0][i]=-v[1][i];
		u_adv_NS[0][i]=0.0;
		v_adv_NS[0][i]=0.0;
		
		v_NS[J][i]=0.0;	//top boundary
		u[J+1][i]=-u[J][i];
		v[J+1][i]=-v[J][i];
		u_adv_NS[J][i]=0.0;
		v_adv_NS[J][i]=0.0;
	}
	for(int j=1;j<=J;j++)	//no slip and no penetration
	{
		u_EW[j][0]=0.0;	//left boundary
		u[j][0]=-u[j][1];
		v[j][0]=2.0-v[j][1];
		u_adv_EW[j][0]=0.0;
		v_adv_EW[j][0]=1.0;
		
		u_EW[j][I]=0.0;	//right boundary
		u[j][I+1]=-u[j][I];
		v[j][I+1]=-v[j][I];
		u_adv_EW[j][I]=0.0;
		v_adv_EW[j][I]=0.0;
	}
	u[0][0]=-u[1][0]; u[J+1][0]=-u[J][0]; u[0][I+1]=-u[1][I+1]; u[J+1][I+1]=-u[J][I+1];
	v[0][0]=-v[1][0]; v[J+1][0]=-v[J][0]; v[0][I+1]=-v[1][I+1]; v[J+1][I+1]=-v[J][I+1];
}
void NS::P_bc()
{
	for(int j=1;j<=J;j++)	//Neumann
	{
		P[j][0]=P[j][1];
		P[j][I+1]=P[j][I];
	}
	for(int i=1;i<=I;i++)	//Neumann
	{
		P[0][i]=P[1][i];
		P[J+1][i]=P[J][i];
	}
}
void NS::advection()	//QUICK scheme is used
{
	for(int j=1;j<=J;j++)	//calculation in inner domain
	{
		for(int i=1;i<=I;i++)
		{
			if((i>1)&&(i<I-1))	//QUICK scheme for EW faces
			{
				u_adv_EW[j][i]=(u_EW[j][i]>0.0) ? (0.75*u[j][i]+0.375*u[j][i+1]-0.125*u[j][i-1])
						: (0.75*u[j][i+1]+0.375*u[j][i]-0.125*u[j][i+2]);
				v_adv_EW[j][i]=(u_EW[j][i]>0.0) ? (0.75*v[j][i]+0.375*v[j][i+1]-0.125*v[j][i-1])
						: (0.75*v[j][i+1]+0.375*v[j][i]-0.125*v[j][i+2]);
			}
			else if((i==1)||(i==(I-1)))	//FOU scheme for EW faces
			{
				u_adv_EW[j][i]=(u_EW[j][i]>0.0) ? u[j][i] : u[j][i+1];
				v_adv_EW[j][i]=(u_EW[j][i]>0.0) ? v[j][i] : v[j][i+1];
			}
			if((j>1)&&(j<J-1))	//QUICK scheme for NS faces
			{
				u_adv_NS[j][i]=(v_NS[j][i]>0.0) ? (0.75*u[j][i]+0.375*u[j+1][i]-0.125*u[j-1][i])
						: (0.75*u[j+1][i]+0.375*u[j][i]-0.125*u[j+2][i]);
				v_adv_NS[j][i]=(v_NS[j][i]>0.0) ? (0.75*v[j][i]+0.375*v[j+1][i]-0.125*v[j-1][i])
						: (0.75*v[j+1][i]+0.375*v[j][i]-0.125*v[j+2][i]);
			}
			else if((j==1)||(j==(J-1)))	//FOU scheme for NS faces
			{
				u_adv_NS[j][i]=(v_NS[j][i]>0.0) ? u[j][i] : u[j+1][i];
				v_adv_NS[j][i]=(v_NS[j][i]>0.0) ? v[j][i] : v[j+1][i];
			}
		}
	}
	for(int j=1;j<=J;j++)	//calculate advection value
	{
		for(int i=1;i<=I;i++)
		{
			A_x[j][i]=rho*((Xm[i]*u_EW[j][i]*u_adv_EW[j][i]-Xm[i-1]*u_EW[j][i-1]*u_adv_EW[j][i-1])/(CX[i]*dx)
					+(v_NS[j][i]*u_adv_NS[j][i]-v_NS[j-1][i]*u_adv_NS[j-1][i])/(CX[i]*dy));
			A_y[j][i]=rho*((Xm[i]*u_EW[j][i]*v_adv_EW[j][i]-Xm[i-1]*u_EW[j][i-1]*v_adv_EW[j][i-1])/(CX[i]*dx)
					+(v_NS[j][i]*v_adv_NS[j][i]-v_NS[j-1][i]*v_adv_NS[j-1][i])/(CX[i]*dy));
		}
	}
}
void NS::diff()
{
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			D_x[j][i]=((Xm[i]*2.0*mu*(u[j][i+1]-u[j][i])/dx)
					-(Xm[i-1]*2.0*mu*(u[j][i]-u[j][i-1])/dx))/(CX[i]*dx)
					+(mu*(0.25*(v[j][i+1]+v[j+1][i+1]-v[j][i-1]-v[j+1][i-1])/dx-0.5*(v[j+1][i]+v[j][i])/CX[i]+(u[j+1][i]-u[j][i])/(CX[i]*dy))
					-mu*(0.25*(v[j][i+1]+v[j-1][i+1]-v[j][i-1]-v[j-1][i-1])/dx-0.5*(v[j][i]+v[j-1][i])/CX[i]+(u[j][i]-u[j-1][i])/(CX[i]*dy)))/(CX[i]*dy);
			D_y[j][i]=(Xm[i]*mu*((v[j][i+1]-v[j][i])/dx-0.5*(v[j][i]+v[j][i+1])/Xm[i]+0.25*(u[j+1][i]+u[j+1][i+1]-u[j-1][i]-u[j-1][i+1])/(Xm[i]*dy))
					-Xm[i-1]*mu*((v[j][i]-v[j][i-1])/dx-0.5*(v[j][i]+v[j][i-1])/Xm[i-1]+0.25*(u[j+1][i]+u[j+1][i-1]-u[j-1][i]-u[j-1][i-1])/(Xm[i-1]*dy)))/(CX[i]*dx)
					+(2.0*mu*((v[j+1][i]-v[j][i])/(CX[i]*dy)+0.5*(u[j+1][i]+u[j][i])/CX[i])
					-2.0*mu*((v[j][i]-v[j-1][i])/(CX[i]*dy)+0.5*(u[j][i]+u[j-1][i])/CX[i]))/(CX[i]*dy);
		}
	}
}
void NS::mom()
{
	advection();	//calculate advection term
	diff();	//calculate diffusion term
	for(int j=1;j<=J;j++)	//explicit calculation
	{
		for(int i=1;i<=I;i++)
		{
			us[j][i]=u[j][i]+dt*(-A_x[j][i]+D_x[j][i]+rho*v[j][i]*v[j][i]/CX[i]
					-(2.0*mu*(0.5*(v[j+1][i]-v[j-1][i])/(CX[i]*dy)+u[j][i]/CX[i]))/CX[i])/rho;
			vs[j][i]=v[j][i]+dt*(-A_y[j][i]+D_y[j][i]-rho*u[j][i]*v[j][i]/CX[i]
					+mu*(0.5*(v[j][i+1]-v[j][i-1])/dx-v[j][i]/CX[i]+0.5*(u[j+1][i]-u[j-1][i])/(CX[i]*dy)))/rho;
		}
	}
}
void NS::calc_fc()
{
	for(int j=1;j<=J;j++)	//calculation of face centered values excluding the boundaries
	{
		for(int i=1;i<=I;i++)
		{
			if(i<I) u_EW[j][i]=0.5*(us[j][i]+us[j][i+1]);
			if(j<J) v_NS[j][i]=0.5*(vs[j][i]+vs[j+1][i]);
		}
	}
}
void NS::Press()
{
	if(COUNT==0)	//setup the equations in the first iteration
	{
		int cnt=0; R0[cnt]=0;
		double alpha,beta;
		for(int j=1;j<=J;j++)
		{
			for(int i=1;i<=I;i++)
			{
				b0[(j-1)*I+i-1]=rho*dx*dx/(CX[i]*dt)*((Xm[i]*u_EW[j][i]-Xm[i-1]*u_EW[j][i-1])/dx+(v_NS[j][i]-v_NS[j-1][i])/dy);
				alpha=0.5*dx/CX[i]; beta=pow((dx/(CX[i]*dy)),2.0);
				if((i==1)&&(j==1))	//bottom left corner cell
				{
					C0[cnt]=(j-1)*I+i-1;	//1,1 term
					A0[cnt]=-(1.0+alpha+beta);
					cnt++;
					C0[cnt]=(j-1)*I+i;	//2,1 term
					A0[cnt]=1.0+alpha;
					cnt++;
					C0[cnt]=j*I+i-1;	//1,2 term
					A0[cnt]=beta;
					cnt++;
				}
				else if((i==1)&&(j>1)&&(j<J))	//left boundary
				{
					C0[cnt]=(j-2)*I+i-1;	//1,j-1 term
					A0[cnt]=beta;
					cnt++;
					C0[cnt]=(j-1)*I+i-1;	//1,j term
					A0[cnt]=-(1.0+alpha+2.0*beta);
					cnt++;
					C0[cnt]=(j-1)*I+i;	//2,j term
					A0[cnt]=1.0+alpha;
					cnt++;
					C0[cnt]=j*I+i-1;	//1,j+1 term
					A0[cnt]=beta;
					cnt++;
				}
				else if((i==1)&&(j==J))	//top left corner cell
				{
					C0[cnt]=(j-2)*I+i-1;	//1,J-1 term
					A0[cnt]=beta;
					cnt++;
					C0[cnt]=(j-1)*I+i-1;	//1,J term
					A0[cnt]=-(1.0+alpha+beta);
					cnt++;
					C0[cnt]=(j-1)*I+i;	//2,J term
					A0[cnt]=1.0+alpha;
					cnt++;
				}
				else if((j==J)&&(i>1)&&(i<I))	//top boundary
				{	
					C0[cnt]=(j-2)*I+i-1;	//i,J-1 term
					A0[cnt]=beta;
					cnt++;
					C0[cnt]=(j-1)*I+i-2;	//i-1,J term
					A0[cnt]=1.0-alpha;
					cnt++;
					C0[cnt]=(j-1)*I+i-1;	//i,J term
					A0[cnt]=-(2.0+beta);
					cnt++;
					C0[cnt]=(j-1)*I+i;	//i+1,J term
					A0[cnt]=1.0+alpha;
					cnt++;
				}
				else if((j==J)&&(i==I))	//top right corner cell
				{
					C0[cnt]=(j-2)*I+i-1;	//I,J-1 term
					A0[cnt]=beta;
					cnt++;
					C0[cnt]=(j-1)*I+i-2;	//I-1,J term
					A0[cnt]=1.0-alpha;
					cnt++;
					C0[cnt]=(j-1)*I+i-1;	//I,J term
					A0[cnt]=-(1.0-alpha+beta);
					cnt++;
				}
				else if((i==I)&&(j>1)&&(j<J))	//right boundary
				{
					C0[cnt]=(j-2)*I+i-1;	//I,j-1 term
					A0[cnt]=beta;
					cnt++;
					C0[cnt]=(j-1)*I+i-2;	//I-1,j term
					A0[cnt]=1.0-alpha;
					cnt++;
					C0[cnt]=(j-1)*I+i-1;	//I,j term
					A0[cnt]=-(1.0-alpha+2.0*beta);
					cnt++;
					C0[cnt]=j*I+i-1;	//I,j+1 term
					A0[cnt]=beta;
					cnt++;
				}
				else if((i==I)&&(j==1))	//bottom right corner cell
				{
					C0[cnt]=(j-1)*I+i-2;	//I-1,1 term
					A0[cnt]=1.0-alpha;
					cnt++;
					C0[cnt]=(j-1)*I+i-1;	//I,1 term
					A0[cnt]=-(1.0-alpha+beta);
					cnt++;
					C0[cnt]=j*I+i-1;	//I,2 term
					A0[cnt]=beta;
					cnt++;
				}
				else if((j==1)&&(i>1)&&(i<I))	//bottom boundary
				{
					C0[cnt]=(j-1)*I+i-2;	//i-1,1 term
					A0[cnt]=1.0-alpha;
					cnt++;
					C0[cnt]=(j-1)*I+i-1;	//i,1 term
					A0[cnt]=-(2.0+beta);
					cnt++;
					C0[cnt]=(j-1)*I+i;	//i+1,1 term
					A0[cnt]=1.0+alpha;
					cnt++;
					C0[cnt]=j*I+i-1;	//i,2 term
					A0[cnt]=beta;
					cnt++;
				}
				else	//inside domain
				{
					C0[cnt]=(j-2)*I+i-1;	//i,j-1 term
					A0[cnt]=beta;
					cnt++;
					C0[cnt]=(j-1)*I+i-2;	//i-1,j term
					A0[cnt]=1.0-alpha;
					cnt++;
					C0[cnt]=(j-1)*I+i-1;	//i,j term
					A0[cnt]=-2.0*(1.0+beta);
					cnt++;
					C0[cnt]=(j-1)*I+i;	//i+1,j term
					A0[cnt]=1.0+alpha;
					cnt++;
					C0[cnt]=j*I+i-1;	//i,j+1 term
					A0[cnt]=beta;
					cnt++;
				}
				R0[(j-1)*I+i]=cnt;
			}
		}
		GMG2::setup(I/2,J/2,C0,R0,A0,C1,R1,A1);	//create the coarse grid coefficient matrices
		GMG2::setup(I/4,J/4,C1,R1,A1,C2,R2,A2);
		GMG2::setup(I/8,J/8,C2,R2,A2,C3,R3,A3);
		GMG2::setup(I/16,J/16,C3,R3,A3,C4,R4,A4);
	}
	else	//create only RHS vector for the subsequent iterations
	{
		for(int j=1;j<=J;j++)
			for(int i=1;i<=I;i++)
				b0[(j-1)*I+i-1]=rho*dx*dx/(CX[i]*dt)*((Xm[i]*u_EW[j][i]-Xm[i-1]*u_EW[j][i-1])/dx+(v_NS[j][i]-v_NS[j-1][i])/dy);
	}
	MG_FDQGMRES::solve(C0,R0,A0,X0,b0);
	//MG_BICGSTAB::solve(C0,R0,A0,X0,b0);
	for(int j=1;j<=J;j++)	//updation of the calculated P
		for(int i=1;i<=I;i++)
			P[j][i]=X0[(j-1)*I+i-1];
}
void NS::update()
{
	double rdx=dt/dx,rdy=dt/dy;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			if(i<I) u_EW[j][i]=u_EW[j][i]-rdx*(P[j][i+1]-P[j][i])/rho;
			if(j<J) v_NS[j][i]=v_NS[j][i]-rdy*(P[j+1][i]-P[j][i])/(CX[i]*rho);
			u[j][i]=us[j][i]-0.5*rdx*((P[j][i+1]-P[j][i])/rho+(P[j][i]-P[j][i-1])/rho);
			v[j][i]=vs[j][i]-0.5*rdy*((P[j+1][i]-P[j][i])/rho+(P[j][i]-P[j-1][i])/rho)/CX[i];
		}
	}
}
void NS::continuity()
{
	double cont=0.0;
	for(int j=1;j<=J;j++)
		for(int i=1;i<=I;i++)
			cont+=(Xm[i]*u_EW[j][i]-Xm[i-1]*u_EW[j][i-1])/dx+(v_NS[j][i]-v_NS[j-1][i])/dy;	//inner domain
	cout<<"NS: count = "<<COUNT<<" cont = "<<abs(cont)<<endl;
}
void NS::solve()
{
	vel_bc();
	mom();
	calc_fc();
	Press(); P_bc();
	update();
	COUNT++;
	if((COUNT%10)==0) { continuity(); max_CFL(); }
}
void NS::max_CFL()
{
	double u_max=0.0,v_max=0.0;
	double cfl_x,cfl_y;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			if(abs(u[j][i])>u_max) u_max=abs(u[j][i]);
			if(abs(v[j][i])>v_max) v_max=abs(v[j][i]);
		}
	}
	cfl_x=u_max*dt/dx; cfl_y=v_max*dt/dy;
	if(cfl_x>cfl_y) cout<<"NS: CFL_cell = "<<cfl_x<<endl;
	else cout<<"NS: CFL_cell = "<<cfl_y<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			if(abs(u_EW[j][i])>u_max) u_max=abs(u_EW[j][i]);
			if(abs(v_NS[j][i])>v_max) v_max=abs(v_NS[j][i]);
		}
	}
	cfl_x=u_max*dt/dx; cfl_y=v_max*dt/dy;
	if(cfl_x>cfl_y) cout<<"NS: CFL_f = "<<cfl_x<<endl;
	else cout<<"NS: CFL_f = "<<cfl_y<<endl;
}
void NS::write_s(int t)
{
	string fname="uvs_"+to_string(t)+".dat";
	ofstream p_out(fname);
	p_out<<"TITLE = \"STAR FIELD FLOW\""<<endl;
	p_out<<"FILETYPE = SOLUTION"<<endl;
	p_out<<"VARIABLES = \"us\",\"vs\""<<endl;
	p_out<<"ZONE T=\""<<t*dt<<"\", I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK, VARLOCATION=([1,2]=CELLCENTERED), SOLUTIONTIME="<<t*dt<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<us[j][i];
		p_out<<endl;
	}
	p_out<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<vs[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"NS: STAR FIELD FILE OUTPUT SUCCESSFULL AT n = "<<t<<endl;
}
void NS::write_bin(int count)
{
	string fname="inter_"+to_string(count);
	ofstream p_out(fname);
	MBASE::write_bin(p_out);
	p_out.close();
	cout<<"NS: INTERMEDIATE FILE OUTPUT SUCCESSFULL AT n = "<<count<<endl;
}
void NS::read_bin(string fname)
{
	ifstream p_in(fname);
	if(p_in.fail()) throw(0);	//file does not exist
	MBASE::read_bin(p_in);
	p_in.close();
	cout<<"NS: SOLUTION INITIALIZED SUCCESSFULLY"<<endl;
}
