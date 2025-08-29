/*
ABSTRACT BASE CLASS
CONVENTION USED ARE AS FOLLOWS
(R,TH)->(X,Y)->(i,j)->(I,J)
V=(V_r,V_th)->(u,v)
*/
class MBASE
{
	protected:
	double *Xm, *Ym;	//the mesh variables
	double *CX, *CY;	//the node variables
	double **u,**v,**P;	//cell centered R velocity, TH velocity, pressure
	double **u_EW,**v_NS;	//advection velocities at the cell faces
	double **A_x,**A_y;	//cell centered advection term for R and TH momentum equations
	double rho,mu;	//density and viscosity
	double dx,dy,dt;	//grid spacing, time step
	int COUNT;
	void grid_gen();	//uniform grid generator
	public:
			MBASE(); ~MBASE();
			void ini(int count,double Rho,double Mu,double Dt);	//initialize the problem
			void grid_write();	//export grid for Tecplot360
			void write(int t);	//export velocity and pressure field for Tecplot360
			void write_adv(int t);	//export advection field
			void write_bin(ofstream &p_out);	//export intermediate data in binary format
			void read_bin(ifstream &p_in);	//import intermediate data
};
MBASE::MBASE()
{
	Xm=new double[I+1];
	Ym=new double[J+1];
	CX=new double[I+1];
	CY=new double[J+1];
	P=new double*[J+2];
	u_EW=new double*[J+2];
	u=new double*[J+2];
	v=new double*[J+2];
	v_NS=new double*[J+1];
	A_x=new double*[J+1];
	A_y=new double*[J+1];
	for(int j=0;j<J+2;j++)
	{
		P[j]=new double[I+2];
		u[j]=new double[I+2];
		v[j]=new double[I+2];
		u_EW[j]=new double[I+1];
		if(j<J+1)
		{
			v_NS[j]=new double[I+2];
			A_x[j]=new double[I+1];
			A_y[j]=new double[I+1];
		}
	}
	cout<<"MBASE: MEMORY ALLOCATED"<<endl;
}
MBASE::~MBASE()
{
	for(int j=0;j<J+2;j++)
	{
		delete[] P[j];
		delete[] u[j];
		delete[] v[j];
		delete[] u_EW[j];
		if(j<J+1)
		{
			delete[] v_NS[j];
			delete[] A_x[j];
			delete[] A_y[j];
		}
	}
	delete[] Xm;
	delete[] Ym;
	delete[] CX;
	delete[] CY;
	delete[] P;
	delete[] u;
	delete[] v;
	delete[] u_EW;
	delete[] v_NS;
	delete[] A_x;
	delete[] A_y;
	cout<<"MBASE: MEMORY RELEASED"<<endl;
}
void MBASE::ini(int count,double Rho,double Mu,double Dt)
{
	rho=Rho; mu=Mu;
	dx=(R_out-R_in)/I; dy=(TH_e-TH_s)/J;
	dt=Dt;
	COUNT=count;
	cout<<"MBASE: Rho = "<<rho<<", Mu = "<<mu<<endl;
	cout<<"MBASE: dt = "<<dt<<endl;
	grid_gen();
}
void MBASE::grid_gen()	//uniform grid generation done here
{
	Xm[0]=R_in;
	CX[0]=Xm[0]-0.5*dx;
	Ym[0]=TH_s;
	CY[0]=Ym[0]-0.5*dy;
	for(int i=1;i<=I;i++)	//generation of mesh and cell centers
	{
		Xm[i]=Xm[i-1]+dx;
		CX[i]=Xm[i]-0.5*dx;
	}
	for(int j=1;j<=J;j++)
	{
		Ym[j]=Ym[j-1]+dy;
		CY[j]=Ym[j]-0.5*dy;
	}
	cout<<"MBASE: GRID GENERATED"<<endl;
}
void MBASE::grid_write()
{
	ofstream p_out("mesh.dat");
	p_out<<"TITLE = \"MESH\""<<endl;
	p_out<<"FILETYPE = GRID"<<endl;
	p_out<<"VARIABLES = \"X\",\"Y\""<<endl;
	p_out<<"ZONE I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK"<<endl;
	for(int j=0;j<=J;j++)	//print X co-ordinates of mesh
	{
		for(int i=0;i<=I;i++)
			p_out<<" "<<Xm[i]*cos(Ym[j]);
		p_out<<endl;
	}
	p_out<<endl<<endl;
	for(int j=0;j<=J;j++)	//print Y co-ordinates of mesh
	{
		for(int i=0;i<=I;i++)
			p_out<<" "<<Xm[i]*sin(Ym[j]);
		p_out<<endl;
	}
	p_out.close();
	cout<<"MBASE: GRID WRITE SUCCESSFULL"<<endl;
}
void MBASE::write(int t)
{
	string fname="uvp_"+to_string(t)+".dat";
	ofstream p_out(fname);
	p_out<<"TITLE = \"FLOW AND PRESSURE FIELD\""<<endl;
	p_out<<"FILETYPE = SOLUTION"<<endl;
	p_out<<"VARIABLES = \"u\",\"v\",\"P\""<<endl;
	p_out<<"ZONE T=\""<<t*dt<<"\", I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK, VARLOCATION=([1,2,3]=CELLCENTERED), SOLUTIONTIME="<<t*dt<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<u[j][i];
		p_out<<endl;
	}
	p_out<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<v[j][i];
		p_out<<endl;
	}
	p_out<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<P[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"MBASE: FILE OUTPUT SUCCESSFULL AT n = "<<t<<endl;
}
void MBASE::write_adv(int t)
{
	string fname="adv_"+to_string(t)+".dat";
	ofstream p_out(fname);
	p_out<<"TITLE = \"ADVECTION FIELD\""<<endl;
	p_out<<"FILETYPE = SOLUTION"<<endl;
	p_out<<"VARIABLES = \"adv_x\",\"adv_y\""<<endl;
	p_out<<"ZONE T=\""<<t*dt<<"\", I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK, VARLOCATION=([1,2]=CELLCENTERED), SOLUTIONTIME="<<t*dt<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<A_x[j][i];
		p_out<<endl;
	}
	p_out<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<A_y[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"MBASE: ADVECTION FILE OUTPUT SUCCESSFULL AT n = "<<t<<endl;
}
void MBASE::write_bin(ofstream &p_out)
{
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I+1;i++)
			p_out.write((char *) &u[j][i],sizeof(u[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I+1;i++)
			p_out.write((char *) &v[j][i],sizeof(v[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I+1;i++)
			p_out.write((char *) &P[j][i],sizeof(P[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I;i++)
			p_out.write((char *) &u_EW[j][i],sizeof(u_EW[j][i]));
	for(int j=0;j<=J;j++)
		for(int i=0;i<=I+1;i++)
			p_out.write((char *) &v_NS[j][i],sizeof(v_NS[j][i]));
	cout<<"MBASE: INTERMEDIATE FILE OUTPUT SUCCESSFULL"<<endl;
}
void MBASE::read_bin(ifstream &p_in)
{
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I+1;i++)
			p_in.read((char *) &u[j][i],sizeof(u[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I+1;i++)
			p_in.read((char *) &v[j][i],sizeof(v[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I+1;i++)
			p_in.read((char *) &P[j][i],sizeof(P[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I;i++)
			p_in.read((char *) &u_EW[j][i],sizeof(u_EW[j][i]));
	for(int j=0;j<=J;j++)
		for(int i=0;i<=I+1;i++)
			p_in.read((char *) &v_NS[j][i],sizeof(v_NS[j][i]));
	cout<<"MBASE: SOLUTION INITIALIZED SUCCESSFULLY"<<endl;
}
