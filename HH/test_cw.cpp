#include "../lbm/Domain.h"


struct myUserData
{
    std::ofstream oss_ss;
    double vmax;
    Vec3_t g0;
    double H;
    double sy;
    double Kn;
    double Gn;
    double Kt;
    double A;
    double Z;
    double kappa;
    double Mu;
    double Eta;
    double Beta;
    double rhos;
    double g;
    double D;
};

double random(double a, double b)
{
    double rn = std::rand()/((double) RAND_MAX);
    return (b-a)*rn+a;
}
  

void Setup(LBM::Domain &dom, void *UD)
{
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);
    myUserData &dat = (*static_cast<myUserData *> (UD));
    Vec3_t vb(dat.vmax,0,0);
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for(size_t ix=0; ix<nx; ++ix)
    {
        double *f = dom.F[ix][ny-1][0];
        double *f1 = dom.F[ix][ny-2][0];
        double rho1 = dom.Rho[ix][ny-2][0];
        Vec3_t vel1 = dom.Vel[ix][ny-2][0];
        for(size_t k=0; k<dom.Nneigh; ++k)
        {
            f[k] = dom.Feq(k,rho1,vb) + f1[k] - dom.Feq(k,rho1,vel1);
        }
        Vec3_t idx(ix,ny-1,0);
        dom.CalcPropsForCell(idx);
    }

    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for(size_t ix=0; ix<nx; ++ix)
    {
    	double *f = dom.F[ix][0][0];
    	double *f1 = dom.F[ix][1][0];
    	for(size_t k=0; k<dom.Nneigh; ++k)
    	{
    		f[k] = f1[k];
    	}
    	Vec3_t idx(ix,0,0);
        dom.CalcPropsForCell(idx);
    }


}

void Initial(LBM::Domain &dom, void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);

    // double py = dat.Ny*2*dat.R + (dat.Ny-1)*dat.gap;
    double py = dat.sy;
    double H = dat.H;
    std::cout<<"py = "<<py<<" H = "<<H<<std::endl;
    for(size_t ix=0; ix<nx; ++ix)
    for(size_t iy=0; iy<ny; ++iy)
    {
        if(iy<py)
        {
            Vec3_t vtemp(0, 0, 0);
            // Vec3_t vtemp((double) dat.vb, 0, 0);
            dom.Rho[ix][iy][0] = 1.0;
            dom.Vel[ix][iy][0] = vtemp;
            dom.BForce[ix][iy][0] = dat.g0;
            for(size_t k=0; k<dom.Nneigh; ++k)
            {
                dom.F[ix][iy][0][k] = dom.Feq(k,1.0,vtemp);            
                dom.Ftemp[ix][iy][0][k] = dom.Feq(k,1.0,vtemp);            
            }
        }else{
            double yy = (double) iy;
            double Y = yy-py;
            double uy = Y/H*dat.vmax; 
            Vec3_t vtemp(uy, 0, 0);
            // Vec3_t vtemp((double) dat.vb, 0, 0);
            dom.Rho[ix][iy][0] = 1.0;
            dom.Vel[ix][iy][0] = vtemp;
            dom.BForce[ix][iy][0] = dat.g0;
            for(size_t k=0; k<dom.Nneigh; ++k)
            {
                dom.F[ix][iy][0][k] = dom.Feq(k,1.0,vtemp);            
                dom.Ftemp[ix][iy][0][k] = dom.Feq(k,1.0,vtemp);            
            }
        }
        

    }
    

}


int main (int argc, char **argv) try
{
    std::srand((unsigned)time(NULL)); 
    char const *infilename = "test_w_ini.txt";
    String fn;
    fn.Printf("%s",infilename);
    std::fstream ifile(fn.CStr(),std::ios::in);
    if(ifile.fail()) throw new Fatal("FILE WRONG");
   
    
    size_t Nproc = 4;
    if(argc>=2) Nproc = atoi(argv[1]); 
    
    //fluid and particles 
    double nu, Re, R, RR, ppl,pppl, pdx, pdy, ratiol, ratiot, ratiom, rho, rhos, Ga;
    double Kn, Gn, Kt, Mu, Eta, Beta, A, kappa, Z, D, pl, sy, d, mag;
    int GX, GY,pinit, PN;
    int Pnx, Pny, pnx, pny, H;
    size_t nx, ny;
    ifile>>nu;  ifile.ignore(200,'\n');        
    ifile>>Re;  ifile.ignore(200,'\n');
    std::cout<<" nu "<<nu<<" Re "<<Re<<std::endl;
    ifile>>R;  ifile.ignore(200,'\n');
    ifile>>RR;  ifile.ignore(200,'\n');
    ifile>>d;  ifile.ignore(200,'\n');
    ifile>>ppl;  ifile.ignore(200,'\n');
    ifile>>pppl;  ifile.ignore(200,'\n');
    std::cout<<" R "<<R<<" RR "<<RR<<" d "<<d<<" ppl "<<ppl<<" pppl "<<pppl<<std::endl;
    ifile>>Pnx;  ifile.ignore(200,'\n');
    ifile>>Pny;  ifile.ignore(200,'\n');
    ifile>>pnx;  ifile.ignore(200,'\n');
    ifile>>pny;  ifile.ignore(200,'\n');
    ifile>>pdx;  ifile.ignore(200,'\n');
    ifile>>pdy;  ifile.ignore(200,'\n');
    ifile>>mag;  ifile.ignore(200,'\n');
    std::cout<<" Pnx "<<Pnx<<" Pny "<<Pny<<" pnx "<<pnx<<" pny "<<pny<<" pdx "<<pdx<<" pdy "<<pdy<<" mag "<<mag<<std::endl;
    ifile>>nx;  ifile.ignore(200,'\n');
    ifile>>ny;  ifile.ignore(200,'\n');
    ifile>>H;  ifile.ignore(200,'\n');
    ifile>>pl;  ifile.ignore(200,'\n');
    std::cout<<"nx "<<nx<<" ny "<<ny<<" H "<<" pl "<<pl<<std::endl;
    ifile>>ratiol;  ifile.ignore(200,'\n');
    ifile>>ratiot;  ifile.ignore(200,'\n');
    ifile>>ratiom;  ifile.ignore(200,'\n');
    std::cout<<" ratiol "<<ratiol<<" ratiot "<<ratiot<<" ratiom "<<ratiom<<std::endl;
    ifile>>rho;  ifile.ignore(200,'\n');
    ifile>>rhos;  ifile.ignore(200,'\n');
    ifile>>Ga;  ifile.ignore(200,'\n');
    std::cout<<" rho "<<rho<<" rhos "<<rhos<<" Ga "<<Ga<<std::endl;
    ifile>>Kn;  ifile.ignore(200,'\n');
    ifile>>Gn;  ifile.ignore(200,'\n');
    ifile>>Kt;  ifile.ignore(200,'\n');
    ifile>>Mu;  ifile.ignore(200,'\n');
    ifile>>Eta;  ifile.ignore(200,'\n');
    ifile>>Beta;  ifile.ignore(200,'\n');
    ifile>>A;  ifile.ignore(200,'\n');
    ifile>>kappa;  ifile.ignore(200,'\n');
    ifile>>Z;  ifile.ignore(200,'\n');
    ifile>>D;  ifile.ignore(200,'\n');
    std::cout<<" Kn "<<Kn<<" Gn "<<Gn<<" Kt "<<Kt<<" Mu "<<Mu<<" Eta "<<Eta<<" Beta "<<Beta<<" A "<<A<<" kappa "<<kappa<<" Z "<<Z<<" D "<<D<<std::endl;;
    ifile>>GX;  ifile.ignore(200,'\n');
    ifile>>GY;  ifile.ignore(200,'\n');
    ifile>>pinit;  ifile.ignore(200,'\n');
    std::cout<<" GX "<<GX<<" GY "<<GY<<" pinit "<<pinit<<std::endl;
    ifile>>sy;  ifile.ignore(200,'\n');
    ifile>>PN;  ifile.ignore(200,'\n');
    std::cout<<" sy "<<sy<<" PN "<<PN<<std::endl;

    

    
    
    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;
    double vmax = nu*Re/H*1.5;
    std::cout<<"vmax "<<vmax<<std::endl;
    
    std::cout<<"H = "<<H<<std::endl;
    std::cout<<"sy = "<<sy<<std::endl;

    double gy = Ga*Ga*nu*nu/((8*R*R*R)*(rhos/rho-1));
    std::cout<<"gy = "<<gy<<std::endl;

    std::cout<<nx<<" "<<ny<<" "<<nz<<std::endl;
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);
    myUserData my_dat;
    dom.UserData = &my_dat;
    
    my_dat.g0 = 0,0,0;
   
    my_dat.H = H;
    my_dat.vmax = vmax;
    my_dat.sy = sy;
    
    my_dat.Kn = Kn;
    my_dat.Gn = Gn;
    my_dat.Kt = Kt;
    my_dat.Mu = Mu;
    my_dat.Eta = Eta;
    my_dat.Beta = Beta;
    my_dat.D = D;
    my_dat.g = gy;

    dom.Nproc = Nproc;       
    my_dat.rhos = rhos;
    Vec3_t v0(0.0,0.0,0.0);
    dom.dtdem = 1*dt;//if don't have moving particle this value is better set as 1, if not, 0.01 
    
    //initial
    
    // dom.Initial(rho,v0,g0);

    // Initial(dom,dom.UserData);//together
    // bool iscontinue = false;

    dom.InitialFromH5("test_w1_0090.h5",my_dat.g0);//together    
    bool iscontinue = true;

    
    

    if(!iscontinue)
    {
        Vec3_t pos(0.0,0.0,0.0);
        Vec3_t v(0.0,0.0,0.0);
        Vec3_t w(0.0,0.0,0.0);
        

        
        for(int i=0; i<PN; ++i)
        {
            double x;
            double y;
            double r;
            int tag;
            ifile>>x>>y>>r>>tag; ifile.ignore(200,'\n');
            std::cout<<" x "<<x<<" y "<<y<<" r "<<r<<" t "<<tag<<std::endl;
            pos = x,y,0;
            dom.Particles.push_back(DEM::Disk(i, pos, v, w, rhos, r, dom.dtdem));
            if(tag<0) dom.Particles[i].FixVeloc();
        }
        
       
       

    }
    

    std::cout<<"Particles number = "<<dom.Particles.size()<<std::endl;
    for(size_t ip=0; ip<dom.Particles.size(); ++ip)
    {
        
        // dom.Particles[ip].Ff = 0.0, 0.0, 0.0;
        dom.Particles[ip].Kn = my_dat.Kn;
        dom.Particles[ip].Gn = my_dat.Gn;
        dom.Particles[ip].Kt = my_dat.Kt;
        dom.Particles[ip].Mu = my_dat.Mu;
        dom.Particles[ip].Eta = my_dat.Eta;
        dom.Particles[ip].Beta = my_dat.Beta;
        dom.Particles[ip].A = my_dat.A;
        dom.Particles[ip].kappa = my_dat.kappa;
        dom.Particles[ip].Z = my_dat.Z;
        dom.Particles[ip].VdwCutoff = 1.57e-10/ratiol;
        dom.Particles[ip].D = my_dat.D;
        dom.Particles[ip].M = M_PI*dom.Particles[ip].R*dom.Particles[ip].R*rhos;
        dom.Particles[ip].Ff =  0.0, -M_PI*dom.Particles[ip].R*dom.Particles[ip].R*(rhos/rho-1)*my_dat.g,0.0;
        if(dom.Particles[ip].IsFree())
        {
            
            dom.Particles[ip].Rh = dom.Particles[ip].R - d;

        }else{
            dom.Particles[ip].Rh = dom.Particles[ip].R - d;
        }
        


    }

    //fluid and particles fisnish 

    
    //rw Particles start
    double Dm = 1e-2*nu;
    int RWP = 5e3;//adsorption limit
    dom.IsRW = true;
    Vec3_t xt(0,0,0);
    double pi = 3.1415926;
    bool is_rwp_release_from_pa = true;
    if(is_rwp_release_from_pa)
    {	
   
		for(size_t ip=0; ip<dom.Particles.size(); ++ip)
	    {    
	        DEM::Disk *Pa = &dom.Particles[ip];
	        if(Pa->IsFree())//release 
	        {
	            for(int ir=0; ir<RWP; ++ir)
	            {
	                double d = ((double) ir)/((double) RWP);
	                Vec3_t e(Pa->Rh*std::cos(d*2*pi),Pa->Rh*std::sin(d*2*pi),0);
	                xt = Pa->X + e;
	                dom.RWParticles.push_back(RW::Particle(xt,Dm));
	                
	                RW::Particle *RWP = &dom.RWParticles.back();
	                RWP->AD = true;//if this rwp is adsorped it is need to be set before solve
	                RWP->ip = ip;
                    RWP->Leave2 = &RW::Particle::LeaveReflect;
	            }
	            Pa->Alimit = RWP;//adsorption limit 
	            Pa->Alimit0 = RWP;//initial adsorption amount
	            Pa->Pa = 0.5;// adsorption probablitty
	            Pa->Pd = 0.1;// desorption probablitty 
	            // note: Pa and Pd can vary from particle to particle this might be something interesting to simulate
	        }else{
	            Pa->Alimit = RWP*RR/R; 
	            Pa->Alimit0 = 0;
	            Pa->Pa = 0.5;
	            Pa->Pd = 0.1;
	        }
	    }
			


    	

    }else{
    	//to do if needed 
    }
    
    
    std::cout<<"RW particles complete "<<std::endl;
    std::cout<<"RW particles NUM "<<dom.RWParticles.size()<<std::endl;


    double Tf = 1e4;
    dom.IsF = true;    
    double dtout = 1e2;
    dom.Box = 0.0, nx-1, 0.0;
    dom.modexy = 0;
    dom.Box1 = 0.0, ny-1, 0.0;
    dom.modexy1 = 1;
    //solving
    dom.SolveIBM( Tf, dtout, "test_c1", Setup, NULL);
    
}MECHSYS_CATCH
