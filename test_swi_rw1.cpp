#include "./lbm/Domain.h"


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
    
    size_t Nproc = 8;
    if(argc>=2) Nproc = atoi(argv[1]); 
    

    double nu = 0.1;
    double Dm = 1e-2*nu;
    int Rn = 10;
    double R = Rn*1.0;
    double ratio = 1.5;
    double RR = ratio*R;
    int Pnx = 2;
    int Pny = 1; 
    int pnx = 2;
    int pny = 1;
    double pdx = 10;//gap between the large particle
    double pdy = 10; 
    double Re = 1e1;
    
    

    size_t nx = std::ceil(Pnx*(RR*2 + pdx));
    double pl = (double) nx/pnx;
    double ppl = 2*R; //distance above the swi 
    double pppl = 2*R; //distance below the top boundary
    // std::cout<<"pl "<<pl<<std::endl;
    // std::cout<<"pll "<<ppl + (pny-1)*pl<<std::endl;
    double H =  ppl + (pny-1)*pl + pppl + 1;
    double sy = (Pny-1)*(std::sqrt(3)*RR+pdy) + 2*RR;
    size_t ny = std::ceil(sy + H);
    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;
    // double ratiol = 1e-6; //Lr/Ll
    // double ratiot = 1e-8; //Tr/Tl
    // double ratiom = 1e-15; //Mr/Ml
    double vmax = nu*Re/H*1.5;
    std::cout<<"vmax "<<vmax<<std::endl;
    
    std::cout<<"H = "<<H<<std::endl;
    std::cout<<"sy = "<<sy<<std::endl;

    
    
    
    double rho = 1.0;
    double rhos = 2.7;
    double Ga = 5.0;
    double gy = Ga*Ga*nu*nu/((8*R*R*R)*(rhos/rho-1));
    std::cout<<"gy = "<<gy<<std::endl;
    //nu = 1.0/30.0;
    std::cout<<nx<<" "<<ny<<" "<<nz<<std::endl;
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);
    myUserData my_dat;
    dom.UserData = &my_dat;
    
    my_dat.g0 = 0,0,0;
   
    my_dat.H = H;
    my_dat.vmax = vmax;
    my_dat.sy = sy;
    
    my_dat.Kn = 10;
    my_dat.Gn = 0;
    my_dat.Kt = 0;
    my_dat.Mu = 0;
    my_dat.Eta = 0;
    my_dat.Beta = 0;
    my_dat.D = 1;
    my_dat.g = gy;

    dom.Nproc = Nproc;       
    my_dat.rhos = rhos;
    Vec3_t v0(0.0,0.0,0.0);
    dom.dtdem = 0.01*dt;
    
    //initial
    // dom.InitialFromH5("test_swi1_0007.h5",g0);
    // dom.Initial(rho,v0,g0);
    Initial(dom,dom.UserData);
    bool iscontinue = false;
    // bool iscontinue = true;

    
    

    if(!iscontinue)
    {
        //fix
        Vec3_t pos(0.0,0.0,0.0);
        Vec3_t v(0.0,0.0,0.0);
        Vec3_t w(0.0,0.0,0.0);
        double py = 0;
        double px = 0;
        int pnum = 0;
        for(int j=0; j<Pny; ++j)
        {
            py = 0.0 + j*(std::sqrt(3)*RR+pdy) + RR + 1;
            int temp = j%2==0 ? Pnx : Pnx+1;
            for(int i = 0; i<temp; ++i)
            {
                px = j%2==0 ? (0.5*pdx+RR+i*(2*RR+pdx)) : i*(2*RR+pdx);
                pos = px, py, 0;
                dom.Particles.push_back(DEM::Disk(pnum, pos, v, w, rhos, RR, dom.dtdem));
                dom.Particles[pnum].FixVeloc();
                
                pnum++;
            }
        }
        

        //move
        for(int j=0; j<pny; ++j)
        {
            py = sy + j*pl + ppl;
            for(int i=0; i<pnx; ++i)
            {
                Vec3_t dxr(random(-0.2*R,0.2*R),random(-0.2*R,0.2*R),0.0);
                // Vec3_t dxr(0.0,0.0,0.0);
                px = 0.5*pl+i*pl;
                pos = px, py , 0;
                pos = pos+dxr;
                v = (pos(1)-sy)/H*vmax, 0,0;
                dom.Particles.push_back(DEM::Disk(-pnum, pos, v, w, rhos, R, dom.dtdem));
                // std::cout<<pos(0)<<" "<<pos(1)<<std::endl;
                pnum++;
            }
        }

        
        
        
            
            // for(int i=0; i< pnx*pny; ++i)
            // {
            //     double x;
            //     double y;
            //     ifile>>x>>y;
            //     pos = x,y+sy,0;
            //     v = Norm(g0)/(2.0*nu)*(-1.0*((pos(1)-sy-H)*(pos(1)-sy-H)-H*H)), 0,0;
            //     dom.Particles.push_back(DEM::Disk(-pnum, pos, v, w, rhos, R, dom.dtdem));
            //     pnum++;
            //     // std::cout<<x<<" "<<y<<std::endl;

            // }
       

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
        // dom.Particles[ip].Beta = my_dat.Beta;
        // dom.Particles[ip].A = my_dat.A;
        // dom.Particles[ip].kappa = my_dat.kappa;
        // dom.Particles[ip].Z = my_dat.Z;
        // dom.Particles[ip].A = 2e-20/((ratiol/ratiot)*(ratiol/ratiot));
        // dom.Particles[ip].kappa = 1e8*ratiol;
        // dom.Particles[ip].Z = 1e-11*ratiot*ratiot/ratiol;
        // dom.Particles[ip].bbeta = 0.3;
        // dom.Particles[ip].epsilon = 2.05769e-20/((ratiol/ratiot)*(ratiol/ratiot));
        // dom.Particles[ip].s = 200e-9/ratiol;
        // dom.Particles[ip].Lc = 100e-9/ratiol;
        // dom.Particles[ip].l = 3.04e-10/ratiol;
        //dom.Particles[ip].VdwCutoff = std::sqrt(dom.Particles[ip].A/(12.0*dom.Particles[ip].Z*dom.Particles[ip].kappa));
        // dom.Particles[ip].VdwCutoff = 1.57e-10/ratiol;
        dom.Particles[ip].D = my_dat.D;
        // dom.Particles[ip].R = R;
        dom.Particles[ip].M = M_PI*dom.Particles[ip].R*dom.Particles[ip].R*rhos;
        dom.Particles[ip].Ff =  0.0, -M_PI*dom.Particles[ip].R*dom.Particles[ip].R*(rhos/rho-1)*my_dat.g,0.0;
        if(dom.Particles[ip].IsFree())
        {
            
            dom.Particles[ip].Rh = 0.8*R;

        }else{
            dom.Particles[ip].Rh = RR-(R-0.8*R);
        }
        


    }
    
    //rw Particles
    int RWP = 5e3;
    dom.IsRW = true;
    Vec3_t xt(0,0,0);
    double pi = 3.1415926;
    for(size_t ip=0; ip<dom.Particles.size(); ++ip)
    {    
        DEM::Disk *Pa = &dom.Particles[ip];
        if(Pa->IsFree())
        {
            for(int ir=0; ir<RWP; ++ir)
            {
                double d = ((double) ir)/((double) RWP);
                Vec3_t e(Pa->Rh*std::cos(d*2*pi),Pa->Rh*std::sin(d*2*pi),0);
                xt = Pa->X + e;
                dom.RWParticles.push_back(RW::Particle(xt,Dm));
                
                RW::Particle *RWP = &dom.RWParticles.back();
                RWP->AD = true;
                RWP->ip = ip;
            }
            Pa->Alimit = RWP; 
            Pa->Alimit0 = RWP;
            Pa->Pa = 0.5;
            Pa->Pd = 0.1;
        }else{
            Pa->Alimit = RWP*RR/R; 
            Pa->Alimit0 = 0;
            Pa->Pa = 0.5;
            Pa->Pd = 0.1;
        }
    }
    
    std::cout<<"RW particles complete "<<std::endl;
    std::cout<<"RW particles NUM "<<dom.RWParticles.size()<<std::endl;


    double Tf = 1e5;
    dom.IsF = true;    
    double dtout = 1e2;
    dom.Box = 0.0, nx-1, 0.0;
    dom.modexy = 0;
    dom.Box1 = 0.0, ny-1, 0.0;
    dom.modexy1 = 1;
    //solving
    dom.SolveIBM( Tf, dtout, "test_swi_rw1", Setup, NULL);
    
}MECHSYS_CATCH
