#include "./lbm/Domain.h"


struct myUserData
{
    std::ofstream oss_ss;
    double g;
    double nu;
    double R;
    double rhos;
    double ratiol;
    double ratiot;
};
void Setup(LBM::Domain &dom, void *UD)
{
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);
    myUserData &dat = (*static_cast<myUserData *> (UD));

   

}

void Initial(LBM::Domain &dom, void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);
    for(size_t ix=0; ix<nx; ++ix)
    for(size_t iy=0; iy<ny; ++iy)
    {
        // Vec3_t vtemp((double) dat.vb*iy/(ny-1), 0, 0);
        // Vec3_t vtemp((double) dat.vb, 0, 0);
        Vec3_t vtemp(0, -0.1, 0);
        
        dom.Rho[ix][iy][0] = 1.0;
        dom.Vel[ix][iy][0] = vtemp;
        dom.BForce[ix][iy][0] = 0.0, 0.0, 0.0;
        for(size_t k=0; k<dom.Nneigh; ++k)
        {
            dom.F[ix][iy][0][k] = dom.Feq(k,1.0,vtemp);            
            dom.Ftemp[ix][iy][0][k] = dom.Feq(k,1.0,vtemp);            
        }

    }
}

double random(double a, double b)
{
    double rn = std::rand()/((double) RAND_MAX);
    return (b-a)*rn+a;
}
  
int main (int argc, char **argv) try
{
    std::srand((unsigned)time(NULL));    
    
    size_t Nproc = 1;
    double nu = 0.01;
    int Rn = 20;
    double R = Rn*1.0;
    double RR = 10;
    int Pnx = 5;//big particle number
    int pnx = 7;
    int pny = 3;
    double spy = 10.0;
    // double vb = 0.01; 
    double pdx = 0;//gap between big particle
    
    if(argc>=2) Nproc = atoi(argv[1]); 

    size_t nx = std::ceil(Pnx*(R*2+pdx));
    size_t ny = 200;
    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;
    double ratiol = 10e-6/(2.0*R); //Lr/Ll
    double ratiot = 1/1; //Tr/Tl     
    std::cout<<"R = "<<R<<std::endl;
    double rho = 1.0;
    double rhos = 2.7;
    double gy = 1e-5;
    std::cout<<"gy = "<<gy<<std::endl;

    //nu = 1.0/30.0;
    std::cout<<nx<<" "<<ny<<" "<<nz<<std::endl;
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.nu = nu;
    my_dat.g = gy;
    my_dat.R = R;
    my_dat.rhos = rhos;
    Vec3_t g0(0.0,0.0,0.0);
    dom.Nproc = Nproc;       

    //dom.Isq = true;
    // dom.IsF = false;
    // dom.IsFt = false;
    my_dat.rhos = rhos;
    Vec3_t v0(0.0,0.0,0.0);
    dom.dtdem = 0.01*dt;
    //initial
    // dom.InitialFromH5("test_wu_1_0001.h5",g0);
    dom.Initial(rho,v0,g0);
    // Initial(dom,dom.UserData);
    bool iscontinue = false;
    // bool iscontinue = true;
    
    
    
    if(!iscontinue)
    {
        //fix
        Vec3_t pos(0.0,0.0,0.0);
        Vec3_t pos1(0.0,0.0,0.0);
        Vec3_t dxp(0.0,0.0,0.0);
        Vec3_t v(0.0,0.0,0.0);
        Vec3_t w(0.0,0.0,0.0);
        double py = 0;
        double px = 0;
        int pnum = 0;
        for(int i=0; i<Pnx; ++i)
        {
            px = 20 + 2*i*R;
            pos = px, 20.0, 0;
            dom.Particles.push_back(DEM::Disk(pnum, pos, v, w, rhos, R, dom.dtdem));
            dom.Particles[pnum].FixVeloc();
            pnum++;
        }

        for(int i=0; i<Pnx; ++i)
        {
            px = 20 + 2*i*R;
            pos = px, ny-1-20.0, 0;
            dom.Particles.push_back(DEM::Disk(-1, pos, v, w, rhos, R, dom.dtdem));
        }
        //move
        double plx = ((double) nx - pnx*2*RR)/pnx;
        double ply = ((double) ny - 2*2*R - 2*spy - pny*2*RR)/pny;
        std::cout<<plx<<" "<<pnx<<std::endl;
        std::cout<<ply<<" "<<pny<<std::endl;
        for(int j=0; j<pny; ++j) 
        {
            py = spy + 2*R + 2*j*(RR+ply);
            for(int i=0; i<pnx; ++i)
            {
                px = RR + 0.5*plx + 2*i*RR + i*plx; 
                pos = px, py, 0;
                dom.Particles.push_back(DEM::Disk(0, pos, v, w, rhos, RR, dom.dtdem));
            }

        }


    }

    std::cout<<"Particles number = "<<dom.Particles.size()<<std::endl;
    for(size_t ip=0; ip<dom.Particles.size(); ++ip)
    {
        
        dom.Particles[ip].Kn = 1.0;
        dom.Particles[ip].Gn = 1.0;
        dom.Particles[ip].Kt = 0.0;
        dom.Particles[ip].Mu = 0.0;
        dom.Particles[ip].Eta = 0.0;
        dom.Particles[ip].Beta = 0.0;
        dom.Particles[ip].A = 2e-20/((ratiol/ratiot)*(ratiol/ratiot));
        dom.Particles[ip].kappa = 1e9*ratiol;
        dom.Particles[ip].Z = 1e-11*ratiot*ratiot/ratiol;
        // dom.Particles[ip].A = 2e-20/((ratiol/ratiot)*(ratiol/ratiot));
        // dom.Particles[ip].kappa = 1e8*ratiol;
        // dom.Particles[ip].Z = 1e-11*ratiot*ratiot/ratiol;
        // dom.Particles[ip].bbeta = 0.3;
        // dom.Particles[ip].epsilon = 2.05769e-20/((ratiol/ratiot)*(ratiol/ratiot));
        // dom.Particles[ip].s = 200e-9/ratiol;
        // dom.Particles[ip].Lc = 100e-9/ratiol;
        // dom.Particles[ip].l = 3.04e-10/ratiol;
        dom.Particles[ip].VdwCutoff = std::sqrt(dom.Particles[ip].A/(12.0*dom.Particles[ip].Z*dom.Particles[ip].kappa));
        dom.Particles[ip].D = 2;
        if(dom.Particles[ip].Tag == 0)
        {
            dom.Particles[ip].Ff =  0.0, 0.0, 0.0;
        }else{
            dom.Particles[ip].Ff =  0.0, -M_PI*dom.Particles[ip].R*dom.Particles[ip].R*(rhos/rho-1)*my_dat.g,0.0;
        }
        
        dom.Particles[ip].Rh = 0.8*dom.Particles[ip].R;

        
        


    }
    
    

    

    double Tf = 1e4;
    dom.IsF = true;    
    double dtout = 1e2;
    // periodic in x
    dom.Box = 0.0, nx-1, 0.0;
    dom.modexy = 0;
    // periodic in y
    // dom.Box = 0.0, ny-1, 0.0;
    // dom.modexy = 1;
    
    //solving
    dom.SolveIBM1( Tf, dtout, "test_wu_3", NULL, -1, 0, NULL);
    
}MECHSYS_CATCH