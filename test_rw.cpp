#include "./lbm/Domain.h"
#include<time.h>

double random(double a, double b)
{
    double rn = std::rand()/((double) RAND_MAX);
    return (b-a)*rn+a;
}

int main (int argc, char **argv) try
{
    std::srand((unsigned)time(NULL)); 
    size_t nx = 300;
    size_t ny = 100;
    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;
    double nu = 0.1;
    double Dm = 1e-3;
    std::cout<<nx<<" "<<ny<<" "<<nz<<std::endl;
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);
    dom.Nproc = 4;
    //dem
    double rhos = 2.7;
    double RR = 5;
    Vec3_t pos(0.1*nx,0.5*ny,0.0);

    Vec3_t v(0.0,0.0,0.0);
    Vec3_t w(0.0,0.0,0.0);
    dom.dtdem = 0.01*dt;
    dom.Particles.push_back(DEM::Disk(-1, pos, v, w, rhos, RR, dom.dtdem));
    // dom.Particles[0].FixVeloc();
    pos = 53,50.1,0;
    v = 0,0 ,0;
    dom.Particles.push_back(DEM::Disk(-1, pos, v, w, rhos, RR, dom.dtdem));

    std::cout<<"Particles number = "<<dom.Particles.size()<<std::endl;
    for(size_t ip=0; ip<dom.Particles.size(); ++ip)
    {
        
        // dom.Particles[ip].Ff = 0.0, 0.0, 0.0;
        dom.Particles[ip].Kn = 10;
        dom.Particles[ip].Gn = 0;
        dom.Particles[ip].Kt = 0;
        dom.Particles[ip].Mu = 0;
        dom.Particles[ip].Eta = 0;
        dom.Particles[ip].Beta = 0;
        dom.Particles[ip].A = 0.0;
        dom.Particles[ip].kappa = 0.0;
        dom.Particles[ip].Z = 0.0;
        // dom.Particles[ip].bbeta = 0.3;
        // dom.Particles[ip].epsilon = 2.05769e-20/((ratiol/ratiot)*(ratiol/ratiot));
        // dom.Particles[ip].s = 200e-9/ratiol;
        // dom.Particles[ip].Lc = 100e-9/ratiol;
        // dom.Particles[ip].l = 3.04e-10/ratiol;
        //dom.Particles[ip].VdwCutoff = std::sqrt(dom.Particles[ip].A/(12.0*dom.Particles[ip].Z*dom.Particles[ip].kappa));
        dom.Particles[ip].VdwCutoff = 1.57e-10;
        dom.Particles[ip].D = 1;
        dom.Particles[ip].R = RR;
        dom.Particles[ip].M = M_PI*dom.Particles[ip].R*dom.Particles[ip].R*rhos;
        dom.Particles[ip].Ff =  M_PI*dom.Particles[ip].R*dom.Particles[ip].R*(rhos-1)*1e-4,0.0,0.0;
        std::cout<<dom.Particles[ip].Ff<<std::endl;       
        dom.Particles[ip].Rh = RR;
        dom.Particles[ip].Pa = 0.5;
        dom.Particles[ip].Pd = 0.1;

    }
    
    //rw Particles
    int xs = 36;
    int xe = 46;
    int ys = 20;
    int ye = 80;
    int RWP = 20;
    dom.IsRW = true;
    for(int ix=xs; ix<xe; ++ix)
    for(int iy=ys; iy<ye; ++iy)
    {
        // int x1 = startx[i];
        int x1 = ix;
        int y1 = iy;
        int x2 = x1+1;
        int y2 = y1+1;
        for(int ip=0; ip<RWP; ++ip)
        {
            Vec3_t xt(random(x1,x2),random(y1,y2),0);
            dom.RWParticles.push_back(RW::Particle(xt,Dm));
            

        } 
    }
    
    std::cout<<"RW particles complete "<<std::endl;
    std::cout<<"RW particles NUM "<<dom.RWParticles.size()<<std::endl;

    //initial
    Vec3_t v0(0,0,0);
    Vec3_t g0(0,0,0);
    dom.Initial(1.0,v0,g0);

    //solve
    double Tf = 1e3;
    dom.IsF = true;    
    double dtout = 1;
    dom.Box = 0.0, nx-1, 0.0;
    dom.modexy = 0;

    dom.Box1 = 0.0, ny-1, 0.0;
    dom.modexy1 = 1;
    
    //solving
    dom.SolveIBM( Tf, dtout, "test_rw1", NULL, NULL);

}MECHSYS_CATCH