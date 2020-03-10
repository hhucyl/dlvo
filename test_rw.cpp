#include "./lbm/Domain.h"
#include <time.h>




int main (int argc, char **argv) try
{
    std::srand((unsigned)time(NULL));    
    size_t Nproc = 4;
    double nu = 5e-4;

    //rwP
    int RWP = 5e4;
    if(argc>=2) Nproc = atoi(argv[0]);
    
    size_t nx = 100;
    size_t ny = 100;
    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;
    double Dm = 0.01;
    //nu = 1.0/30.0;
    std::cout<<nx<<" "<<ny<<" "<<nz<<std::endl;
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);
    
    dom.Nproc = Nproc;       

    //initial
    Vec3_t v0(0.01,0.0,0.0);
    Vec3_t g0(0.0,0.0,0.0);
    dom.IsF = true;
    dom.dtdem = 1*dt;
    dom.Initial(1.0,v0,g0);
   
    //RWParticles

    int nn = 0;
    
    for(int ip=0; ip<RWP; ++ip)
    {
        Vec3_t xt(49,49,0);
        int xx = std::round(xt(0));
        int yy = std::round(xt(1));
        if(dom.Gamma[xx][yy][0]<1e-12)
        {
            dom.RWParticles.push_back(RW::Particle(xt,Dm));
        }
    } 
    dom.IsRW = true;
    
    std::cout<<"RW particles complete "<<std::endl;
    std::cout<<"RW particles NUM "<<dom.RWParticles.size()<<std::endl;

    


    double Tf = 1e4;
    double dtout = 100;
    dom.Box = 0.0,(double) nx-1, 0.0;
    dom.modexy = 0;
    dom.Box1 = 0.0, ny-1, 0.0;
    dom.modexy1 = 1;
    //solving
    dom.SolveIBM( Tf, dtout, "test_rw2", NULL, NULL);
    
    return 0;
}MECHSYS_CATCH
