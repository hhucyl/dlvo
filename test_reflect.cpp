#include "./lbm/Domain.h"
#include <time.h>

double random(double a, double b)
{
    double rn = std::rand()/((double) RAND_MAX);
    return (b-a)*rn+a;
}
  
int main (int argc, char **argv) try
{
    std::srand((unsigned)time(NULL));    
    size_t Nproc = 12;
    size_t Rn = 10;
    double Re = 1e4;
    double nu = 5e-4;
    double vmax = nu*Re/(2*Rn)*1.5;
    std::cout<<"vmax "<<vmax<<std::endl;

    double Dm0 = 2.3e-9;
    double Dm = nu/1e-6*Dm0;
    std::cout<<"Dm0 = "<<Dm0<<" Dm = "<<Dm<<std::endl;
    //rwP
    int RWP = 20;
    if(argc>=2) Nproc = atoi(argv[0]);
    
    double R = (double) Rn;
    size_t nx = std::ceil(4*R);
    size_t ny = std::ceil(4*R);
    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;
    double rho = 1.0;
    double rhos = 2.0;
    std::cout<<"R = "<<R<<std::endl;
    std::cout<<nx<<" "<<ny<<" "<<nz<<std::endl;
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);
    
    
    Vec3_t g0(6e-8,0.0,0.0);
  
    std::cout<<"gx = "<<6e-8<<std::endl;
    dom.Nproc = Nproc;       

    //initial
    
    
    
    Vec3_t pos(0.5*nx,0.5*ny+1,0.0);
    Vec3_t v(0.0,0.0,0.0);
    Vec3_t w(0.0,0.0,0.0);
    //DEM
    dom.dtdem = 1*dt;
    int pnum = 0;
    //fixed
    dom.Particles.push_back(DEM::Disk(pnum, pos, v, w, rhos, R, dom.dtdem));
    dom.Particles[0].FixVeloc();
     
    
    std::cout<<"Particles number = "<<dom.Particles.size()<<std::endl;
    
    for(int ip=0; ip<(int) dom.Particles.size(); ++ip)
    {
        dom.Particles[ip].Ff = 0.0, 0.0, 0.0;
        dom.Particles[ip].Kn = 5;
        dom.Particles[ip].Gn = -0.3;
        dom.Particles[ip].Kt = 2.5;
        dom.Particles[ip].Mu = 0.4;
        dom.Particles[ip].Eta = 0.0;
        dom.Particles[ip].Beta = 0.0;
        dom.Particles[ip].Rh = 0.8*R;
        dom.Particles[ip].FixVeloc();

    }
    
    

    Vec3_t v0(0.0,0.0,0.0);
    dom.IsF = true;
    
    dom.Initial(rho,v0,g0);
    // Initial(dom, dom.UserData);
    // dom.InitialFromH5("test_pbed1_0999.h5",g0);

    //RWParticles
    if(!dom.IsRWContinue)
    {
        
        // std::vector<int> startx{22,65,109,154,199};
        // for(int i=0; i<startx.size(); ++i)
        for(int ix=5; ix<nx-6; ++ix)
        for(int iy=2; iy<ny-3; ++iy)
        {
            // int x1 = startx[i];
            int x1 = ix;
            int y1 = iy;
            int x2 = x1+1;
            int y2 = y1+1;
            for(int ip=0; ip<RWP; ++ip)
            {
                Vec3_t xt(random(x1,x2),random(y1,y2),0);
                int xx = std::round(xt(0));
                int yy = std::round(xt(1));
                Vec3_t xxx(xx,yy,0);
                bool flag = Norm(xxx-dom.Particles[0].X)>(dom.Particles[0].Rh+1);
                if(flag)
                    dom.RWParticles.push_back(RW::Particle(xt,2,Dm));
            } 
        }
        std::cout<<"RW particles complete "<<std::endl;
        std::cout<<"RW particles NUM "<<dom.RWParticles.size()<<std::endl;
    }else{
        for(size_t ip=0; ip<dom.RWParticles.size(); ++ip)
        {
            dom.RWParticles[ip].Dm = Dm;
        }
        std::cout<<"Assign Dm"<<std::endl;
        
    }
    


    double Tf = 1e6;
    double dtout = 1e3;
    dom.Box = 0.0,(double) nx-1, 0.0;
    dom.modexy = 0;
    //solving
    dom.SolvePRW( Tf, dtout, "test_reflect1", NULL, NULL);
    
    return 0;
}MECHSYS_CATCH