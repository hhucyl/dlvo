#include "../lbm/Domain.h"

double random(double a, double b)
{
    double rn = std::rand()/((double) RAND_MAX);
    return (b-a)*rn+a;
}

int main (int argc, char **argv) try
{
    std::srand((unsigned)time(NULL)); 

    
    size_t Nproc = 4;
    if(argc>=2) Nproc = atoi(argv[1]); 
    size_t nx = 400;
    size_t ny = 200;
    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;
    double nu = 0.01;
    double R = 10;
    double Rh = 8;
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);

    double rhos = 2.4;
    double rho = 1.0;

    Vec3_t v0(0,0,0);
    Vec3_t g0(1e-7,0,0);
    double pg = 1e-5;
    dom.Initial(rho,v0,g0);
    dom.dtdem = 1*dt;
    dom.Nproc = Nproc;

    for(size_t i=0; i<nx; ++i)
    {
        dom.IsSolid[i][0][0] = true;
        dom.IsSolid[i][ny-1][0] = true;
    }

    bool iscontinue = false;

    if(!iscontinue)
    {
        Vec3_t pos(0.0,0.0,0.0);
        Vec3_t v(0.0,0.0,0.0);
        Vec3_t w(0.0,0.0,0.0);

        //apply wall
        for(size_t i=0; i<20; ++i)
        {
            pos = 2*R*i+R, R, 0;
            // std::cout<<pos<<std::endl;
            dom.Particles.push_back(DEM::Disk(-1, pos, v, w, rhos, R, dom.dtdem));
            dom.Particles.back().FixVeloc();
            pos = 2*R*i+R, (double)ny-1 - R, 0;
            dom.Particles.push_back(DEM::Disk(-1, pos, v, w, rhos, R, dom.dtdem));
            dom.Particles.back().FixVeloc();
        }

        //apply free particle
        std::vector<double> y_pos{40, 70, 100, 130, 160};
        std::vector<double> x_pos{150};
        for(size_t ipx=0; ipx<x_pos.size(); ++ipx)
        for(size_t ipy=0; ipy<y_pos.size(); ++ipy)
        {
            pos = x_pos[ipx], y_pos[ipy], 0;
            dom.Particles.push_back(DEM::Disk(1, pos, v, w, rhos, R, dom.dtdem));
        }

        
    }

    //particles poperties
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
        dom.Particles[ip].A = 0;
        dom.Particles[ip].kappa = 0;
        dom.Particles[ip].Z = 0;
        dom.Particles[ip].VdwCutoff = 1.57e-10;
        dom.Particles[ip].D = 2;
        dom.Particles[ip].M = M_PI*dom.Particles[ip].R*dom.Particles[ip].R*rhos;
        dom.Particles[ip].Ff =  M_PI*dom.Particles[ip].R*dom.Particles[ip].R*(rhos/rho-1)*pg,0.0,0;

        if(dom.Particles[ip].IsFree())
        {
            
            dom.Particles[ip].Rh = Rh;

        }else{
            dom.Particles[ip].Rh = Rh;
        }
    }

    //fluid and particles fisnish
    //rw Particles start
    double Dm = 1e-2*nu;
    int RWP = 20;//adsorption limit

    dom.IsRW = true;

    //Particle adsorption desorption
    for(size_t ip=0; ip<dom.Particles.size(); ++ip)
    {
        DEM::Disk *Pa = &dom.Particles[ip];
        if(Pa->IsFree())
        {     
            Pa->Alimit = 0;//adsorption limit 
            Pa->Alimit0 = 0;//initial adsorption amount
            Pa->Pa = -1;// adsorption probablitty
            Pa->Pd = -1;// desorption probablitty 
            // note: Pa and Pd can vary from particle to particle this might be something interesting to simulate
        }else{
            Pa->Alimit = 0; 
            Pa->Alimit0 = 0;
            Pa->Pa = -1;
            Pa->Pd = -1;
        }
    }

    //rw particles
    size_t xl = 170;
    size_t xu = 180;
    size_t yl = 30;
    size_t yu = 160;
    Vec3_t xt(0,0,0);
    for(size_t ix=xl; ix<xu; ++ix)
    for(size_t iy=yl; iy<yu; ++iy)
    {
        int x1 = (int) ix;
        int y1 = (int) iy;
        int x2 = (int) ix+1;
        int y2 = (int) iy+1;
        for(int ip=0; ip<RWP; ++ip)
        {
            xt = random(x1,x2),random(y1, y2), 0;
            dom.RWParticles.push_back(RW::Particle(xt,Dm));
            dom.RWParticles[ip].Leave2 = &RW::Particle::LeaveReflect;
        }
    }
    
    std::cout<<"RW particles complete "<<std::endl;
    std::cout<<"RW particles NUM "<<dom.RWParticles.size()<<std::endl;
 
    double Tf = 1e5;
    dom.IsF = true;    
    double dtout = 1e2;
    dom.Box = 0.0, nx-1, 0.0;
    dom.modexy = 0;
    dom.Box1 = 2*Rh-1, ny-1-2*Rh, 0.0;
    dom.modexy1 = 1;
    //solving
    dom.SolveIBM( Tf, dtout, "test_pip", NULL, NULL);


}MECHSYS_CATCH