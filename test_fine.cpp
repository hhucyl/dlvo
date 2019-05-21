#include "./lbm/Domain.h"


struct myUserData
{
    std::ofstream oss_ss;
    double g;
    double nu;
    double R;
    double rhos;
    double vb;
};
void Setup(LBM::Domain &dom, void *UD)
{
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);
    myUserData &dat = (*static_cast<myUserData *> (UD));
    Vec3_t vvb(dat.vb,0.0,0.0);
    // std::cout<<vvb<<std::endl;
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for(size_t ix=0;ix<nx;++ix)
    {
        double *f = dom.F[ix][ny-1][0];
        double *f1 = dom.F[ix][ny-2][0];
        double rho1 = dom.Rho[ix][ny-2][0];
        Vec3_t vel1 = dom.Vel[ix][ny-2][0];
        
        for(size_t k=0; k<dom.Nneigh; ++k)
        {
            f[k] = dom.Feq(k,rho1,vvb) + f1[k] - dom.Feq(k,rho1,vel1);
        }
        Vec3_t idx(ix,ny-1,0);
        dom.CalcPropsForCell(idx);
    }

}

void Initial(LBM::Domain &dom, void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);
    for(size_t ix=0; ix<nx; ++ix)
    for(size_t iy=0; iy<ny; ++iy)
    {
        Vec3_t vtemp((double) dat.vb*iy/(ny-1), 0, 0);
        // Vec3_t vtemp((double) dat.vb, 0, 0);
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
    
    size_t Nproc = 12;
    double nu = 0.01;
    int Rn = 5;
    double R = Rn*1.0;
    double len = 3*R;//60.0/5.0*R;
    double ratio = 10.0;
    double RR = ratio*R;
    int Pnx = 6;
    int Pny = 11; 
    int pnx = 10;
    int pny = 10; 
    double vb = 0.01; 
    if(argc>=2) Nproc = atoi(argv[1]); 

    size_t nx = std::ceil(Pnx*RR*2);
    size_t ny = std::ceil((std::sqrt(3)*Pny+1)*RR)+std::ceil(len*pny)+6*Rn;
    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;
    double ratiol = 10e-6/(2.0*R); //Lr/Ll
    double ratiot = 1/1; //Tr/Tl     
    std::cout<<"R = "<<R<<std::endl;
    std::cout<<"RR = "<<RR<<std::endl;
    double rho = 1.0;
    double rhos = 2.7;
    double Ga = 2.0;
    double gy = Ga*Ga*nu*nu/((8*R*R*R)*(rhos/rho-1));
    std::cout<<"gy = "<<gy<<std::endl;
    //nu = 1.0/30.0;
    std::cout<<nx<<" "<<ny<<" "<<nz<<std::endl;
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.nu = nu;
    my_dat.g = gy;
    my_dat.R = R;
    my_dat.vb = vb;
    Vec3_t g0(0.0,0.0,0.0);
    dom.Nproc = Nproc;       

    //dom.Isq = true;
    // dom.IsF = false;
    // dom.IsFt = false;
   
    //initial
    
    my_dat.rhos = rhos;
    Vec3_t v0(0.0,0.0,0.0);
    dom.dtdem = 0.01*dt;
    //fix
    Vec3_t pos(0.0,0.0,0.0);
    Vec3_t pos1(0.0,0.0,0.0);
    Vec3_t dxp(0.0,0.0,0.0);
    Vec3_t v(0.0,0.0,0.0);
    Vec3_t w(0.0,0.0,0.0);
    int pnum = 0;
    for(int j=0; j<Pny; ++j)
    {
        int kx = 0;
        if(j%2!=0) kx = -1;
        for(int i = kx; i<Pnx; ++i)
        {
            dxp = (2.0*i+j%2+1)*RR,(std::sqrt(3)*j+1)*RR,0.0;
            pos1 = pos+dxp;
            dom.Particles.push_back(DEM::Disk(pnum, pos1, v, w, rhos, 0.8*RR, dom.dtdem));
            dom.Particles[pnum].FixVeloc();
            
            pnum++;
        }
    }
    for(size_t ix=0; ix<nx; ++ix)
    {
        dom.IsSolid[ix][0][0] = true;
    }

    //move
    double llen = ((int)nx - pnx*2*R)/pnx;
    dxp = llen,0.0,0.0;
    double ys = (std::sqrt(3)*Pny+1)*RR;
    for(int iy=0; iy<pny; ++iy)
    {
        pos = 0.5*llen,(iy+1)*len+ys,0.0;
        for(int ix=0; ix<pnx; ++ix)
        {
            Vec3_t dxr(random(-0.3*R,0.3*R),random(-0.3*R,0.3*R),0.0);
            // Vec3_t dxr(0.0,0.0,0.0);
            dom.Particles.push_back(DEM::Disk(-pnum, pos+dxr, v, w, rhos, R, dom.dtdem));
            // std::cout<<pos(0)<<" "<<pos(1)<<std::endl;
            pos = pos+dxp;
            pnum++;
        }
    }
    

    std::cout<<"Particles number = "<<dom.Particles.size()<<std::endl;
    for(size_t ip=0; ip<dom.Particles.size(); ++ip)
    {
        dom.Particles[ip].Ff =  0.0, -M_PI*R*R*(rhos/rho-1)*my_dat.g,0.0, 0.0;
        // dom.Particles[ip].Ff = 0.0, 0.0, 0.0;
        dom.Particles[ip].Kn = 1.0;
        dom.Particles[ip].Gn = 1.0;
        dom.Particles[ip].Kt = 0.0;
        dom.Particles[ip].Mu = 0.0;
        dom.Particles[ip].Eta = 0.0;
        dom.Particles[ip].Beta = 0.0;
        dom.Particles[ip].A = 2e-20/((ratiol/ratiot)*(ratiol/ratiot));
        dom.Particles[ip].kappa = 1e9*ratiol;
        dom.Particles[ip].Z = 1e-11*ratiot*ratiot/ratiol;
        dom.Particles[ip].bbeta = 0.3;
        dom.Particles[ip].epsilon = 2.05769e-20/((ratiol/ratiot)*(ratiol/ratiot));
        dom.Particles[ip].s = 200e-9/ratiol;
        dom.Particles[ip].Lc = 100e-9/ratiol;
        dom.Particles[ip].l = 3.04e-10/ratiol;
        dom.Particles[ip].VdwCutoff = std::sqrt(dom.Particles[ip].A/(12.0*dom.Particles[ip].Z*dom.Particles[ip].kappa));
        dom.Particles[ip].D = 2;
        if(dom.Particles[ip].IsFree())
        {
            
            dom.Particles[ip].Rh = 0.8*R;

        }else{
            dom.Particles[ip].Rh = 0.8*RR;

        }
        


    }
    

    // dom.InitialFromH5("test_cong_0999.h5",g0);
    // dom.Initial(rho,v0,g0);
    Initial(dom,dom.UserData);

    double Tf = 1e6;
    dom.IsF = true;    
    double dtout = 1e3;
    dom.Box = 0.0, nx-1, 0.0;
    dom.modexy = 0;
    //solving
    dom.SolveIBM( Tf, dtout, "test_fine", Setup, NULL);
    
}MECHSYS_CATCH