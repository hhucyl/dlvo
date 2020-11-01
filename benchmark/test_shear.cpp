#include "../lbm/Domain.h"
#include <cstring>
struct myUserData
{
    size_t nx;
    size_t ny;
};

double random(double a, double b)
{
    double rn = std::rand()/((double) RAND_MAX);
    return (b-a)*rn+a;
}

void Initial(LBM::Domain &dom, void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    for(size_t ix=0; ix<dat.nx; ++ix)
    for(size_t iy=0; iy<dat.ny; ++iy)
    {
        Vec3_t vtemp(0.01/(dat.ny-1)*iy, 0, 0);
        dom.Rho[ix][iy][0] = 1.0;
        dom.Vel[ix][iy][0] = vtemp;
        dom.BForce[ix][iy][0] = 0.,0.0,0.0;
        for(size_t k=0; k<dom.Nneigh; ++k)
        {
            dom.F[ix][iy][0][k] = dom.Feq(k,1.0,vtemp);            
            dom.Ftemp[ix][iy][0][k] = dom.Feq(k,1.0,vtemp);            
        }
    }
}

int main (int argc, char **argv) try
{
    std::srand((unsigned)time(NULL));
    size_t Nproc = 1;
    if(argc>=2) Nproc = atoi(argv[1]);  
    size_t nx = 100;
    size_t ny = 80;
    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;
    double nu = 0.01;
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);
    dom.Nproc = Nproc;
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.nx = nx;
    my_dat.ny = ny;
	
    Initial(dom,dom.UserData);


    int RWP = 80;
    double Dm = 1e-2;
    

    double Time = 0;
    double Tf = 1e5;
    double tout = 0;
    double dtout = 1e2;
    int idx_out = 0;

    dom.Box = 0.0, nx-1, 0.0;
    dom.modexy = 0;
    dom.Box1 = 0.0, ny-1, 0.0;
    dom.modexy1 = 1;

    //trace particle
    Vec3_t xt(0,0,0);
    for(size_t ix=0; ix<nx; ++ix)
    {
        double x1 = (double) ix-0.5;
        double x2 = (double) ix+0.5;
        for(int ip=0; ip<RWP; ++ip)
        {

            xt = random(x1,x2),random((double)ny-1-1, (double)ny-1), 0;
            dom.RWParticles.push_back(RW::Particle(xt,Dm));
            dom.RWParticles.back().Leave2 = &RW::Particle::LeaveReflect;

        }
    }

    int con[nx] = {0};
    std::vector<int> con_num[nx];
    std::vector<int> recycle;
    while(Time<Tf)
    {
        if (Time>=tout)
        {
            
            String fn;
            fn.Printf("%s_%04d", "test_s", idx_out);
            
            dom.WriteXDMF(fn.CStr());
            idx_out++;
            tout += dtout;
        }
        
        dom.SetZero();
        //collide and streaming
        dom.CollideMRTIBM();
        dom.Stream();
        dom.BounceBack(false);
        dom.CalcProps();

        std::memset(con,0,sizeof(con));

        // #ifdef USE_OMP
        // #pragma omp parallel for schedule(static) num_threads(Nproc)
        // #endif
        for(int i=0;i<(int) dom.RWParticles.size();++i)
        {
            RW::Particle *RWP = &dom.RWParticles[i];
            int x1 = std::round(RWP->X(0));
            bool flag2 = RWP->X(1)>=((double)ny-1-1) && RWP->X(1)<=((double)ny-1);
            if(x1<0) x1=nx-1;
            if(flag2)
            {
                con_num[x1].push_back(i);
                // #pragma omp atomic
                con[x1] += 1;
            }
            
        }

        // for(int i=0; i<nx; ++i)
        // {
        //     // std::cout<<"i "<<i<<" "<<std::endl;
        //     std::cout<<con[i]<<" ";

        // }
        // std::cout<<std::endl;
        // std::cout<<" finish count "<<std::endl;

        //add and erase
        for(int ix=0; ix<nx; ++ix)
        {
            if(con[ix]<RWP)
            {
                for(int ip=con[ix]; ip<RWP; ++ip)
                {
                    xt = random(ix-0.5,ix+0.5),random((double)ny-1-1, (double)ny-1), 0;
                    // if(recycle.size()>0.1)
                    // {
                    //     int index = recycle.back();
                    //     dom.RWParticles[index].X =xt;
                    //     dom.RWParticles[index].Xb =xt;
                    //     recycle.pop_back();
                    //     // std::cout<<recycle.size()<<std::endl;
                    // }else{
                        
                        dom.RWParticles.push_back(RW::Particle(xt,Dm));
                        dom.RWParticles.back().Leave2 = &RW::Particle::LeaveReflect;
                    // }
                    
                }
            }else{
                for(int ip=RWP; ip<con[ix]; ++ip)
                {
                    int c_size = con_num[ix].size(); 
                    int con_num_index = std::round(random(0,(double)c_size-1));
                    int index = con_num[ix][con_num_index];
                    // recycle.push_back(index);
                    dom.RWParticles[index].X(1) = -100;
                }
            }
            con_num[ix].clear();
        }

        // std::memset(con,0,sizeof(con));
        // #ifdef USE_OMP
        // #pragma omp parallel for schedule(static) num_threads(Nproc)
        // #endif
        // for(int i=0;i<(int) dom.RWParticles.size();++i)
        // {
        //     RW::Particle *RWP = &dom.RWParticles[i];
        //     int x1 = std::round(RWP->X(0));
        //     bool flag2 = RWP->X(1)>=((double)ny-1-1) && RWP->X(1)<=((double)ny-1);
        //     if(x1<0) x1=nx-1;
        //     if(flag2)
        //     {
        //         #pragma omp atomic
        //             con[x1] += 1;
        //     }
        // }
        // for(int i=0; i<nx; ++i)
        // {
        //     std::cout<<con[i]<<" ";
        // }
        // std::cout<<std::endl;
        // std::cout<<Time<<std::endl;


        #ifdef USE_OMP
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        #endif
        for(int i=0;i< (int) dom.RWParticles.size();++i)
        {
            
            RW::Particle *RWP = &dom.RWParticles[i];
            std::vector<int> idc{-1,-1};
            
            

            if(RWP->X(1)<-10) continue;
            
            Vec3_t v0(0,0,0); 
            std::vector<Vec3_t> VV{v0,v0,v0,v0};
            std::vector<int> idx{-1,-1,-1,-1};
            // std::cout<<111<<std::endl;

            dom.Pa2GridV(RWP,idx,VV);

            // std::cout<<222<<" "<<idx[0]<<" "<<idx[1]<<" "<<idx[2]<<" "<<idx[3]<<std::endl;
            RWP->Move(VV,idx,dt);
            // std::cout<<333<<std::endl;

            // RWP->Leave(modexy,Box);
            // RWP->LeaveReflect(modexy1,Box1);
            RWP->Boundary(dom.modexy, dom.Box, dom.modexy1, dom.Box1);
            // std::cout<<444<<std::endl;
            
            dom.Pa2Grid(RWP,idc);
            // std::cout<<RWP->X(0)<<" "<<RWP->X(1)<<" "<<idc[0]<<" "<<idc[1]<<std::endl;
            // std::cout<<555<<std::endl;

            int ix = idc[0];
            int iy = idc[1];

            #pragma omp atomic
                dom.Con[ix][iy][0] += 1;
            // std::cout<<666<<std::endl;

        }
        Time += 1;
        dom.Time = Time;
    }

}MECHSYS_CATCH