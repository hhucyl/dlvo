#include "../lbm/Domain.h"


struct myUserData
{
    std::ofstream oss_ss;
    int pinit;
    int GX;
    int GY;
    double g;
    double vmax;
    double nu;
    double R;
    double rhos;
    double rwl;
    double fpl;
    double bbl;
    double sy;
    double ll;
    double tpl;
    double ratiol;
    double ratiot;
    double ratiom;
    double Kn;
    double Gn;
    double Kt;
    double Mu;
    double Eta;
    double Beta;
    double A;
    double kappa;
    double Z;
    double D;
    Vec3_t g0;
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

    //free slip for the top and vmax for the rw top
    // int rwtk = std::ceil(ny-1-dat.fpl);
    // std::cout<<"rwtk "<<rwtk<<std::endl;
    // Vec3_t vb(dat.vmax,0,0);
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for(size_t ix=0; ix<nx; ++ix)
    {
        double *f = dom.F[ix][ny-1][0];
        f[4] = f[2];
        f[7] = f[6];
        f[8] = f[5];
        Vec3_t idx(ix,ny-1,0);
        dom.CalcPropsForCell(idx);
     //    double *f_rwt = dom.F[ix][rwtk][0];
     //    double *f_rwt1 = dom.F[ix][rwtk-1][0];
     //    double rho1 = dom.Rho[ix][rwtk-1][0];
     //    Vec3_t vel1 = dom.Vel[ix][rwtk-1][0];
     //    for(size_t k=0; k<dom.Nneigh; ++k)
    	// {
    	// 	f_rwt[k] = dom.Feq(k,rho1,vb) + f_rwt1[k] - dom.Feq(k,rho1,vel1);
    	// }
    	// idx = ix, rwtk, 0 ;
     //    dom.CalcPropsForCell(idx);

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

    //fix the moving particle if leave the domain in y direction
	int ipp = 0;
	double pvy = 0;
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for(size_t ip=0; ip<dom.Particles.size(); ++ip)
    {
        if(!dom.Particles[ip].IsFree()) continue;
        if(dom.Particles[ip].X(1)<dat.bbl-dat.R)
        {
            dom.Particles[ip].V = 0.0;
            dom.Particles[ip].W = 0.0;
			dom.Particles[ip].X(1) = -100;
			dom.Particles[ip].Xb(1) = -100;
            dom.Particles[ip].FixVeloc();
        }



		if(dom.Particles[ip].X(1)>ny-1-dat.fpl - dat.R - 2)
        {
            #pragma omp atomic
            ipp++;
		}
        #pragma omp atomic
		pvy += dom.Particles[ip].V(1);
	}
	// std::cout<<"ipp "<<ipp<<std::endl;
	// std::cout<<"pinit "<<dat.pinit<<std::endl;
	// dat.pinit = dat.pinit + 1;

	//GY = 1 has bug
	if(ipp<dat.pinit)
    {
        // std::cout<<ipp<<std::endl;

        // double py = dat.sy + (dat.pny-1)*dat.pl;
        // double py = ny-1-3*dat.R;
        Vec3_t pos(0,0,0);
        Vec3_t v(0,0,0);
        Vec3_t w(0,0,0);
        int pnum = dom.Particles.size();
        int GridCount[3] = {0};
        int starty = std::ceil((double)ny-1-(double)dat.fpl);
        int lx = std::ceil((double)nx/(double)dat.GX);
        // int ly = std::ceil((double)dat.fpl/(double)dat.GY);
        // std::cout<<"starty "<<starty<<std::endl;
        // std::cout<<dat.GX<<" "<<dat.GY<<std::endl;

        for(int igx=0; igx< 3; ++igx)
        {
            // GridCount[igx*(dat.GX-1)+igy] = 0;
            int ixs = std::max(0,igx*lx);
            int ixe = std::min((double)nx,(double)(igx+1)*lx);
            int iys = starty;
            int iye = ny;
            // std::cout<<"x "<<ixs<<" "<<ixe<<std::endl;
            // std::cout<<"y "<<iys<<" "<<iye<<std::endl;
            for(int ix=ixs; ix<ixe; ++ix)
            for(int iy=iys; iy<iye; ++iy)
            {
                if(dom.Check[ix][iy][0][0]>0)
                {
                    GridCount[igx] += 1;
                }
                
            }
            // std::cout<<igx<<" "<<igy<<" "<<igx*(dat.GX-1)+igy<<" "<<GridCount[igx*(dat.GX-1)+igy]<<std::endl;
        }

        for(int i=ipp; i<dat.pinit; ++i)
        {
        
            
            int index = std::distance(GridCount,std::min_element(GridCount,GridCount+3));
            // std::cout<<index<<std::endl;
            // for(int i=0; i<3; ++i)
            // {
            // 	std::cout<<GridCount[i]<<" ";
            // }
            // std::cout<<std::endl;

            // int igy = index%(dat.GX-1);
            // int igx = std::floor(index/(dat.GX-1));
            int igx = index;

            double ll = dat.ll;
            // std::cout<<"ll "<<ll<<std::endl;
            int ixs = std::max(ll,(double)igx*lx);
            int ixe = std::min((double)nx-1-ll,(double)(igx+1)*lx);
            int iys = std::ceil(starty+ll);
            int iye = std::floor((double)ny-1-ll);
            // std::cout<<"starty "<<starty<<std::endl;
            // std::cout<<"igx "<<igx<<" x "<<ixs<<" "<<ixe<<std::endl;
            // std::cout<<"igy "<<0<<" y "<<iys<<" "<<iye<<std::endl;
            int N = 0;
            bool overlap = false;
            double px,py;
            while(true)
            {
                overlap = false;
                px = random((double)ixs,(double)ixe);
                py = random((double)iys,(double)iye);
                pos = px, py , 0;
                for(int ip=0; ip<(int)dom.Particles.size(); ++ip)
                {
                    if(pos(0)<0 || pos(1)<0) continue;
        			if(!dom.Particles[ip].IsFree()) continue;

                    if(Norm(pos - dom.Particles[ip].X)<2.5*dat.R || Norm(pos - dom.GhostParticles[ip].X)<2.5*dat.R)
                    {
                        overlap = true;
                        break;
                    }
                }
                N += 1;
                if(!overlap) break;
                if(N>1000) throw new Fatal("Cannot find!!!!");
            }
            
            GridCount[index] += dat.R*4*dat.R;
			
            // Vec3_t dxr(random(-0.5*dat.R,0.5*dat.R),random(-0.5*dat.R,0.5*dat.R),0.0);
            // double px = 0.5*dat.pl+i*dat.pl;
            pos = px, py , 0;
            // pos = pos +dxr;
            // v = Norm(dat.g0)/(2.0*dat.nu)*(-1.0*((pos(1)-dat.sy-dat.H)*(pos(1)-dat.sy-dat.H)-dat.H*dat.H)), 0,0;
            int ipx = std::round(px);
            int ipy = std::round(py);
            ipy = (ipy>(int)ny-1) ? ny-1 : ipy; 
            ipx = (ipx>(int)nx-1) ? nx-1 : ipx; 
            ipx = (ipx<0) ? 0 : ipx;
            v = dom.Vel[ipx][ipy][0](0), 0, 0; 
            dom.Particles.push_back(DEM::Disk(-pnum, pos, v, w, dat.rhos, dat.R, dom.dtdem));
            // std::cout<<pos(0)<<" "<<pos(1)<<std::endl;
            pnum++;
            dom.Particles.back().Ff =  0.0, -M_PI*dat.R*dat.R*(dat.rhos/1.0-1)*dat.g,0.0, 0.0;
            // dom.Particles.back().Ff = 0.0, 0.0, 0.0;
            dom.Particles.back().Kn = dat.Kn;
            dom.Particles.back().Gn = dat.Gn;
            dom.Particles.back().Kt = dat.Kt;
            dom.Particles.back().Mu = dat.Mu;
            dom.Particles.back().Eta = dat.Eta;
            dom.Particles.back().Beta = dat.Beta;
            dom.Particles.back().A = dat.A;
            dom.Particles.back().kappa = dat.kappa;
            dom.Particles.back().Z = dat.Z;
			dom.Particles.back().VdwCutoff = 1.57e-10/dat.ratiol;
            dom.Particles.back().D = dat.D;
        
            
            dom.Particles.back().Rh = 0.80*dat.R;
        }
    
    }


}

void Initial(LBM::Domain &dom, void *UD)
{
     myUserData &dat = (*static_cast<myUserData *> (UD));
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);

    // double py = dat.Ny*2*dat.R + (dat.Ny-1)*dat.gap;
    double sy = dat.sy;
    double H = (double)ny-1-sy;
    double vmax = dat.vmax;
    std::cout<<"sy = "<<sy<<" H = "<<H<<std::endl;
    // std::cout<<"max vel "<<vmax<<std::endl;
    std::cout<<"max vel "<<Norm(dat.g0)/(2.0*dat.nu)*H*H<<std::endl;
    for(size_t ix=0; ix<nx; ++ix)
    for(size_t iy=0; iy<ny; ++iy)
    {
        if(iy<sy)
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
            double Y = yy-sy;
            // double uy = dat.g/(2.0*dat.nu)*(H*(yy-py) - (yy-py)*(yy-py)); 
            double uy = Norm(dat.g0)/(2.0*dat.nu)*(-1.0*((Y-H)*(Y-H)-H*H));
            // double uy = (yy>sy+H) ? vmax : vmax*Y/H;
            
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
    

    char const *infilename = "test_swi_ini.txt";
    String fn;
    fn.Printf("%s",infilename);
    std::fstream ifile(fn.CStr(),std::ios::in);
    if(ifile.fail()) throw new Fatal("FILE WRONG");

    size_t Nproc = 4;
    if(argc>=2) Nproc = atoi(argv[1]); 
    double nu, Re, R, RR, d, rwl, fpl, bbl, pdx, pdy, ratiol, ratiot, ratiom, rho, rhos, Ga, ll,sy,tpl, vmax;
    double Kn, Gn, Kt, Mu, Eta, Beta, A, kappa, Z, D;
    int GX, GY,pinit;
    int Pnx, Pny, fpn;
    size_t nx, ny;
    ifile>>nu;  ifile.ignore(200,'\n');        
    ifile>>Re;  ifile.ignore(200,'\n');
    ifile>>vmax;  ifile.ignore(200,'\n');
    std::cout<<" nu "<<nu<<" Re "<<Re<<" vmax "<<vmax<<std::endl;
    ifile>>R;  ifile.ignore(200,'\n');
    ifile>>RR;  ifile.ignore(200,'\n');
    ifile>>d;  ifile.ignore(200,'\n');
    ifile>>rwl;  ifile.ignore(200,'\n');
    ifile>>fpl;  ifile.ignore(200,'\n');
    ifile>>bbl;  ifile.ignore(200,'\n');
    ifile>>tpl;  ifile.ignore(200,'\n');
    std::cout<<" R "<<R<<" RR "<<RR<<" d "<<d<<" rwl "<<rwl<<" fpl "<<fpl<<" bbl "<<bbl<<" tpl "<<tpl<<std::endl;
    ifile>>Pnx;  ifile.ignore(200,'\n');
    ifile>>Pny;  ifile.ignore(200,'\n');
    ifile>>fpn;  ifile.ignore(200,'\n');
    ifile>>pdx;  ifile.ignore(200,'\n');
    ifile>>pdy;  ifile.ignore(200,'\n');
    std::cout<<" Pnx "<<Pnx<<" Pny "<<Pny<<" fpn "<<fpn<<" pdx "<<pdx<<" pdy "<<pdy<<std::endl;
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
    std::cout<<" GX "<<GX<<" GY "<<GY<<std::endl;

    // double nu = 0.01;
    // int Rn = 10;
    // double R = Rn*1.0;
    // double ppl = 2*R;
    // double ratio = 5;
    // double RR = ratio*R;
    // int Pnx = 3;
    // int Pny = 3; 
    // int pnx = 10;
    // int pny = 8;
    // double pdx = 60;//gap between the large particle
    // double pdy = 60; 
    // double Re = 1e3;

    ifile>>nx;  ifile.ignore(200,'\n');
    ifile>>ny;  ifile.ignore(200,'\n');
    ifile>>ll;  ifile.ignore(200,'\n');
    ifile>>sy;  ifile.ignore(200,'\n');
    ifile>>pinit;  ifile.ignore(200,'\n');
    std::cout<<"nx "<<nx<<" ny "<<ny<<" ll "<<ll<<" sy "<<sy<<" pinit "<<pinit<<std::endl;

    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;


    double gy = Ga*Ga*nu*nu/((8*R*R*R)*(rhos/rho-1));
    std::cout<<"gy = "<<gy<<std::endl;

    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.pinit = fpn;
    my_dat.GX = GX;
    my_dat.GY = GY;
    my_dat.nu = nu;
    my_dat.g = gy;

    double H = (double)ny-1 - sy;
    std::cout<<"H "<<H<<std::endl; 
    Vec3_t g0(2.0*nu*vmax/(H*H),0.0,0.0);
    std::cout<<"gx = "<<g0<<std::endl;

    my_dat.g0 = g0;
    my_dat.R = R;
    my_dat.vmax = vmax;
    my_dat.rwl = rwl;
    my_dat.fpl = fpl;
    my_dat.bbl = bbl;
    my_dat.tpl = tpl;
    my_dat.ll = ll;
    my_dat.sy = sy;
    my_dat.ratiol = ratiol;
    my_dat.ratiot = ratiot;
    my_dat.ratiom = ratiom;
    my_dat.Kn = Kn;
    my_dat.Gn = Gn;
    my_dat.Kt = Kt;
    my_dat.Mu = Mu;
    my_dat.Eta = Eta;
    my_dat.Beta = Beta;
    my_dat.A = A/((ratiol/ratiot)*(ratiol/ratiot)*ratiom);
    my_dat.kappa = kappa*ratiol;
    my_dat.Z = Z*ratiot*ratiot/(ratiol*ratiom);
    my_dat.D = D;
    // std::cout<<ipp<<std::endl;


    dom.Nproc = Nproc;       
    my_dat.rhos = rhos;
    Vec3_t v0(0.0,0.0,0.0);
    dom.dtdem = 0.01*dt;
    
    //initial
    // dom.Initial(rho,v0,g0);
    
    dom.InitialFromH5("test_swi1_0007.h5",g0);
    bool iscontinue = true;
	// Initial(dom,dom.UserData);
    // bool iscontinue = false;
    
    

    if(!iscontinue)
    {
        //fix
        Vec3_t pos(0.0,0.0,0.0);
        Vec3_t v(0.0,0.0,0.0);
        Vec3_t w(0.0,0.0,0.0);
        for(int i=0; i<pinit; ++i)
        {
            double x;
            double y;
            double r;
            int tag;
            ifile>>x>>y>>r>>tag; ifile.ignore(200,'\n');
            // std::cout<<" x "<<x<<" y "<<y<<" r "<<r<<" t "<<tag<<std::endl;
            pos = x,y,0;
            dom.Particles.push_back(DEM::Disk(i, pos, v, w, rhos, r, dom.dtdem));
            if(tag<0){
            	dom.Particles[i].FixVeloc();	
            }else{
            	dom.Particles[i].V(0) = vmax;
            } 
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
        dom.Particles[ip].Rh = dom.Particles[ip].R -d;
        
        // if(dom.Particles[ip].IsFree())
        // {
            
        //     dom.Particles[ip].Rh = dom.Particles[ip].R -d;

        // }else{
        //     dom.Particles[ip].Rh = dom.Particles[ip].R -d;
        // }
        


    }

    //rwp
    double Dm = nu/1e3;
    int sy_k = std::ceil(sy);
    int ey_k = std::floor(sy_k+rwl);
    int RWPN = 20;

    for(int ix=0; ix<nx; ++ix)
    for(int iy=sy_k+1; iy<ey_k+1;++iy)
    {
        for(int ir=0; ir<RWPN; ++ir)
        {
            Vec3_t xr(random(ix-1,ix),random(iy,iy+1),0);
            dom.RWParticles.push_back(RW::Particle(xr,Dm));
        }
    }

    std::cout<<"RW particles complete "<<std::endl;
    std::cout<<"RW particles NUM "<<dom.RWParticles.size()<<"Dm "<<Dm<<std::endl;

    double Tf = 1e6;
    dom.IsF = true;    
    dom.IsRW = true;    
    double dtout = 2e3;
    dom.Box = 0.0, nx-1, 0.0;
    dom.modexy = 0;
    dom.Box1 = bbl, ey_k, 0.0;
    dom.modexy1 = 1;
    //solving
    dom.SolveIBM( Tf, dtout, "test_swi_r1", Setup, NULL);
    
}MECHSYS_CATCH