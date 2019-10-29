#include "./lbm/Domain.h"


struct myUserData
{
    std::ofstream oss_ss;
    int pinit;
    int GX;
    int GY;
    double g;
    double nu;
    double R;
    double rhos;
    double H;
    Vec3_t g0;
    double pl;
    double sy;
    double pny;
    double pnx;
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

    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for(size_t ix=0; ix<nx; ++ix)
    {
        double *f = dom.F[ix][ny-1][0];
        f[4] = f[2];
        f[7] = f[6];
        f[8] = f[5];
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

    //fix the moving particle if leave the domain in y direction
	int ipp = 0;
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for(size_t ip=0; ip<dom.Particles.size(); ++ip)
    {
        if(!dom.Particles[ip].IsFree()) continue;
        if(dom.Particles[ip].X(1)<2*(double) dat.R)
        {
            dom.Particles[ip].V = 0.0;
            dom.Particles[ip].W = 0.0;
			dom.Particles[ip].X(1) = -100;
			dom.Particles[ip].Xb(1) = -100;
            dom.Particles[ip].FixVeloc();

        }
		if(dom.Particles[ip].X(1)>ny-1-dat.H*0.5)
        {
            #pragma omp atomic
            ipp++;
		}

	}

	if(ipp<dat.pinit)
    {
        // std::cout<<ipp<<std::endl;

        // double py = dat.sy + (dat.pny-1)*dat.pl;
        // double py = ny-1-3*dat.R;
        Vec3_t pos(0,0,0);
        Vec3_t v(0,0,0);
        Vec3_t w(0,0,0);
        int pnum = dom.Particles.size();
        int GridCount[dat.GX*dat.GY];
        int starty = std::ceil((double)ny-1-(double)dat.H*0.5);
        int lx = std::ceil((double)nx/(double)dat.GX);
        int ly = std::ceil((double)dat.H*0.5/(double)dat.GY);
        // std::cout<<"starty "<<starty<<std::endl;
        // std::cout<<dat.GX<<" "<<dat.GY<<std::endl;

        for(int igx=0; igx< dat.GX; ++igx)
        for(int igy=0; igy< dat.GY; ++igy)
        {
            GridCount[igx*(dat.GX-1)+igy] = 0;
            int ixs = std::max(0,igx*lx);
            int ixe = std::min((double)nx,(double)(igx+1)*lx);
            int iys = std::max((double)starty,(double)igy*ly+starty);
            int iye = std::min((double)ny,(double)(igy+1)*ly+starty);
            // std::cout<<"x "<<ixs<<" "<<ixe<<std::endl;
            // std::cout<<"y "<<iys<<" "<<iye<<std::endl;
            for(int ix=ixs; ix<ixe; ++ix)
            for(int iy=iys; iy<iye; ++iy)
            {
                if(dom.Check[ix][iy][0]>0)
                {
                    GridCount[igx*(dat.GX-1)+igy] += 1;
                }
                
            }
            // std::cout<<igx<<" "<<igy<<" "<<igx*(dat.GX-1)+igy<<" "<<GridCount[igx*(dat.GX-1)+igy]<<std::endl;
        }

        for(int i=ipp; i<dat.pinit; ++i)
        {
        
            
            int index = std::distance(GridCount,std::min_element(GridCount,GridCount+dat.GX*dat.GY));
            // std::cout<<index<<std::endl;
            // for(int i=0; i<dat.GX*dat.GY; ++i)
            // std::cout<<GridCount[i]<<std::endl;

            int igy = index%(dat.GX-1);
            int igx = std::floor(index/(dat.GX-1));
            double ll = dat.R*1.2;
            int ixs = std::max(ll,(double)igx*lx);
            int ixe = std::min((double)nx-ll,(double)(igx+1)*lx);
            int iys = std::max((double)starty+ll,(double)igy*ly+starty);
            int iye = std::min((double)ny-ll,(double)(igy+1)*ly+starty);
            // std::cout<<igx<<" x "<<ixs<<" "<<ixe<<std::endl;
            // std::cout<<igy<<" y "<<iys<<" "<<iye<<std::endl;
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
                    if(Norm(pos - dom.Particles[ip].X)<3.0*dat.R || Norm(pos - dom.GhostParticles[ip].X)<3.0*dat.R)
                    {
                        overlap = true;
                        N += 1;
                        break;
                    }
                }
                if(!overlap) break;
            }
            if(N>1000) throw new Fatal("Cannot find!!!!");
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
            v = dom.Vel[ipx][ipy][0]; 
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
//			dom.Particles[ip].bbeta = 0.3;
//			dom.Particles[ip].epsilon = 2.05769e-20/((ratiol/ratiot)*(ratiol/ratiot));
//			dom.Particles[ip].s = 200e-9/ratiol;
//			dom.Particles[ip].Lc = 100e-9/ratiol;
//			dom.Particles[ip].l = 3.04e-10/ratiol;
            // dom.Particles.back().VdwCutoff = std::sqrt(dom.Particles.back().A/(12.0*dom.Particles.back().Z*dom.Particles.back().kappa));
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
    double py = dat.sy;
    double H = dat.H;
    std::cout<<"py = "<<py<<" H = "<<H<<std::endl;
    std::cout<<"max vel "<<Norm(dat.g0)/(2.0*dat.nu)*H*H<<std::endl;
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
            // double uy = dat.g/(2.0*dat.nu)*(H*(yy-py) - (yy-py)*(yy-py)); 
            double uy = Norm(dat.g0)/(2.0*dat.nu)*(-1.0*((Y-H)*(Y-H)-H*H)); 
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
    
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for(size_t ip=0; ip<dom.Particles.size(); ++ip)
    {
        
		if(dom.Particles[ip].X(1)>ny-1-dat.H*0.5)
        {
            #pragma omp atomic
            ++dat.pinit;
		}

	}
    std::cout<<"buffer particle num "<<dat.pinit<<std::endl;
}


int main (int argc, char **argv) try
{
    std::srand((unsigned)time(NULL));    
    

    char const *infilename = "test_swi_ini.txt";
    String fn;
    fn.Printf("%s",infilename);
    std::fstream ifile(fn.CStr(),std::ios::in);
    if(ifile.fail()) throw new Fatal("FILE WRONG");

    size_t Nproc = 12;
    if(argc>=2) Nproc = atoi(argv[1]); 
    double nu, Re, R, RR, ppl,pppl, pdx, pdy, ratiol, ratiot, ratiom, rho, rhos, Ga;
    double Kn, Gn, Kt, Mu, Eta, Beta, A, kappa, Z, D;
    int GX, GY,pinit;
    int Pnx, Pny, pnx, pny;
    ifile>>nu;  ifile.ignore(200,'\n');        
    ifile>>Re;  ifile.ignore(200,'\n');
    std::cout<<" nu "<<nu<<" Re "<<Re<<std::endl;
    ifile>>R;  ifile.ignore(200,'\n');
    ifile>>RR;  ifile.ignore(200,'\n');
    ifile>>ppl;  ifile.ignore(200,'\n');
    ifile>>pppl;  ifile.ignore(200,'\n');
    std::cout<<" R "<<R<<" RR "<<RR<<" ppl "<<ppl<<" pppl "<<pppl<<std::endl;
    ifile>>Pnx;  ifile.ignore(200,'\n');
    ifile>>Pny;  ifile.ignore(200,'\n');
    ifile>>pnx;  ifile.ignore(200,'\n');
    ifile>>pny;  ifile.ignore(200,'\n');
    ifile>>pdx;  ifile.ignore(200,'\n');
    ifile>>pdy;  ifile.ignore(200,'\n');
    std::cout<<" Pnx "<<Pnx<<" Pny "<<Pny<<" pnx "<<pnx<<" pny "<<pny<<" pdx "<<pdx<<" pdy "<<pdy<<std::endl;
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

    

    size_t nx = std::ceil(Pnx*(RR*2 + pdx));
    double pl = (double) nx/pnx;
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
    
    
    std::cout<<"H = "<<H<<std::endl;
    std::cout<<"sy = "<<sy<<std::endl;

    
    double vmax = nu*Re/H*1.5;
    std::cout<<"vmax "<<vmax<<std::endl;
    
    // double rho = 1.0;
    // double rhos = 2.7;
    // double Ga = 100.0;
    double gy = Ga*Ga*nu*nu/((8*R*R*R)*(rhos/rho-1));
    std::cout<<"gy = "<<gy<<std::endl;
    //nu = 1.0/30.0;
    std::cout<<nx<<" "<<ny<<" "<<nz<<std::endl;
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.pinit = pnx*pny*0.5;
    my_dat.GX = GX;
    my_dat.GY = GY;
    my_dat.nu = nu;
    my_dat.g = gy;
    my_dat.R = R;
    Vec3_t g0(2.0*nu*vmax/((double)H*(double)H),0.0,0.0);
    std::cout<<"gx = "<<g0<<std::endl;
    my_dat.H = H;
    my_dat.g0 = g0;
    my_dat.pnx = pnx;
    my_dat.pny = pny;
    my_dat.pl = pl;
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
        
        
        // for(int j=0; j<pny; ++j)
        // {
        //     py = sy + j*pl;
        //     for(int i=0; i<pnx; ++i)
        //     {
        //         Vec3_t dxr(random(-1.0*R,1.0*R),random(-1.0*R,1.0*R),0.0);
        //         // Vec3_t dxr(0.0,0.0,0.0);
        //         px = 0.5*pl+i*pl;
        //         pos = px, py , 0;
        //         pos = pos+dxr;
        //         v = Norm(g0)/(2.0*nu)*(-1.0*((pos(1)-sy-H)*(pos(1)-sy-H)-H*H)), 0,0;
        //         dom.Particles.push_back(DEM::Disk(-pnum, pos, v, w, rhos, R, dom.dtdem));
        //         // std::cout<<pos(0)<<" "<<pos(1)<<std::endl;
        //         pnum++;
        //     }
        // }

        
        
        
            
            for(int i=0; i< pnx*pny; ++i)
            {
                double x;
                double y;
                ifile>>x>>y;
                pos = x,y+sy,0;
                v = Norm(g0)/(2.0*nu)*(-1.0*((pos(1)-sy-H)*(pos(1)-sy-H)-H*H)), 0,0;
                dom.Particles.push_back(DEM::Disk(-pnum, pos, v, w, rhos, R, dom.dtdem));
                pnum++;
                // std::cout<<x<<" "<<y<<std::endl;

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
        // dom.Particles[ip].A = 2e-20/((ratiol/ratiot)*(ratiol/ratiot));
        // dom.Particles[ip].kappa = 1e8*ratiol;
        // dom.Particles[ip].Z = 1e-11*ratiot*ratiot/ratiol;
        // dom.Particles[ip].bbeta = 0.3;
        // dom.Particles[ip].epsilon = 2.05769e-20/((ratiol/ratiot)*(ratiol/ratiot));
        // dom.Particles[ip].s = 200e-9/ratiol;
        // dom.Particles[ip].Lc = 100e-9/ratiol;
        // dom.Particles[ip].l = 3.04e-10/ratiol;
        //dom.Particles[ip].VdwCutoff = std::sqrt(dom.Particles[ip].A/(12.0*dom.Particles[ip].Z*dom.Particles[ip].kappa));
        dom.Particles[ip].VdwCutoff = 1.57e-10/ratiol;
        dom.Particles[ip].D = my_dat.D;
        // dom.Particles[ip].R = R;
        dom.Particles[ip].M = M_PI*dom.Particles[ip].R*dom.Particles[ip].R*rhos;
        dom.Particles[ip].Ff =  0.0, -M_PI*dom.Particles[ip].R*dom.Particles[ip].R*(rhos/rho-1)*my_dat.g,0.0;
        if(dom.Particles[ip].IsFree())
        {
            
            dom.Particles[ip].Rh = 0.8*R;

        }else{
            dom.Particles[ip].Rh = RR;
        }
        


    }
    

    

    double Tf = 1e5;
    dom.IsF = true;    
    double dtout = 1e2;
    dom.Box = 0.0, nx-1, 0.0;
    dom.modexy = 0;
    //solving
    dom.SolveIBM( Tf, dtout, "test_swi1", Setup, NULL);
    
}MECHSYS_CATCH