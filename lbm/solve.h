#ifndef LBM_SOLVE_H
#define LBM_SOLVE_H

inline void Domain::Solve(double Tf, double dtout, char const * TheFileKey, ptDFun_t ptSetup, ptDFun_t ptReport)
{
    StartSolve();
    std::cout<<"Box "<<Box<<std::endl;
    std::cout<<"modexy "<<modexy<<std::endl;
    std::cout<<"dt of LBM "<<dt<<" dt of DEM "<<dtdem<<std::endl;
    double tout = 0;
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1   , TERM_RST);
    GhostParticles.assign(Particles.begin(), Particles.end()); 
    while(Time<Tf)
    {
        if (Time>=tout)
        {
            
            String fn;
            fn.Printf("%s_%04d", TheFileKey, idx_out);
            
            WriteXDMF(fn.CStr());
            idx_out++;
            
            if (ptReport!=NULL) (*ptReport) ((*this), UserData); 
            
            

            tout += dtout;
        }
        

        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData); 
        


        for(int i=0; i<std::floor(dt/dtdem); ++i)
        {
            // std::cout<<i<<std::endl;
            bool flag = i==0 || i==(std::floor(dt/dtdem)-1);
        
            if(flag){
                SetZero();
            }
            //set added force and check leave particles
            LeaveAndForcedForce();

            GhostParticles.clear();
            GhostParticles.assign(Particles.begin(), Particles.end()); 
        
            GhostPeriodic();
        
         
            //set fluid force
            if(flag){
                AddDisksG();
            }

            //update particles contact
            if(flag){
                UpdateParticlesContacts();
            }
            
            UpdateParticlePairForce();
        
            //move
            MoveParticles();
        }
        

        //collide and streaming
        (this->*ptr2collide)();
        Stream();
        BounceBack(false);
        CalcProps();

        Time += 1;
    }
    printf("%s  Final CPU time       = %s\n",TERM_CLR2, TERM_RST);
}

inline void Domain::SolveIBM(double Tf, double dtout, char const * TheFileKey, ptDFun_t ptSetup, ptDFun_t ptReport)
{
    StartSolve();
    std::cout<<"Box "<<Box<<std::endl;
    std::cout<<"modexy "<<modexy<<std::endl;
    std::cout<<"dt of LBM "<<dt<<" dt of DEM "<<dtdem<<std::endl;
    double tout = 0;
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1   , TERM_RST);
    GhostParticles.assign(Particles.begin(), Particles.end()); 
    while(Time<Tf)
    {
        if (Time>=tout)
        {
            
            String fn;
            fn.Printf("%s_%04d", TheFileKey, idx_out);
            
            WriteXDMF(fn.CStr());
            idx_out++;
            
            if (ptReport!=NULL) (*ptReport) ((*this), UserData); 
            
            

            tout += dtout;
        }
        
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData); 
        

        


        for(int i=0; i<std::floor(dt/dtdem); ++i)
        {
            // std::cout<<i<<std::endl;
            bool flag = i==0 || i==(std::floor(dt/dtdem)-1);
        
            if(flag){
                SetZero();
            }
            //set added force and check leave particles
            LeaveAndForcedForce();

            GhostParticles.clear();
            GhostParticles.assign(Particles.begin(), Particles.end()); 
        
            GhostPeriodic();
            
        
         
            //set fluid force
            if(flag){
                AddDisksIBM();
            }

            //update particles contact
            if(flag){
                UpdateParticlesContacts();
                // UpdateParticlesContactsIBM();
            }
            
            UpdateParticlePairForce();
        
            //move
            MoveParticles();
            // GhostParticles.clear();
            // GhostParticles.assign(Particles.begin(), Particles.end());
            // GhostPeriodic(); 
        }
        

        //collide and streaming
        CollideMRTIBM();
        Stream();
        BounceBack(false);
        CalcProps();
        // std::cout<<std::boolalpha<<GhostParticles[0].Ghost<<std::endl;

        //trace particle
        If(IsRW) rwsolve_sub(dt);

        Time += 1;
    }
    printf("%s  Final CPU time       = %s\n",TERM_CLR2, TERM_RST);
}

inline void Domain::SolveIBM1(double Tf, double dtout, char const * TheFileKey, ptDFun_t ptSetup, int id, size_t xy, ptDFun_t ptReport)
{
    StartSolve();
    std::cout<<"Box "<<Box<<std::endl;
    std::cout<<"modexy "<<modexy<<std::endl;
    std::cout<<"dt of LBM "<<dt<<" dt of DEM "<<dtdem<<std::endl;
    double tout = 0;
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1   , TERM_RST);
    GhostParticles.assign(Particles.begin(), Particles.end()); 
    while(Time<Tf)
    {
        if (Time>=tout)
        {
            
            String fn;
            fn.Printf("%s_%04d", TheFileKey, idx_out);
            
            WriteXDMF(fn.CStr());
            idx_out++;
            
            if (ptReport!=NULL) (*ptReport) ((*this), UserData); 
            
            

            tout += dtout;
        }
        
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData); 
        

        


        for(int i=0; i<std::floor(dt/dtdem); ++i)
        {
            // std::cout<<i<<std::endl;
            bool flag = i==0 || i==(std::floor(dt/dtdem)-1);
        
            if(flag){
                SetZero();
            }
            //set added force and check leave particles
            LeaveAndForcedForce();

            GhostParticles.clear();
            GhostParticles.assign(Particles.begin(), Particles.end()); 
        
            GhostPeriodic();
            
        
         
            //set fluid force
            if(flag){
                AddDisksIBM();
            }

            //update particles contact
            if(flag){
                UpdateParticlesContacts();
                // UpdateParticlesContactsVL();
            }
            
            UpdateParticlePairForce();
        
            //move
            MoveParticlesWithRestriction(id,xy);

        }
        

        //collide and streaming
        CollideMRTIBM();
        Stream();
        BounceBack(false);
        CalcProps();
        Time += 1;
    }
    printf("%s  Final CPU time       = %s\n",TERM_CLR2, TERM_RST);
}


inline void Domain::SolveRW(double Tf, double dtout, char const * TheFileKey, ptDFun_t ptSetup, ptDFun_t ptReport)
{
    StartSolve();
    std::cout<<"Box "<<Box<<std::endl;
    std::cout<<"modexy "<<modexy<<std::endl;
    std::cout<<"dt of LBM "<<dt<<" dt of DEM "<<dtdem<<std::endl;
    double tout = 0;
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1   , TERM_RST);
    GhostParticles.assign(Particles.begin(), Particles.end()); 
    while(Time<Tf)
    {
        if (Time>=tout)
        {
            
            String fn;
            fn.Printf("%s_%04d", TheFileKey, idx_out);
            
            WriteXDMF(fn.CStr());
            idx_out++;
            
            if (ptReport!=NULL) (*ptReport) ((*this), UserData); 
            
            

            tout += dtout;
        }
        

        


        for(int i=0; i<std::floor(dt/dtdem); ++i)
        {
            // std::cout<<i<<std::endl;
            bool flag = i==0 || i==(std::floor(dt/dtdem)-1);
        
            if(flag){
                SetZero();
            }
            //set added force and check leave particles
            LeaveAndForcedForce();

            GhostParticles.clear();
            GhostParticles.assign(Particles.begin(), Particles.end()); 
        
            GhostPeriodic();
        
         
            //set fluid force
            if(flag){
                AddDisksG();
            }

            //update particles contact
            if(flag){
                UpdateParticlesContacts();
            }
            
            UpdateParticlePairForce();
        
            //move
            MoveParticles();
        }
        

        //collide and streaming
        (this->*ptr2collide)();
        Stream();
        BounceBack(false);
        CalcProps();

        //trace particle
        rwsolve_sub(dt);

        Time += 1;
    }
    printf("%s  Final CPU time       = %s\n",TERM_CLR2, TERM_RST);
}

inline void Domain::rwsolve_sub(double dt)
{
    // std::cout<<1<<std::endl;
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(int i=0;i<(int) RWParticles.size();++i)
    {
        RW::Particle *RWP = &RWParticles[i]; 
        int x1 = std::floor(RWP->X(0));
        int y1 = std::floor(RWP->X(1));
        int x2 = x1+1;
        int y2 = y1+1;
        if(x1<0 || y1<0 || x2>Ndim(0)-1 || y2>Ndim(1)-1)
        {
            RWP->Leave(modexy,Box);
            x1 = std::floor(RWP->X(0));
            y1 = std::floor(RWP->X(1));
            x2 = x1+1;
            y2 = y1+1;
        }
        std::vector<Vec3_t> VV{Vel[x1][y1][0],Vel[x2][y1][0],Vel[x1][y2][0],Vel[x2][y2][0]};
        std::vector<int> idx{x1,x2,y1,y2};
        RWP->Move(VV,idx,dt);
        RWP->Leave(modexy,Box);
        int ix = std::round(RWP->X(0));
        int iy = std::round(RWP->X(1));
        if(Check[ix][iy][0]>0) 
        {
            int ip = Check[ix][iy][0];
            DEM::Disk* Pa = &Particles[ip];
            DEM::Disk* GPa = &GhostParticles[ip];
            if(Norm(Pa->X-RWP->X)<=Pa->Rh) 
            {
                Vec3_t Xi(0,0,0);
                RWP->FindIntersectV(Pa->X,Pa->V,Pa->Rh,RWP->X,RWP->Xb,Xi);
                double dt = Norm(Xi-RWP->Xb)/Norm(RWP->X-RWP->Xb);
                // Vec3_t tmp;
                // Rotation(Pa->W,Pa->Q,tmp);     
                // Vec3_t VelPt   = Pa->V + cross(tmp,Xi);
                Vec3_t VelPt   = Pa->V;
                bool IsIn = false;
                do{
                    IsIn = false;
                    RWP->Move(VelPt,dt);
                    RWP->Leave(modexy,Box);
                    if(Norm(Pa->X-RWP->X)<=Pa->Rh && !IsIn)
                    {
                        IsIn = true;
                    }
                    if(Norm(GPa->X-RWP->X)<=GPa->Rh && !IsIn)
                    {
                        IsIn = true;
                    }
                    int ix = std::round(RWP->X(0));
                    int iy = std::round(RWP->X(1));
                    std::pair<int,int> temp(ix,iy);
                    if(GridPair.count(temp)>0 && !IsIn)
                    {
                        int ip1 = GridPair[temp].first;
                        int ip2 = GridPair[temp].second;
                        DEM::Disk *P1 = &Particles[ip1];
                        DEM::Disk *GP1 = &GhostParticles[ip1];
                        DEM::Disk *P2 = &Particles[ip2];
                        DEM::Disk *GP2 = &GhostParticles[ip2];
                        if(Norm(P1->X-RWP->X)<=P1->Rh && !IsIn)
                        {
                            IsIn = true;
                        }
                        if(Norm(GP1->X-RWP->X)<=GP1->Rh && !IsIn)
                        {
                            IsIn = true;
                        }
                        if(Norm(P2->X-RWP->X)<=P2->Rh && !IsIn)
                        {
                            IsIn = true;
                        }
                        if(Norm(GP2->X-RWP->X)<=GP2->Rh && !IsIn)
                        {
                            IsIn = true;
                        }
                         
                    }
                        
                    
                }while(IsIn)
            }

            if(Norm(GPa->X-RWP->X)<=GPa->Rh) 
            {
                Vec3_t Xi(0,0,0);
                RWP->FindIntersectV(GPa->X,GPa->V,GPa->Rh,RWP->X,RWP->Xb,Xi);
                double dt = Norm(Xi-RWP->Xb)/Norm(RWP->X-RWP->Xb);
                // Vec3_t tmp;
                // Rotation(GPa->W,GPa->Q,tmp);     
                // Vec3_t VelPt   = GPa->V + cross(tmp,Xi);
                Vec3_t VelPt   = GPa->V;
                bool IsIn = false;
                do{
                    IsIn = false;
                    RWP->Move(VelPt,dt);
                    RWP->Leave(modexy,Box);
                    if(Norm(Pa->X-RWP->X)<=Pa->Rh && !IsIn)
                    {
                        IsIn = true;
                    }
                    if(Norm(GPa->X-RWP->X)<=GPa->Rh && !IsIn)
                    {
                        IsIn = true;
                    }
                    int ix = std::round(RWP->X(0));
                    int iy = std::round(RWP->X(1));
                    std::pair<int,int> temp(ix,iy);
                    if(GridPair.count(temp)>0 && !IsIn)
                    {
                        int ip1 = GridPair[temp].first;
                        int ip2 = GridPair[temp].second;
                        DEM::Disk *P1 = &Particles[ip1];
                        DEM::Disk *GP1 = &GhostParticles[ip1];
                        DEM::Disk *P2 = &Particles[ip2];
                        DEM::Disk *GP2 = &GhostParticles[ip2];
                        if(Norm(P1->X-RWP->X)<=P1->Rh && !IsIn)
                        {
                            IsIn = true;
                        }
                        if(Norm(GP1->X-RWP->X)<=GP1->Rh && !IsIn)
                        {
                            IsIn = true;
                        }
                        if(Norm(P2->X-RWP->X)<=P2->Rh && !IsIn)
                        {
                            IsIn = true;
                        }
                        if(Norm(GP2->X-RWP->X)<=GP2->Rh && !IsIn)
                        {
                            IsIn = true;
                        }
                         
                    }
                        
                    
                }while(IsIn)
            }
        }
    }

}

inline void Domain::CheckInside()
{
    for(size_t i=0; i<RWParticles.size(); ++i)
    {
        RW::Particle *RWP = &RWParticles[i]; 
        for(size_t ip=0; ip<Particles.size(); ++ip)
        {
            DEM::Disk *Pa = &Particles[ip];
            if(Norm(Pa->X-RWP->X)<Pa->Rh)
            {
                RWParticles.erase(RWParticles.begin()+i);
                break;
            }
        }

    }
}

#endif