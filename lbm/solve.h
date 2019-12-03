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
    std::cout<<"Box1 "<<Box1<<std::endl;
    std::cout<<"modexy1 "<<modexy1<<std::endl;
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
        if(IsRW) rwsolve_sub(dt);

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
    std::cout<<"Box1 "<<Box<<std::endl;
    std::cout<<"modexy1 "<<modexy<<std::endl;
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



#endif