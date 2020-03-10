#ifndef LBM_SOLVE_RW_H
#define LBM_SOLVE_RW_H

inline void Domain::rwsolve_sub_nopar(double dt)
{
    if(Time<0.5) std::cout<<"--- rwsolve_sub_nopar ---"<<std::endl;
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(int i=0;i<(int) RWParticles.size();++i)
    {
        RW::Particle *RWP = &RWParticles[i];
        std::vector<int> idc{-1,-1};
        Vec3_t v0(0,0,0); 
        std::vector<Vec3_t> VV{v0,v0,v0,v0};
        std::vector<int> idx{-1,-1,-1,-1};
        // std::cout<<111<<std::endl;
        Pa2GridV(RWP,idx,VV);
        // std::cout<<222<<std::endl;
        RWP->Move(VV,idx,dt);
        RWP->Leave(modexy,Box);
        RWP->LeaveReflect(modexy1,Box1);
        Pa2Grid(RWP,idc);
        int ix = idc[0];
        int iy = idc[1];
        #pragma omp atomic
            Con[ix][iy][0] += 1;
    }
}

inline void Domain::rwsolve_sub(double dt)
{
    if(Time<0.5) std::cout<<"--- rwsolve_sub ---"<<std::endl;
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(int i=0;i<(int) RWParticles.size();++i)
    {
        RW::Particle *RWP = &RWParticles[i];
        std::vector<int> idc{-1,-1};
        if(RWP->AD)
        {
            int ip = RWP->ip;
            if(ip<0)
            {
                std::cout<<"ip = "<<ip<<std::endl;
                throw new Fatal("ERROR IN RWPARTICLE IP!!!");
            }
            DEM::Disk* Pa = &Particles[ip]; 
            Vec3_t VelPt = Pa->V;
            RWP->Move2(VelPt,dt);
            RWP->Leave(modexy,Box);
            RWP->LeaveReflect(modexy1,Box1);
            RWP->Desorption(Pa->Pd);
            if(!RWP->AD)
            {
                #pragma omp atomic
                    Pa->Alimit0 -= 1;
            }
        }

        if(!RWP->AD)
        {
            Vec3_t v0(0,0,0); 
            std::vector<Vec3_t> VV{v0,v0,v0,v0};
            std::vector<int> idx{-1,-1,-1,-1};
            // std::cout<<111<<std::endl;
            Pa2GridV(RWP,idx,VV);
            // std::cout<<222<<std::endl;
            RWP->Move(VV,idx,dt);
            RWP->Leave(modexy,Box);
            RWP->LeaveReflect(modexy1,Box1);

            
            Pa2Grid(RWP,idc);
            int ix = idc[0];
            int iy = idc[1];
            // std::cerr<<RWP->X(0)<<" "<<RWP->X(1)<<" "<<ix<<" "<<iy<<std::endl;
            // if(ix<0 || ix>Ndim(0)-1) std::cout<<"x "<<ix<<std::endl;
            // if(iy<0 || iy>Ndim(1)-1) std::cout<<"y "<<iy<<std::endl;
            if(Check[ix][iy][0][0]>-0.5) 
            {

                if(Check[ix][iy][0][1]>-0.5)
                {
                    int ip1 = Check[ix][iy][0][0];
                    int ip2 = Check[ix][iy][0][1];
                    DEM::Disk *P1 = &Particles[ip1];
                    DEM::Disk *GP1 = &GhostParticles[ip1];
                    if(Norm(P1->X-RWP->X)<=P1->Rh || Norm(GP1->X-RWP->X)<=P1->Rh ) 
                    {
                        SurfaceReaction(ip1, P1, GP1, RWP);
                        CheckInside(P1, GP1, RWP);
                    }
                    DEM::Disk *P2 = &Particles[ip2];
                    DEM::Disk *GP2 = &GhostParticles[ip2];
                    if(Norm(P2->X-RWP->X)<=P2->Rh || Norm(GP2->X-RWP->X)<=P2->Rh )
                    {
                        SurfaceReaction(ip2, P2, GP2, RWP);
                        CheckInside(P2, GP2, RWP);
                    }
                }else{
                    
                    // std::cout<<333<<std::endl;
                    int ip = Check[ix][iy][0][0];
                    DEM::Disk* Pa = &Particles[ip];
                    DEM::Disk* GPa = &GhostParticles[ip];
                    SurfaceReaction(ip, Pa, GPa, RWP);
                    CheckInside(Pa, GPa, RWP);
                }

                
            }
        }
        // std::cout<<444<<std::endl;
        Pa2Grid(RWP,idc);
        int ix = idc[0];
        int iy = idc[1];
        #pragma omp atomic
            Con[ix][iy][0] += 1;

    }

}



inline void Domain::SurfaceReaction(int ip, DEM::Disk *Pa, DEM::Disk *GPa, RW::Particle *RWP)
{
    bool flag1 = Norm(Pa->X-RWP->X)<=Pa->Rh;
    bool flag2 = Norm(GPa->X-RWP->X)<=Pa->Rh;
    if( flag1 || flag2 ) 
    {
        Vec3_t O(0,0,0);
        if(flag1)
        {
            O = Pa->X;
        }
        else
        {
            O = GPa->X;
        }

        if(Pa->Alimit0 < Pa->Alimit) RWP->Adsorption(Pa->Pa, ip);
        // RWP->Adsorption(Pa->Pa, ip);
        if(RWP->AD)
        {

            Vec3_t Xi(0,0,0);
            RWP->FindIntersectV(O,Pa->V,Pa->Rh,RWP->X,RWP->Xb,Xi);
            double q = Norm(Xi-O)/Norm(RWP->X-RWP->Xb);
            if(q>1.1)
            {   
                Xi = RWP->X;
            } 
            RWP->er = (Xi-O)/Norm(Xi-O);
            RWP->X = O + Pa->Rh*RWP->er;

            RWP->Leave(modexy,Box);
            RWP->LeaveReflect(modexy1,Box1);
            #pragma omp atomic
                Pa->Alimit0 += 1;

        }
    }
}

inline void Domain::CheckInside(DEM::Disk *Pa, DEM::Disk *GPa, RW::Particle *RWP)
{
    if((Norm(Pa->X-RWP->X)<=Pa->Rh || Norm(GPa->X-RWP->X)<=Pa->Rh) && !RWP->AD) 
    {

        Vec3_t Xi(0,0,0);
        Vec3_t O(0,0,0);
        Vec3_t Ob(0,0,0);
        if(Norm(Pa->X-RWP->X)<=Pa->Rh)
        {
            RWP->FindIntersectV(Pa->X,Pa->V,Pa->Rh,RWP->X,RWP->Xb,Xi);
            O = Pa->X;
            Ob = Pa->Xb;
        }else{
            RWP->FindIntersectV(GPa->X,GPa->V,GPa->Rh,RWP->X,RWP->Xb,Xi);
            O = GPa->X;
            Ob = GPa->Xb;
        }

        double dtt = Norm(Xi-RWP->Xb)/Norm(RWP->X-RWP->Xb);
        
        if(dtt>1.1)
        { 
            Xi = RWP->X;
            Vec3_t er = (Xi-O)/Norm(Xi-O);
            RWP->X = O + Pa->Rh*er;

            RWP->Leave(modexy,Box);
            RWP->LeaveReflect(modexy1,Box1);
            Xi = RWP->X;
        } 
        // std::cout<<dt<<std::endl;
        dtt = 1.0;
        // Vec3_t tmp;
        // Rotation(Pa->W,Pa->Q,tmp);     
        // Vec3_t VelPt   = Pa->V + cross(tmp,Xi);
        Vec3_t VelPt   = Pa->V;

        // Vec3_t v0(0,0,0); 
        // std::vector<Vec3_t> VV{v0,v0,v0,v0};
        // std::vector<int> idx{-1,-1,-1,-1};
        
        bool IsIn = false;
        int nn = 0;
        bool fflag = false;
        do{
            nn++;
            // std::cout<<"step "<<nn<<std::endl;
            IsIn = false;
            RWP->X = Xi;
            RWP->Move1(VelPt,dtt);
            RWP->Leave(modexy,Box);
            RWP->LeaveReflect(modexy1,Box1);
            
            if(Norm(Pa->X-RWP->X)<=Pa->Rh )
            {
                IsIn = true;
            }
            if(Norm(GPa->X-RWP->X)<=GPa->Rh )
            {
                IsIn = true;
            }
            std::vector<int> idc{-1,-1};
            Pa2Grid(RWP,idc);
            int ix = idc[0];
            int iy = idc[1];
            if(Check[ix][iy][0][1]>-0.5 )
            {
                int ip1 = Check[ix][iy][0][0];
                int ip2 = Check[ix][iy][0][1];
                DEM::Disk *P1 = &Particles[ip1];
                DEM::Disk *GP1 = &GhostParticles[ip1];
                DEM::Disk *P2 = &Particles[ip2];
                DEM::Disk *GP2 = &GhostParticles[ip2];
                if(Norm(P1->X-RWP->X)<=P1->Rh )
                {
                    IsIn = true;
                }
                if(Norm(GP1->X-RWP->X)<=GP1->Rh )
                {
                    IsIn = true;
                }
                if(Norm(P2->X-RWP->X)<=P2->Rh )
                {
                    IsIn = true;
                }
                if(Norm(GP2->X-RWP->X)<=GP2->Rh )
                {
                    IsIn = true;
                }
                    
            }

            if(nn>1e4)
            {
                fflag = true;
                break;
            }
            
        }while(IsIn);
        // std::cout<<nn<<std::endl;
        if(fflag)
        {
            // std::cout<<Xi<<std::endl;
            Vec3_t er = (Xi-O)/Norm(Xi-O);
            RWP->X = O + Pa->Rh*er;
            RWP->Xb = RWP->X;
            RWP->Leave(modexy,Box);
            RWP->LeaveReflect(modexy1,Box1);
            // throw new Fatal("N!!!");
        }
    }
    
}


inline void Domain::Pa2Grid(RW::Particle *RWP,std::vector<int> &idc)
{
    int xx = std::round(RWP->X(0));
    int yy = std::round(RWP->X(1));
    if(xx<0) xx = (int)Ndim(0) + xx;
    if(xx>(int)Ndim(0)-1) xx = xx - (int)Ndim(0);
    if(yy<0) yy = (int)Ndim(1) + yy;
    if(yy>(int)Ndim(0)-1) yy = yy - (int)Ndim(1);
    idc[0] = xx;
    idc[1] = yy;  
}
inline void Domain::Pa2GridV(RW::Particle *RWP, std::vector<int> &idx, std::vector<Vec3_t> &VV)
{
    //VV 0=>11 1=>21 2=>12 3=>22
    //*12   *22
    //*11   *21
    //idx 0=>x1 1=>x2 2=>y1 3=>y2
    //!! The other dimension of periodic dimension is considered as no flux boundary using specular reflection
    // the relative location is not change only the velocity was changed and idc
    int x1 = std::floor(RWP->X(0));
    int y1 = std::floor(RWP->X(1));
    int x2 = x1+1;
    int y2 = y1+1;
    idx[0] = x1;
    idx[1] = x2;
    idx[2] = y1;
    idx[3] = y2;
    if(modexy ==0)
    {
        if(x1<0)
        {
            int x11 = (int)Ndim(modexy)+x1;
            VV[0] = Vel[x11][y1][0];
            VV[1] = Vel[x2][y1][0];
            VV[2] = Vel[x11][y2][0];
            VV[3] = Vel[x2][y2][0];
        }else{
            if(x2>(int)Ndim(modexy)-1)
            {
                int x22 = x2 - (int)Ndim(modexy);
                VV[0] = Vel[x1][y1][0];
                VV[1] = Vel[x22][y1][0];
                VV[2] = Vel[x1][y2][0];
                VV[3] = Vel[x22][y2][0];
            }else{
                VV[0] = Vel[x1][y1][0];
                VV[1] = Vel[x2][y1][0];
                VV[2] = Vel[x1][y2][0];
                VV[3] = Vel[x2][y2][0];
            }
            
        }
        

    }else{
        if(y1<0)
        {
            int y11 = (int)Ndim(modexy)+y1;
            VV[0] = Vel[x1][y11][0];
            VV[1] = Vel[x2][y11][0];
            VV[2] = Vel[x1][y2][0];
            VV[3] = Vel[x2][y2][0];
        }else{
            if(y2>(int)Ndim(modexy)-1)
            {
                int y22 = y2 - (int)Ndim(modexy);
                VV[0] = Vel[x1][y1][0];
                VV[1] = Vel[x2][y1][0];
                VV[2] = Vel[x1][y22][0];
                VV[3] = Vel[x2][y22][0];
            }else{
                VV[0] = Vel[x1][y1][0];
                VV[1] = Vel[x2][y1][0];
                VV[2] = Vel[x1][y2][0];
                VV[3] = Vel[x2][y2][0];
            }
            
        }

    }
}


#endif