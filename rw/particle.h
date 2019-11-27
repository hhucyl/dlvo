#ifndef RW_PARTICLE_H
#define RW_PARTICLE_H

#include <mechsys/linalg/matvec.h>
#include<random>

namespace RW
{
class Particle
{
    public:
    Particle(Vec3_t &X, double dm);
    Vec3_t X;
    Vec3_t Xb;
    Vec3_t er; //for adsorption
    // Vec3_t O; //for adsorption
    // double Rh; //for adsorption
    double Dm;

    bool AD; //whether adsorped by the solid surface
    // double Pa; //adsorption probability
    // double Pd; //desorption probability

    int ip;
    void Move(std::vector<Vec3_t> &VV, std::vector<int> &idx, double dt);
    void Move1(Vec3_t &V, double dt);
    void Move2(Vec3_t &V, double dt);
    void Leave(int modexy, Vec3_t &Box);
    void LeaveReflect(int modexy1, Vec3_t &Box1);
    void Reflect(Vec3_t &C, double R);
    void FindIntersectV(Vec3_t &C, Vec3_t &V, double R, Vec3_t &X, Vec3_t &Xb, Vec3_t &Xi); //with v
    void FindIntersect1(Vec3_t &C, double R, Vec3_t &X, Vec3_t &Xb, Vec3_t &Xi); //x in. xb out
    void FindIntersect2(Vec3_t &C, double R, Vec3_t &X, Vec3_t &Xb, Vec3_t &Xi); //x in. xb in

    void Adsorption(double Pa, int iip);
    void Desorption(double Pd);

    //random
    
    std::mt19937 gen;
};
inline Particle::Particle(Vec3_t &X0, double dm): gen(std::random_device()())
{
    X = X0;
    Xb = X;
    Dm = dm;
    AD = false;
    ip = -1;
    // Pa = -1.0;
    // Pd = -1.0;
    // O = 0,0,0;
    // er = 0,0,0;
    // Rh = 0,0,0;
}

inline void Particle::Adsorption(double Pa, int iip)
{
    // std::cout<<Pa<<std::endl;
    if(Pa>-0.5)
    {
        std::uniform_real_distribution<double> uniform(0,1);
        if(uniform(gen)<Pa)
        {
            AD = true;
            ip = iip;
        }
    }
    
}

inline void Particle::Desorption(double Pd)
{
    if(AD && Pd>-0.5)
    {
        std::uniform_real_distribution<double> uniform(0,1);
        if(uniform(gen)<Pd)
        {
            AD = false;
            ip = -1;
            // X = O + Rh*er;
            // std::cout<<"X "<<X<<" O "<<O(0)<<" "<<O(1)<<" "<<O(2)<<" Rh "<<Rh<<" er "<<er<<std::endl;
            // Xb = X;
        }
    }
    
}

inline void Particle::Move(std::vector<Vec3_t> &VV, std::vector<int> &idx, double dt)
{
    //VV 0=>11 1=>21 2=>12 3=>22
    //*12   *22
    //*11   *21
    //idx 0=>x1 1=>x2 2=>y1 3=>y2
    double x1 = (double) idx[0];
    double x2 = (double) idx[1];
    double y1 = (double) idx[2];
    double y2 = (double) idx[3];
    double x = X(0);
    double y = X(1);
    Vec3_t V(0,0,0);
    V = 1.0/((x2-x1)*(y2-y1))*(VV[0]*(x2-x)*(y2-y)+VV[1]*(x-x1)*(y2-y)+VV[2]*(x2-x)*(y-y1)+VV[3]*(x-x1)*(y-y1));
    // std::cout<<VV[0]<<" "<<VV[1]<<" "<<VV[2]<<" "<<VV[3]<<std::endl;
    // std::cout<<x1<<" "<<x2<<" "<<y1<<" "<<y2<<std::endl;
    // std::cout<<V(0)<<"    "<<V(1)<<std::endl;
    Xb = X;
    std::normal_distribution<double> normal(0,1);
    Vec3_t e(normal(gen),normal(gen),0.);
    X = X+V*dt+std::sqrt(2*Dm*dt)*e;    
}

inline void Particle::Move1(Vec3_t &V, double dt)
{
    std::normal_distribution<double> normal(0,1);
    Vec3_t e(normal(gen),normal(gen),0.);
    X = X+V*dt+std::sqrt(2*Dm*dt)*e;    
}

inline void Particle::Move2(Vec3_t &V, double dt)
{
    X = X+V*dt;    
}

inline void Particle::Leave(int modexy, Vec3_t &Box)
{
    if(1000*X(modexy)<1000*Box(0))
    {   
        double dist = std::fabs(X(modexy)-Box(0));
        X(modexy) = Box(1)-dist+1.0;
        double distb = Xb(modexy)-Box(0);
        Xb(modexy) = Box(1)+distb+1.0;   

    }
    if(X(modexy)>Box(1))
    {
        double dist = std::fabs(X(modexy)-Box(1));
        X(modexy) = Box(0)+dist-1.0;
        double distb = Xb(modexy)-Box(1);
        Xb(modexy) = Box(0)+distb-1.0;

    }
    
}

inline void Particle::LeaveReflect(int modexy1, Vec3_t &Box1)
{
    if(X(modexy1)<Box1(0))
    {
        X(modexy1) = Box1(0) + std::abs(X(modexy1) - Box1(0));
    }
    if(X(modexy1)>Box1(1))
    {
        X(modexy1) = Box1(1) - std::abs(X(modexy1) - Box1(1));
    }
}

inline void Particle::FindIntersectV(Vec3_t &C, Vec3_t &V, double R, Vec3_t &X, Vec3_t &Xb, Vec3_t &Xi)
{
    double a = X(0)-C(0);
    double b = Xb(0)-X(0)+V(0);
    double c = X(1)-C(1);
    double d = Xb(1)-X(1)+V(1);
    double AA = b*b + d*d;
    double BB = 2*(a*b+c*d);
    double CC = a*a + c*c - R*R;
    double Delta = BB*BB-4.0*AA*CC;
    if(Delta>0)
    {
        double q=0;
        double q1 = (-BB+std::sqrt(Delta))/(2.0*AA);
        double q2 = (-BB-std::sqrt(Delta))/(2.0*AA);
        bool flag1 = q1>=0 && q1-1<1e-12;
        bool flag2 = q2>=0 && q2-1<1e-12;
        if(flag1)
        {
            q = q1;
        }else{
            if(flag2)
            {
                q = q2;
            }else{
                q = std::max(q1,q2);
                // std::cout<<q1<<" "<<q2<<" "<<X<<" "<<Xb<<" "<<C<<" "<<V<<" "<<R<<std::endl;
                // throw new Fatal("ERROR IN RWPARTICLE FindIntersectV!!!");
            }
        }
        // std::cout<<"q = "<<q<<std::endl;
        Xi = X(0)+q*(Xb(0)-X(0)),X(1)+q*(Xb(1)-X(1)),0.0;
    }else{
        throw new Fatal("WRONG IN DELTA in FindIntersectV!!!!!");    
    }

}

inline void Particle::FindIntersect1(Vec3_t &C, double R, Vec3_t &X, Vec3_t &Xb, Vec3_t &Xi)
{
    double L1 = Norm(X-Xb);
    double L2 = Norm(X-C);
    double Dot = (X(0)-C(0))*(Xb(0)-X(0))+(X(1)-C(1))*(Xb(1)-X(1));
    double AA = L1*L1;
    double BB = 2.0*Dot;
    double CC = L2*L2 - R*R;
    double Delta = BB*BB-4.0*AA*CC;
    if(Delta>0)
    {
        double q;
        double q1 = (-BB+std::sqrt(Delta))/(2.0*AA);
        double q2 = (-BB-std::sqrt(Delta))/(2.0*AA);
        bool flag1 = q1>=0 && q1-1<1e-12;
        bool flag2 = q2>=0 && q2-1<1e-12;
        if(flag1)
        {
            q = q1;
        }else{
            if(flag2)
            {
                q = q2;
            }else{
                std::cout<<q2<<" "<<X<<" "<<Xb<<" "<<"ERROR IN RWPARTICLE REFLECT!!!"<<std::endl;
            }
        }
        // std::cout<<"q = "<<q<<std::endl;
        Xi = X(0)+q*(Xb(0)-X(0)),Xb(1)+q*(Xb(1)-X(1)),0.0;
    }else{
        throw new Fatal("WRONG IN DELTA!!!!!");    
    }
}

inline void Particle::FindIntersect2(Vec3_t &C, double R, Vec3_t &X, Vec3_t &Xb, Vec3_t &Xi)
{
    double L1 = Norm(X-Xb);
    double L2 = Norm(X-C);
    double Dot = (X(0)-C(0))*(Xb(0)-X(0))+(X(1)-C(1))*(Xb(1)-X(1));
    double AA = L1*L1;
    double BB = 2.0*Dot;
    double CC = L2*L2 - R*R;
    double Delta = BB*BB-4.0*AA*CC;
    if(Delta>0)
    {
        double q;
        double q1 = (-BB+std::sqrt(Delta))/(2.0*AA);
        double q2 = (-BB-std::sqrt(Delta))/(2.0*AA);
        bool flag1 = q1>=0;
        bool flag2 = q2>=0;
        if(flag1)
        {
            q = q1;
        }else{
            if(flag2)
            {
                q = q2;
            }else{
                std::cout<<q2<<" "<<X<<" "<<Xb<<" "<<"ERROR IN RWPARTICLE REFLECT!!!"<<std::endl;
            }
        }
        // std::cout<<"q = "<<q<<std::endl;
        Xi = Xb(0)+q*(X(0)-Xb(0)),Xb(1)+q*(X(1)-Xb(1)),0.0;
    }else{
        throw new Fatal("WRONG IN DELTA!!!!!");    
    }
}

inline void Particle::Reflect(Vec3_t &C, double R)
{
    double L1 = Norm(X-Xb);
    double L2 = Norm(Xb-C);
    // std::cout<<"L1 = "<<L1<<" L2 ="<<L2<<std::endl;
    double AA = L1*L1;
    double BB = -2.0*((Xb(1)-C(1))*(Xb(1)-X(1))+(Xb(0)-C(0))*(Xb(0)-X(0)));
    double CC = L2*L2 - R*R;
    double Delta = BB*BB-4.0*AA*CC;
    // std::cout<<"Delta = "<<Delta<<std::endl;
    
    if(Delta>0)
    {
        double q;
        double q1 = (-BB+std::sqrt(Delta))/(2.0*AA);
        double q2 = (-BB-std::sqrt(Delta))/(2.0*AA);
        bool flag1 = q1>=0 && q1-1<1e-9;
        bool flag2 = q2>=0 && q2-1<1e-9;
        if(flag1)
        {
            q = q1;
        }else{
            if(flag2)
            {
                q = q2;
            }else{
                std::cout<<q2<<" "<<X<<" "<<Xb<<" "<<"ERROR IN RWPARTICLE REFLECT!!!"<<std::endl;
            }
        }
        // std::cout<<"q = "<<q<<std::endl;
        Vec3_t Xi(Xb(0)-q*(Xb(0)-X(0)),Xb(1)-q*(Xb(1)-X(1)),0.0);
        // std::cout<<"D = "<<Xi<<std::endl;
        double y = -((X(0)-2.*Xi(0))*(C(0)-Xi(0))*(C(1)-Xi(1)) + (X(1)-2.*Xi(1))*(C(1)-Xi(1))*(C(1)-Xi(1)) + X(0)*(C(0)-Xi(0))*(C(1)-Xi(1)) - X(1)*(C(0)-Xi(0))*(C(0)-Xi(0)))/(R*R);
        double x;
        if(std::fabs(C(1)-Xi(1))>0)
        {
            x = ((C(0)-Xi(0))*(y-X(1)))/(C(1)-Xi(1)) + X(0);
        }else{
            x = -((C(1)-Xi(1))*y + (X(0)-2.*Xi(0))*(C(0)-Xi(0)) + (X(1)-2.*Xi(1))*(C(1)-Xi(1)))/(C(0)-Xi(0));
        }
        // std::cout<<x<<" "<<y<<std::endl;
        X = x,y,0.0;
    }
    
}


}


#endif