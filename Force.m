clear
clc
R = 10;
ratiol = 10e-6/(2*R);
ratiot = 1/1;
A = 2e-20/(ratiol/ratiot)^2;
kappa = 1e9*ratiol;
Z = 1e-11*ratiot*ratiot/ratiol;
x = 5e-4:0.0001:0.02;
fvdw = Fvdw(A,R,R,x);
fe = Fe(kappa,Z,R,R,x);
subplot(121)
plot(x,fvdw,'r')
hold on
plot(x,fe,'b')
legend('Vdw','Elec')
xlabel('\itdx(\rm Unit of Lattice Boltzmann)')
ylabel('\itF(\rm Unit of Lattice Boltzmann)')
subplot(122)
F = fvdw-fe;
kkk = find(F<0);
plot(x,F,'r')
hold on
plot(x(kkk),F(kkk),'b')
grid on
plot([5e-4,5e-4],[0,Fvdw(A,R,R,5e-4)-Fe(kappa,Z,R,R,5e-4)],'k')
plot([5e-3,5e-3],[0,Fvdw(A,R,R,5e-3)-Fe(kappa,Z,R,R,5e-3)],'k')
xlabel('\itdx(\rm Unit of Lattice Boltzmann)')
ylabel('\itF(\rm Unit of Lattice Boltzmann)')
FF = fvdw'-fe';
FFe = fe'.*ratiol/ratiot/ratiot;
FFvdw = fvdw'.*ratiol/ratiot/ratiot;
xx = x'.*ratiol;
