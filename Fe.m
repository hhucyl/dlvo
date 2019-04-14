function F = Fe(kappa,Z,R1,R2,x)
    F = abs(kappa*(R1*R2/(R1+R2))*Z.*exp(-kappa.*x));
end