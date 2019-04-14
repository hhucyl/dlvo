function F = Fvdw(A,R1,R2,x)
    F = abs(A./(6.0.*x.^2).*(R1*R2/(R1+R2)));
end