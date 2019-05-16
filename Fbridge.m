function F = Fbridge(beta,R1,R2,epsilon,s,Lc,l,x)
    F = zeros(numel(x),1);
    kkk = find(x<Lc);
    F(kkk) = abs(beta*4.0*pi*(R1*R2/(R1+R2))*epsilon/(s*s).*(Lc-x(kkk))./l);
end