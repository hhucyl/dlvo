clear
clc
prefix = {'/media/user/My Book/'};
middle = {'10_f',...
    '10_f_e_5','10_f_e_25','10_f_e_35',...
    '10_f_mu_0.0','10_f_mu_0.8',...
    '10_f_r_0.6','10_f_r_0.8'};
for i = 1:numel(middle)
    name = strcat(prefix,middle(i),'/test_swi1_0999.h5');
    nx = double(h5read(char(name),'/Nx'));
    ny = double(h5read(char(name),'/Ny'));
    p = h5read(char(name),'/Pposition');
    pr = h5read(char(name),'/PR');
    np = numel(p)/6;
    px = p(1:3:3*np-2);
    py = p(2:3:3*np-1);
    ppr = pr(1:np);
    figure(i)
    viscircles([px,py],ppr)
    xlim([0 nx])
    axis equal
    fn = strcat(middle(i),'.jpg');
    saveas(gcf,char(fn))
    
end