clear
clc
prefix = {'/media/user/My Book/'};
middle = {'10_f',...
    '10_f_e_5','10_f_e_25','10_f_e_35',...
    '10_f_mu_0.0','10_f_mu_0.8',...
    '10_f_r_0.6','10_f_r_0.8'};
% middle = {'10_f'};
num = [0:999];
ppy = 400;
for ii = 1:numel(middle)
    for i = 1:numel(num)
        name = strcat(prefix,middle(ii),'/test_swi1_',num2str(num(i),'%04d'),'.h5');
        nx = double(h5read(char(name),'/Nx'));
        ny = double(h5read(char(name),'/Ny'));
        v  = h5read(char(name),'/Velocity_0');
        vx = reshape(v(1:3:end-2),[nx,ny]);
        vy = reshape(v(2:3:end-1),[nx,ny]);
        vyy = vy(:,ppy);
        F(i,ii) = sum(vyy(vyy<0));
        i
    end
    ii
end