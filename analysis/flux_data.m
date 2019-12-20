clear
clc
prefix = {'/media/user/My Book/'};
middle = {'10_f',...
    '10_f',...
    '10_f_r_0.6','10_f_r_0.8'};
% middle = {'10_f'};
num = [50:999];
H = 377;
c = {'','r*','k*','b*'};
for ii = 2:4
    for i = 1:numel(num)
        name = strcat(prefix,middle(ii),'/test_swi1_',num2str(num(i),'%04d'),'.h5');
        nx = double(h5read(char(name),'/Nx'));
        ny = double(h5read(char(name),'/Ny'));
        v  = h5read(char(name),'/Velocity_0');
        vx = reshape(v(1:3:end-2),[nx,ny]);
        vy = reshape(v(2:3:end-1),[nx,ny]);
        
        sy = floor(ny-0.5*H);
        vyy = vy(:,sy);
        f = sum(vyy(vyy<0));
        F(i,ii) = abs(f);
        h(ii) = semilogy(i,F(i,ii),char(c(ii)));
        hold on
        drawnow
        
    end
    ii
end
xlabel('\itt')
ylabel('\itFlux')
legend(h(2:4),'\itr\rm = 0.3','\itr\rm = 0.6','\itr\rm = 0.8','location','southwest')
legend boxoff