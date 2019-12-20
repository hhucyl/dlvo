clear
clc
prefix = {'/media/user/My Book/'};
middle = {'10_f'};
num = [0:999];
H = 377;
for i = 1:numel(num)
    name = strcat(prefix,middle(1),'/test_swi1_',num2str(num(i),'%04d'),'.h5');
    nx = double(h5read(char(name),'/Nx'));
    ny = double(h5read(char(name),'/Ny'));
    p = h5read(char(name),'/Pposition');
    pr = h5read(char(name),'/PR');
    v  = h5read(char(name),'/Velocity_0');
    vx = reshape(v(1:3:end-2),[nx,ny]);
    vy = reshape(v(2:3:end-1),[nx,ny]);
    np = numel(p)/3;
    px = p(1:3:end-2);
    py = p(2:3:end-1);
    ppr = pr(1:np)-2;
    vv = sqrt(vx.^2+vy.^2);
    pcolor(vv')
    caxis([0 0.04])
    grid off
    shading interp
    hold on
    viscircles([px,py],ppr-2);
    hold on    
    axis equal
    xlim([0-8 nx+8])
    ylim([0 ny])
    plot([0 nx], [ny-0.5*H,ny-0.5*H],'k')
    title(i)

    drawnow
    
    fn = strcat(num2str(i),'.jpg');
    saveas(gcf,char(fn))
    clf
end