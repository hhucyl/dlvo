clear
clc
prefix = {'/media/user/My Book/'};
middle = {'','10_f',...
    '10_f_r_0.6','10_f_r_0.8'};
num = [0:999];
H = 377;
r = 8;
R = 50;
c = {'','r*','k*','b*'};
for ii = 2:4
    for i = 1:numel(num)
        name = strcat(prefix,middle(ii),'/test_swi1_',num2str(num(i),'%04d'),'.h5');
        nx = double(h5read(char(name),'/Nx'));
        ny = double(h5read(char(name),'/Ny'));
        yflag1 = ny - 0.5*H;
        yflag2 = ny - H;
        pos = h5read(char(name),'/Pposition');
        np = numel(pos)/6;
        px = pos(1:3:3*np-2);
        py = pos(2:3:3*np-1);
        k1 = find(py<=yflag1);
        k2 = find(py>=yflag2);
        kk{i,ii} = intersect(k1,k2);
        Ap = numel(kk{i,ii})*pi*r*r;%particle area
        TA = (yflag1-yflag2)*nx;
        e(i,ii) = (TA-Ap)/TA;
        h(ii) = plot(i,e(i,ii),char(c(ii)));
        hold on
        drawnow
    end
end
xlabel('\itt')
ylabel('\itn')
legend(h(2:4),'\itr\rm = 0.3','\itr\rm = 0.6','\itr\rm = 0.8','location','southwest')
legend boxoff