clear
clc
prefix = {'/home/user/dlvo/test_rw1_'};
num = [0:999];
R = 5;
for i=1:numel(num)
    name = strcat(prefix,num2str(num(i),'%04d'),'.h5');
    rp = h5read(char(name),'/RWPposition');
    rpad = h5read(char(name),'/RWPIsAD');
    p = h5read(char(name),'/Pposition');
    np = numel(p)/6;
    px = p(1:3:3*np-2);
    gpx = p(3*np+1:3:end-2);
    gpy = p(3*np+2:3:end-1);
    py = p(2:3:3*np-1);
    rpx = rp(1:3:end-2);
    rpy = rp(2:3:end-1);
    RR = zeros(np,1)+R;
    viscircles([px,py],RR,'Color','b');
    hold on
    viscircles([gpx,gpy],RR,'Color','b');
    kkk1 = find(rpad<0);
    kkk2 = find(rpad>0);
    ADN(i,1) = numel(kkk2);
    plot(rpx(kkk1),rpy(kkk1),'k.');
    plot(rpx(kkk2),rpy(kkk2),'r.');
    axis equal
    xlim([0 90])
    drawnow
%     fn = strcat(num2str(i),'.jpg');
%     saveas(gcf,char(fn))
    clf
end
figure
plot(ADN)