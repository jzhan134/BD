clear
clc
addpath('./Matlab/');
set(0,'DefaultFigureWindowStyle','docked')
figure(10)
movie = false;
data = load('./traj/xyz1.dat');
pnum = 300;
frame = size(data,1)/pnum;
x = [];
y = [];

for f = 50:56
    x = cat(2,x,data((f-1)*pnum+1:f*pnum,2));
    y = cat(2,y,data((f-1)*pnum+1:f*pnum,3));
end
x = [x, nan(300,1)];
y = [y, nan(300,1)];

clf
hold on
colormap(flipud(jet));
c = 1:8;
for i = 1:pnum
    patch(y(i,:), x(i,:), c, 'EdgeColor','interp','linewidth',3)
end
axis([-60 60 -60 60])
pbaspect([1 1 1])
box on
xticks(-60:30:60)
yticks(-60:30:60)
set(gca,'fontsize',30,'FontWeight','bold')