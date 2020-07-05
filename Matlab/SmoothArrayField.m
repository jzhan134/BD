clear
clc
set(0,'DefaultFigureWindowStyle','normal')
x = -100:0.25:100;
[xx,yy] = meshgrid(x,x);
file = '../library/fields/array/3b3_5v3_I.txt';

%% load and format map data
data = load(file);
Fx = reshape(data(:,5),[length(x),length(x)]);
Fy = reshape(data(:,6),[length(x),length(x)]);

%% smooth using gaussian filter
Fx_s = 2*imgaussfilt(Fx, 5);
Fy_s = 2*imgaussfilt(Fy, 5);

%% compare before and after smooth
figure(10)
colormap(jet)
clf
contourf(xx,yy,sqrt(Fx.^2 + Fy.^2),'LineColor','none');
colorbar
pbaspect([1 1 1])
axis([-100 100 -100 100])
title('before')
get_caxis = caxis;

figure(11)
clf
colormap(jet)
contourf(xx,yy,sqrt(Fx_s.^2 + Fy_s.^2),'LineColor','none');colorbar
pbaspect([1 1 1])
axis([-100 100 -100 100])
title('after')
caxis(get_caxis)

%% reformat map data for storage
Fx_new = reshape(Fx_s,[length(x)*length(x),1]);
Fy_new = reshape(Fy_s,[length(x)*length(x),1]);
data_new = data;
data_new(:,[5,6]) = [Fx_new,Fy_new];
% save('../library/fields/array/3b3_5v3_I_corr.txt','data_new','-ascii')