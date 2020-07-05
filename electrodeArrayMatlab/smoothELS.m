file ="../library/fields/array/3x3_quad.txt";
set(0,'DefaultFigureWindowStyle','normal')
x = -100:0.25:100;
[xx,yy] = meshgrid(x,x);
kT = 1.38e-23 * 298;
epsilon = 80 * 8.85e-12;
a = 1.5e-6;
fcm = 0.4667;
pref = 2*pi*epsilon*a^3*fcm/kT;

%% read electric field file
data = load(file);
Ex = reshape(data(:,3),[length(x),length(x)]);
Ey = reshape(data(:,4),[length(x),length(x)]);
Fx = reshape(data(:,5),[length(x),length(x)]);
Fx = reshape(data(:,5),[length(x),length(x)]);
E = sqrt(Ex.^2 + Ey.^2)./(1e4);
Ex_filt = 1*imgaussfilt(Ex,5);
Ey_filt = 1*imgaussfilt(Ey,5);
E_filt = sqrt(Ex_filt.^2 + Ey_filt.^2)./(1e4);
for i = 1:801
    for j = 1:801
        if (E(i,j) > 4)
            E(i,j) = 4;
        end
        if (E_filt(i,j) > 4)
            E_filt(i,j) = 4;
        end
    end
end
            
%% energy landscape in the unit of kT
figure(11)
clf
set(gcf,'position',[80 1250 400 400])
sb = subplot(1,1,1);
set(sb,'position',[0 0 1 1])
colormap((jet))
hold on 
contourf(xx,yy, E',50,'linecolor','none');
pbaspect([1 1 1])
axis([-90 90 -90 90])
xticks('')
yticks('')
caxis([0 5])
box on

figure(12)
clf
set(gcf,'position',[500 1250 400 400])
sb = subplot(1,1,1);
set(sb,'position',[0 0 1 1])
colormap((jet))
hold on 
contourf(xx,yy, E_filt',50,'linecolor','none');
pbaspect([1 1 1])
axis([-90 90 -90 90])
xticks('')
yticks('')
caxis([0 5])
box on
data_filt = data;
data(:,3) =reshape(Ex_filt,[801*801,1]);
data(:,4) =reshape(Ey_filt,[801*801,1]);
save(file,'data','-ascii');