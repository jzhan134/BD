
set(0,'DefaultFigureWindowStyle','docked')
x = -100:0.25:100;
[xx,yy] = meshgrid(x,x);
kT = 1.38e-23 * 298;
epsilon = 80 * 8.85e-12;
a = 1.5e-6;
fcm = 0.4667;
pref = 2*pi*epsilon*a^3*fcm/kT;
c = [];
for tt = 2:6
file =['./library/fields/Assembly/elsForAnisotropy/10v',num2str(tt),'.txt'];
data = load(file);
Ex = reshape(data(:,3),[length(x),length(x)]);
Ey = reshape(data(:,4),[length(x),length(x)]);
E = pref*(Ex.^2 + Ey.^2);
for i = 1:801
    for j = 1:801
        if (sqrt((i-401)^2 + (j-401)^2) > 60*4 || E(i,j) > 100)
            E(i,j) = 0;
        end
    end
end
maxGyr = (sum(sum(E.*xx.*xx)))/sum(sum(E));
minGyr = (sum(sum(E.*yy.*yy)))/sum(sum(E));
c = cat(2,c,(sqrt(maxGyr)-sqrt(minGyr))/sqrt(maxGyr + minGyr));
figure(tt)
colormap((jet))
contourf(xx,yy, E','linecolor','none');
pbaspect([1 1 1])
axis([-90 90 -90 90])
xticks(-90:30:90)
yticks(-90:30:90)
axis on
box on
xlabel('x/ um')
ylabel('y/ um')
set(gca,'fontsize',14)
% drawnow
end
