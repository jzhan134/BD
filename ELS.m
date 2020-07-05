file ="./library/fields/Assembly/quadrupole.txt";
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
E = sqrt(Ex.^2 + Ey.^2);
U = pref*(Ex.^2 + Ey.^2);
U(U>100) = NaN;      
%% energy landscape in the unit of kT
figure(1)
clf
colormap((jet))
hold on 
contourf(xx,yy, U',50,'linecolor','none');
pbaspect([1 1 1])
axis([-50 50 -50 50])
axis on
box on
xlabel('x/ um')
ylabel('y/ um')
colorbar
set(gca,'fontsize',14)

figure(2)
clf
E0 = 1/140e-6;
hold on
plot(0:1:100,E(401:4:end,401)/E0,'x')
plot(0:50, (0:50)*4/140)
axis([0 50 0 3])