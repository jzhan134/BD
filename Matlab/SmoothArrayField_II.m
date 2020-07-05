file = '../library/fields/array/assembly/quad.txt';
fig = figure(11);
set(fig, 'Position',  [1100, 100, 700, 700])
sb = subplot(1,1,1);
set(sb, 'position',[0 0 1 1])
x = -100:0.25:100;
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
Fy = reshape(data(:,6),[length(x),length(x)]);
mid_x = Fx(401-120:401+120,401-120:401+120);
mid_y = Fy(401-120:401+120,401-120:401+120);

blankx = zeros(801,801);
blanky = zeros(801,801);
blankEx = zeros(801,801);
blankEy = zeros(801,801);

% (1,1)
blankx(401-360:401-120,401-360:401-120) = mid_x;
blanky(401-360:401-120,401-360:401-120) = mid_y;
% (1,2)
blankx(401-360:401-120,401-120:401+120) = mid_x;
blanky(401-360:401-120,401-120:401+120) = mid_y;
% (1,3)
blankx(401-360:401-120,401+120:401+360) = mid_x;
blanky(401-360:401-120,401+120:401+360) = mid_y;
% (2,1)
blankx(401-120:401+120,401-360:401-120) = mid_x;
blanky(401-120:401+120,401-360:401-120) = mid_y;
% (2,2)
blankx(401-120:401+120,401-120:401+120) = mid_x;
blanky(401-120:401+120,401-120:401+120) = mid_y;
% (2,3)
blankx(401-120:401+120,401+120:401+360) = mid_x;
blanky(401-120:401+120,401+120:401+360) = mid_y;
% (3,1)
blankx(401+120:401+360,401-360:401-120) = mid_x;
blanky(401+120:401+360,401-360:401-120) = mid_y;
% (3,2)
blankx(401+120:401+360,401-120:401+120) = mid_x;
blanky(401+120:401+360,401-120:401+120) = mid_y;
% (3,3)
blankx(401+120:401+360,401+120:401+360) = mid_x;
blanky(401+120:401+360,401+120:401+360) = mid_y;

mid_x = Ex(401-120:401+120,401-120:401+120);
mid_y = Ey(401-120:401+120,401-120:401+120);
% (1,1)
blankEx(401-360:401-120,401-360:401-120) = mid_x;
blankEy(401-360:401-120,401-360:401-120) = mid_y;
% (1,2)
blankEx(401-360:401-120,401-120:401+120) = mid_x;
blankEy(401-360:401-120,401-120:401+120) = mid_y;
% (1,3)
blankEx(401-360:401-120,401+120:401+360) = mid_x;
blankEy(401-360:401-120,401+120:401+360) = mid_y;
% (2,1)
blankEx(401-120:401+120,401-360:401-120) = mid_x;
blankEy(401-120:401+120,401-360:401-120) = mid_y;
% (2,2)
blankEx(401-120:401+120,401-120:401+120) = mid_x;
blankEy(401-120:401+120,401-120:401+120) = mid_y;
% (2,3)
blankEx(401-120:401+120,401+120:401+360) = mid_x;
blankEy(401-120:401+120,401+120:401+360) = mid_y;
% (3,1)
blankEx(401+120:401+360,401-360:401-120) = mid_x;
blankEy(401+120:401+360,401-360:401-120) = mid_y;
% (3,2)
blankEx(401+120:401+360,401-120:401+120) = mid_x;
blankEy(401+120:401+360,401-120:401+120) = mid_y;
% (3,3)
blankEx(401+120:401+360,401+120:401+360) = mid_x;
blankEy(401+120:401+360,401+120:401+360) = mid_y;

clf
[xx,yy] = meshgrid(-100:.25:100);
E = pref * (blankEx.^2 + blankEy.^2);
contour(xx,yy, E, 0:2:50);
colorbar;
pbaspect([1 1 1])
% axis([-40 40 -40 40])
axis([-100 100 -100 100])
axis on
box on
xlabel('x/a')
ylabel('y/a')
set(gca,'fontsize',14)

Fx_new = reshape(blankx,[length(xx)*length(xx),1]);
Fy_new = reshape(blanky,[length(xx)*length(xx),1]);
Ex_new = reshape(blankEx,[length(xx)*length(xx),1]);
Ey_new = reshape(blankEy,[length(xx)*length(xx),1]);
data_new = data;
data_new(:,[3,4,5,6]) = [Ex_new, Ey_new, Fx_new,Fy_new];
save('../library/fields/array/assembly/quad.txt','data_new','-ascii')