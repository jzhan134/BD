clear
clc
addpath('./Matlab/');
file = './library/fields/array/random_oct_0_0.txt';
data = load(file);
Ex = [];
Ey = [];
Fx = [];
Fy = [];
[xx, yy] = meshgrid(-100:.25:100);
while(~isempty(data))
    Ex = cat(2,Ex, data(1:801,3));
    Ey = cat(2,Ey, data(1:801,4));
    Fx = cat(2,Fx, data(1:801,5));
    Fy = cat(2,Fy, data(1:801,6));
    data(1:801,:) = [];
end

%% energy landscape in the unit of kT
figure(11)
clf
colormap(jet)
% contour(xx,yy,sqrt(Ex.^2 + Ey.^2));
contour(xx,yy,Fy);
colorbar
pbaspect([1 1 1])
axis([-100 100 -100 100])
axis on
xlabel('x/\mum')
ylabel('y/\mum')
set(gca,'fontsize',14)

% range to capture
range = [-30, 30, -30, 30];
corr_range = 401 + range .*4;
applied_range = [30 90 -30 30];
corr_app_range = 401 + applied_range .*4;

% captured window
captured_Ex = Ex(corr_range(1):corr_range(2),corr_range(3):corr_range(4));
captured_Ey = Ey(corr_range(1):corr_range(2),corr_range(3):corr_range(4));
captured_Fx = Fx(corr_range(1):corr_range(2),corr_range(3):corr_range(4));
captured_Fy = Fy(corr_range(1):corr_range(2),corr_range(3):corr_range(4));

% add to clear background
newEx = zeros(801, 801);
newEy = zeros(801, 801);
newFx = zeros(801, 801);
newFy = zeros(801, 801);
newEx(corr_app_range(1):corr_app_range(2),corr_app_range(3):corr_app_range(4)) = captured_Ex;
newEy(corr_app_range(1):corr_app_range(2),corr_app_range(3):corr_app_range(4)) = captured_Ey;
newFx(corr_app_range(1):corr_app_range(2),corr_app_range(3):corr_app_range(4)) = captured_Fx;
newFy(corr_app_range(1):corr_app_range(2),corr_app_range(3):corr_app_range(4)) = captured_Fy;

figure(12)
clf
colormap(jet)
% contour(xx,yy,sqrt(newEx.^2 + newEy.^2));
contour(xx,yy,newFx);
% contour(xx,yy,newFy);

colorbar
pbaspect([1 1 1])
axis([-100 100 -100 100])
axis on
xlabel('x/\mum')
ylabel('y/\mum')
set(gca,'fontsize',14)

Ex = [];
Ey = [];
Fx = [];
Fy = [];
while(~isempty(newEx))
    Ex = cat(1,Ex, newEx(:,1));
    Ey = cat(1,Ey, newEy(:,1));
    Fx = cat(1,Fx, newFx(:,1));
    Fy = cat(1,Fy, newFy(:,1));
    newEx(:,1) = [];
    newEy(:,1) = [];
    newFx(:,1) = [];
    newFy(:,1) = [];
end
NewTable = [zeros(801*801,1),zeros(801*801,1),...
    Ex,Ey,Fx,Fy,...
    zeros(801*801,1),zeros(801*801,1),...
    zeros(801*801,1),zeros(801*801,1)];
save('./library/fields/array/oct_at_60_0.txt','NewTable','-ascii')