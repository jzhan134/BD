clear
clc
% close all
tic
corr = load('corr.dat');
load('10v2.mat');
E = E*2*1.5^2;
[xx,yy]=meshgrid(-100:.25:100);

currCorr = getCorrectedCorrelation(corr, 0.903);

currEta = zeros(size(E));
for i = 1:size(E,1)
    for j = 1:size(E,2)
        if ~isnan(E(i,j)) && sqrt(((i-401)/4)^2 + ((j-401)/4)^2) < (55/cos(22.5*pi/180)-7)
            currEta(i,j) = interpolate_eta(currCorr,E(i,j));
        end
    end
end
currN = sum(sum(currEta*.25^2/(pi*1.5^2)))

figure(12)
clf
colormap((jet))
hold on 
contourf(xx,yy, currEta',50,'linecolor','none');
pbaspect([1 1 1])
axis([-90 90 -90 90])
xticks(-90:30:90)
yticks(-90:30:90)
axis on
box on
xlabel('x/ um')
ylabel('y/ um')
set(gca,'fontsize',14)