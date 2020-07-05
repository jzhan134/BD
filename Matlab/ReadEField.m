function [xx,yy,E] = ReadEField(file)
x = -100:0.25:100;
[xx,yy] = meshgrid(x,x);
kT = 1.38e-23 * 298;
epsilon = 80 * 8.85e-12;
a = 1.5e-6;
fcm = 0.4667;
pref = 2*pi*epsilon*a^3*fcm/kT;
data = load(file);
Ex = reshape(data(:,3),[length(x),length(x)]);
Ey = reshape(data(:,4),[length(x),length(x)]);
E = pref*(Ex.^2 + Ey.^2);
end
