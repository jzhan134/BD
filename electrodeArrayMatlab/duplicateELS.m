file ="../library/fields/array/3x3_5v5_raw.txt";
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
data_corr = data;
for i = 3:6
    take = reshape(data(:,i),[length(x),length(x)]);
    take_corr = zeros(801,801);
    for dx = -1:1
        for dy = -1:1
            take_corr((401-120:401+120)+dx*240,(401-120:401+120)+dy*240) = take(401-120:401+120,401-120:401+120);
        end
    end
    data_corr(:,i) = reshape(take_corr,[length(x)*length(x), 1]);
end
save("../library/fields/array/3x3_5v5_corr.txt",'data_corr','-ascii');