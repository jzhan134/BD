Fx = [];
Fy = [];
while(~isempty(Fy2))
    Fy = cat(1,Fy, Fy2(:,1));
    Fx = cat(1,Fx, Fx2(:,1));
    Fx2(:,1) = [];
    Fy2(:,1) = [];
end
file = './library/fields/array/rot70_11.txt';
data = load(file);
data(:,5) = Fx;
data(:,6) = Fy;
save('./library/fields/array/rot70_11_smooth.txt','data','-ascii')