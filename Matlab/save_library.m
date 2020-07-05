idx = 0;
for i = 0:99
    file = ['../traj/xyz',num2str(i),'.dat'];
    if exist(file)
        data = load(file);
        CF = data(end-599:end,:);
        save(['../library/configurations/600/sphere',num2str(idx),'.txt'],'CF','-ascii');
        idx = idx + 1;
    end
end