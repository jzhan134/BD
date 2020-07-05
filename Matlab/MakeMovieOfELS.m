v = VideoWriter('~/Desktop/fan_els.avi');
v.FrameRate = 5;
open(v);
for rep = 0:3
for i = 0:11
    file = ['./library/fields/array/rot70_',num2str(i),'.txt'];
    EnergyLandscapeVisualization;
    caxis([0 3e4])
    this_frame = getframe(gcf);
    writeVideo(v,this_frame.cdata);
end
end
close(v);