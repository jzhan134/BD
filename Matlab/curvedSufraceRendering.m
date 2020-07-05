clear
clc
addpath('./Matlab/');
% set(0,'DefaultFigureWindowStyle','docked')
fig = figure(1);
set(fig, 'Position',  [2100, 100, 500, 500])
sb = subplot(1,1,1);
set(sb, 'position',[0 0 1 1])

[X,Y] = meshgrid(-40:2:40);
R = 40;
rim = R;
Z =  sqrt(max(R^2 - X.^2 - Y.^2,0))-sqrt(R^2-rim^2)-1;

% data = load('./traj/xyz2.dat');
data = load('./191209 coordination number study/900_R40/xyz2.dat');

pnum = max(data(:,1)+1);
frame = size(data,1)/pnum;

dispType = 2; %1: 3D, 2: top, 3: side

saveMovieName = 'changing_coord_num';
v = VideoWriter([saveMovieName,'.avi']);
v.FrameRate = 5;
open(v);

for f = 40:80
    nei = zeros(1,pnum);
    CF = data((f-1)*pnum+1:f*pnum,:);
    clf
    for i = 1:pnum
        for j = i+1:pnum
            if sqrt((CF(i,2)-CF(j,2))^2 ...
                    +(CF(i,3)-CF(j,3))^2 ...
                +(CF(i,4)-CF(j,4))^2) < (1+sqrt(3))
                nei(i) = nei(i) + 1;
                nei(j) = nei(j) + 1;
            end
        end
    end
    hold on
    plot3(CF(:,2),CF(:,3),CF(:,4), 'o','markersize', 6,...
        'MarkerEdgeColor',[.2 .2 .2],'MarkerFaceColor',[.5 .5 .5]);
    for i = 1:pnum
        if CF(i,4) ~= 0
            if nei(i) > 6
                plot3(CF(i,2),CF(i,3),CF(i,4), 'o','markersize', 6,...
                    'MarkerEdgeColor','none','MarkerFaceColor','r');
            elseif nei(i) < 6
                plot3(CF(i,2),CF(i,3),CF(i,4), 'o','markersize', 6,...
                    'MarkerEdgeColor','none','MarkerFaceColor','b');
            end
        end

    end
    if dispType ~= 2 
        surf(X,Y,Z,'EdgeColor',[.5 .5 .5],'FaceColor',[1 1 1])
    end

    if dispType == 1 % 3D view
        view(45,45)
    elseif dispType == 2 % top view
        view(0, 90)
    else % side view
        view(90,0)
    end
    axis([-50 50 -50 50 0 50])
    xticks(-50:25:50)
    yticks(-50:25:50)
%     title(num2str(f/10))
    pbaspect([1,1,.5])
    box on
    set(gca,'fontsize',12,'FontWeight','bold')
    hold off
    drawnow;
    this_frame = getframe(gcf);
    writeVideo(v,this_frame.cdata);
end
close(v)