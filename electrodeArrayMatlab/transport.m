clear
clc
set(0,'DefaultFigureWindowStyle','normal')
figure(10)
set(gcf,'position',[10 10 900 900])
movie = false;
data = load('../library/configurations/3600_quench.txt');
% data = load('../traj/xyz0.dat');
pnum = max(data(:,1))+1;
a = 0.8;
frame = size(data,1)/pnum;
if movie
    v = VideoWriter('traj_demo.avi');
    v.FrameRate = 50;
    open(v);
end

for f = 1:frame
    CF = data((f-1)*pnum+1:f*pnum,:);
    CF_plot = CF*a;
    clf
    hold on
    bd = min(max(f/15,0),20);
    CF_plot(abs(abs(CF_plot(:,2)) - 30) < bd | abs(abs(CF_plot(:,3)) - 30) < bd,:) = [];
    CF_plot(abs(abs(CF_plot(:,2)) - 90) < bd | abs(abs(CF_plot(:,3)) - 90) < bd,:) = [];
    for dx = -90:60:90
        plot([dx dx],[-90 90],'k-')
        plot([-90 90],[dx dx],'k-')
    end
    for dx = -92.5:10:92.5
        for dy = -92.5:10:92.5
            rectangle('position',[dx,dy,5,5],...
                'FaceColor',[1 1 0 0.3],...
                'EdgeColor',[1 1 0 0.3])
        end
    end
%     num
% for i = 1:pnum
%     rectangle('position',...
%         [CF(i,2)*a-a, -CF(i,3)*a-a,2*a,2*a],...
%         'FaceColor','k',...
%         'EdgeColor','k',...
%         'curvature',[1 1])
% end
    plot(CF_plot(:,2),-CF_plot(:,3),'ko',...
        'markersize', 6,...
        'MarkerFaceColor','k');
    idx = round(CF(:,2)*a/60)+ round(-CF(:,3)*a/60)*3;
    for i = 1:3
        for j = 1:3
            annotation('textbox', ...
                [(i-1)/3, (j)/3-0.1 0.1 0.1], ...
            'String', length(find(idx == (j-2)*3+(i-2))), ...
            'EdgeColor', 'k', ...
            'BackgroundColor','w',...
            'HorizontalAlignment', 'left', ...
            'FontSize', 16);
        end
    end

    xticks('')
    yticks('')
    axis([-90 90 -90 90])
    pbaspect([1,1,1])
    sb = subplot(1, 1, 1);
    set(sb,'position',[0 0 1 1])
    box on
    hold off
    annotation('textbox', ...
        [.9 .0 .1 .1], ...
    'String', num2str(f/10), ...
    'EdgeColor', 'none', ...
    'BackgroundColor','w',...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom',...
    'FontSize', 16);
    drawnow;
    if movie
        this_frame = getframe(gcf);
        writeVideo(v,this_frame.cdata);
    end
%     pause(0.2)
end
if movie
    close(v)
end