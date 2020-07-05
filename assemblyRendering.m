clear
clc
addpath('./Matlab/Image_Analysis_v0626');
set(0,'DefaultFigureWindowStyle','normal')
fig = figure(11);
set(fig, 'Position',  [100, 1100, 700, 700])
sb = subplot(1,1,1);
set(sb, 'position',[0 0 1 1])
movie = 0;
colorRendering = 1; % 0:grayscale; 1:full color
% data = load('./library/configurations/library300/sphere10.txt');
data = load('./traj/300/xyz0.dat');
if movie
    v = VideoWriter('10v10.avi');
    v.FrameRate = 10;
    open(v);
end
pnum = max(data(:,1)) + 1;
a = 1.435;
frame = size(data,1)/pnum;

for f = frame
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% image analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%
    CF = data((f-1)*pnum+1:f*pnum,:);
    CF_plot = [(1:pnum)',CF(:,[2,3])*a];
    clf
    hold on
    gyr = [mean((CF_plot(:,2) - mean(CF_plot(:,2))).^2),...
        mean((CF_plot(:,2) - mean(CF_plot(:,2))).*(CF_plot(:,3) - mean(CF_plot(:,3))));...
        mean((CF_plot(:,2) - mean(CF_plot(:,2))).*(CF_plot(:,3) - mean(CF_plot(:,3)))),...
        mean((CF_plot(:,3) - mean(CF_plot(:,3))).^2)];
    S = eig(gyr);
    Rg = sqrt(S(1)) + sqrt(S(2));
    tt = abs(sqrt(S(1)) - sqrt(S(2)));
    c = tt/Rg;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot config. %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~colorRendering
%         contour(xx,yy, E',[0:2:4,10:10:500]);
        plot(CF_plot(:,2),CF_plot(:,3),'ko','markersize', 4, ...
            'MarkerFaceColor','k');
    else
        [P_Type, GB, domain, Psi6, G_C6, ~, L_C6,theta] = IMAGE_ANALYSIS(CF);
        sc = scatter(CF_plot(:,2),CF_plot(:,3), 50,...
            'MarkerFaceColor',[1,0,0],...
            'MarkerEdgeColor','none',...
            'MarkerFaceAlpha',Psi6);
        for pt = 1:pnum
%             if CF(pt,5) == 2  % GB
%                 plot(CF_plot(pt,2), CF_plot(pt,3),'o','markersize',...
%                     4,'MarkerFaceColor','y','MarkerEdgeColor','y');                
%             end
% 
            if P_Type(pt) == -3 && G_C6>=0.8 && Psi6 <= 0.93  % GB
                plot(CF_plot(pt,2), CF_plot(pt,3),'o','markersize',...
                    6.5,'MarkerFaceColor','y','MarkerEdgeColor','y');                
            end
        end
%         if G_C6>=0.8 && Psi6 <= 0.93
            plot((-25:25)*cos(GB(1)*pi/180)+GB(2)*a, ...
                (-25:25)*sin(GB(1)*pi/180)+GB(3)*a,...
                '-','linewidth',4,'color','k')
%         end
    end
        
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% figure setting %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     axis([-40 40 -40 40])
    sb = subplot(1,1,1);
    set(sb,'position',[0 0 1 1])
    axis([-90 90 -90 90])
    pbaspect([1,1,1])
%     title(num2str(f))
    box on
%     set(gca,'fontsize',12,'FontWeight','bold')
    hold off
    drawnow;
    if movie
        this_frame = getframe(gcf);
        writeVideo(v,this_frame.cdata);
    end
end
if movie
    close(v)
end
