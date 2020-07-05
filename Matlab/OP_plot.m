clear
clc
%% constant parameters
file_path = '../traj/';
file_name = '0xyz0.dat';
opfile_name = '0op0.dat';
movie_name = 'OP_traj.avi';
use_op_data = true;
movie_maker = true;
skip_time = 1;


%% constant parameters
set(0,'DefaultFigureWindowStyle','docked')
addpath('./Image_Analysis_v0502');
figure(3)
pnum = 300;
a = 1.435;
data = load([file_path,file_name]);
frames = size(data,1)/pnum;
work_frame = 1:skip_time:frames;
% plot using op data or using traj data only

if use_op_data
    dataop = load([file_path,opfile_name]);
    Psi6 = dataop(:,2);
    C6 = dataop(:,3);
    GB = dataop(:,5);
end

%% movie maker option
if movie_maker
    v = VideoWriter(movie_name);
    v.FrameRate = 10;
    open(v);
end

%% image analysis frame by frame
for f = work_frame
    %% image analysis
    CF = data((f-1)*pnum+1:f*pnum,:);
    if ~use_op_data
        [P_Type, GB(f),GB_c, ~, Psi6(f), C6(f), L_PSI6] = IMAGE_ANALYSIS(CF);
        C6(f) = C6(f)/5.6;
    end
    
    clf;
    hold on
    %% grain boundary
    yyaxis left
    plot(1:skip_time:f,mod(GB(1:skip_time:f)-22.5,180),'r-','linewidth',4)
    plot(f,mod(GB(f)-22.5,180),'ro','markerfacecolor','r','markersize',10)
    axis([0 frames 0 180])
    xlabel('time (s)')
    ylabel('Grain Boundary Orientation')
    set(gca,'YColor','k')
    
    %% Order parameters
    yyaxis right
    region = floor(frames/50);
    for tt = 1:region
        plot([tt*50, tt*50],[0 1],'k--','linewidth',2)
    end
    % psi6
    plot(1:skip_time:f,Psi6(1:skip_time:f),'g-','linewidth',4)
    plot(f,Psi6(f),'go','markerfacecolor','g','markersize',10)
    % C6
    plot(1:skip_time:f,C6(1:skip_time:f),'b-','linewidth',4)
    plot(f,C6(f),'bo','markerfacecolor','b','markersize',10)
    xlabel('time (s)')
    ylabel('\psi_6(green) & C_6(blue)')
    axis([0 frames-1 0 1])
    pbaspect([1,1,1])
    set(gca,'fontsize',12)
    set(gca,'YColor','k')
    drawnow;
    
    %% save movie
    if movie_maker
        this_frame = getframe(gcf);
        writeVideo(v,this_frame.cdata);
    end
end
if movie_maker
    close(v)
end