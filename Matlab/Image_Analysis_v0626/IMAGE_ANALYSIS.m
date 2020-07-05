%%%%%%%%%%%%%% Automatic Image Analysis - real-time version %%%%%%%%%%%%%
%{ 
data structure
    CurrFrame: [particle index, x( in unit of a), y]

    P_Type: partical category (<0 perpheral particles (unused), =-3 grain
    boundary, >0 crystalline particles differentiate by domains indice

    GB_INFO: [grain boundary orientation (-67,112), grain boundary center x,
    grain boundary center y]. GB orientation between +-90 from 22.5 field
    direction

    DOMAIN_INFO: [domain center x, domain center y, domain orientation 
    (-30 ~ +30), domain size] captured up to two domains

    G_PSI6: global psi_6
    G_C6: global C_6
    
    L_PSI6: local psi_6 of each particle
    L_C6: local C_6 of each particle
    L_theta6: local lattice orientation of each particle
%}
function [P_Type, GB_INFO, DOMAIN_INFO, G_PSI6, G_C6, L_PSI6, L_C6,L_theta6] = IMAGE_ANALYSIS(CurrFrame)
%% format the input data
if CurrFrame(1,1) == 0
    CurrFrame(:,1) = (1:size(CurrFrame,1))';
end
    

%% differentiate particle types
[P_Type, neighboridx, L_C6, G_PSI6, G_C6,L_PSI6,L_theta6] = CATEGORIZATION(CurrFrame);
    
%% identify major grain boundary
if ~isempty(CurrFrame(P_Type == 0,:))
    GB = [];
    GB_list = CONNECTIVITY2(CurrFrame(P_Type == 0,1), neighboridx);
    for i = 1:length(GB_list)
        
        if length(GB_list{i})>6
            GB = cat(2,GB,GB_list{i});
        end
    end
%     GB = GB{1};
%     P_Type(P_Type == 0) = -3;
    P_Type(GB) = -3;
end
    
    
%% group crystalline particles into domains.
% major domain = 1; minor domain = 2; other domains = 3
if max(P_Type) >= 1
    Type = CONNECTIVITY(CurrFrame(P_Type > 0,1), neighboridx);
    P_Type(Type{1}) = 1;
    P_Type(Type{1}) = 2;
    for i = 3:size(Type,2)
        P_Type(Type{i}) = 3;
    end
end


%% Format domain output information
DOMAIN_INFO = nan(4,min(max(P_Type),2));
for i = 1:min(max(P_Type),2)
    DOMAIN_INFO(:,i) = [mean(CurrFrame((P_Type==i),2)),...
        mean(CurrFrame((P_Type==i),3)),...
        median(L_theta6((P_Type==i))), ...
        length(CurrFrame((P_Type==i),2))];
end


%% Format grain boundary output information
GB_INFO = nan(1,3);
if ~isempty(CurrFrame(P_Type==-3,:))
    x_mean = mean(CurrFrame(P_Type == -3,2));
    y_mean = mean(CurrFrame(P_Type == -3,3));
    x = CurrFrame(P_Type == -3,2);
    y = CurrFrame(P_Type == -3,3);
    bin = -60:120;
    sig = NaN(1,size(bin,2));
    for i = 1:size(bin,2)
        degree = bin(i)*pi/180;
        fac1 = abs(x.*tan(degree) - y + (-tan(degree)*x_mean + y_mean));
        fac2 = sqrt(tan(degree).^2+1);
        sig(i) = sqrt(sum((fac1./fac2).^2)/size(x,1));
    end
    [~,best_Dirc] = min(sig);
    GB_INFO(1) = bin(best_Dirc);
    GB_INFO(2:3) = [x_mean,y_mean];
end
end
