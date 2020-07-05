%{
    L_PSI6 ranges between [0,1]
    L_theta6 ranges between [-30,30], -30 when no surrounding particle is available
    neighboridx: list of surrounding particles including the center particle itself
    G_PSI6: global psi6 value between [0,1]
    G_C6: global c6 value between [0,1]
%}

function [L_PSI6, L_C6, neighboridx, G_PSI6, G_C6,L_theta6] = ORDER_PARAMETER(CurrFrame)

pnum = size(CurrFrame,1);   % particle number
psir = zeros(pnum,1);       % real part of local psi6
psii = zeros(pnum,1);       % imaginary part of local psi6

% output variables
neighboridx = cell(pnum,1); % list of neighbors indices for each particle
L_theta6 = zeros(1,pnum);   % theta value of each particle
PSI6 = [0 0];               % real and imaginary parts of global psi6
L_PSI6 = zeros(1,pnum);     % psi6 value of each particle
L_C6 = zeros(1,pnum);
dist= zeros(pnum,pnum);
dh = 2.64;


% calculate pair-wise particle distances
for p1 = 1:pnum
    for p2 = 1:pnum
        dx(p1,p2) = CurrFrame(p1,2)-CurrFrame(p2,2);
        dy(p1,p2) = CurrFrame(p1,3)-CurrFrame(p2,3);
        dist(p1,p2) = sqrt(dx(p1,p2).^2 + dy(p1,p2).^2);
    end
end

% sum up local op values
for p1 = 1:pnum
    for p2 = p1+1:pnum
        if dist(p1,p2) < dh
            
            % particle-particle angles
            theta = atan(dy(p1,p2)/dx(p1,p2));
            L_theta6(p1) = L_theta6(p1) + theta;
            L_theta6(p2) = L_theta6(p2) + theta;
            
            % real and imaginary parts of local psi6
            psir(p1) = psir(p1) + cos(6*theta);
            psii(p1) = psii(p1) + sin(6*theta);
            psir(p2) = psir(p2) + cos(6*theta);
            psii(p2) = psii(p2) + sin(6*theta);  
            
            % neighbor list
            neighboridx{p1} = cat(2,neighboridx{p1},p2);
            neighboridx{p2} = cat(2,neighboridx{p2},p1);
        end
    end
end

for p1 = 1:pnum
    
    if size(neighboridx{p1},2) ~= 0
        L_C6(p1) = L_C6(p1)/size(neighboridx{p1},2);
        L_theta6(p1) = L_theta6(p1)/size(neighboridx{p1},2)*180/pi; 
        psii(p1) = psii(p1)/size(neighboridx{p1},2);
        psir(p1) = psir(p1)/size(neighboridx{p1},2);
        L_PSI6(p1) = sqrt(psii(p1)^2 + psir(p1)^2);
        
        PSI6(1) = PSI6(1) + psii(p1)/pnum;
        PSI6(2) = PSI6(2) + psir(p1)/pnum;
    else
        L_theta6(p1) = -30;
    end
    
    neighboridx{p1} = sort(cat(2,neighboridx{p1},p1)); 
    
end


for p1 = 1:pnum
    for p2 = p1+1:pnum
        if dist(p1,p2) < dh
            numer = psii(p1)*psii(p2) + psir(p1)*psir(p2);
            denom = sqrt(numer^2 + (psir(p1)*psii(p2) - psir(p2)*psii(p1))^2);
            if numer/denom > 0.32
                L_C6(p1) = L_C6(p1) + 1/5.6;
                L_C6(p2) = L_C6(p2) + 1/5.6;
            end
        end
    end
end

for p = 1:pnum
    if L_C6(p) >= 1
        L_C6(p) = 1;
    end
end
G_PSI6 = sqrt(PSI6(1)^2 + PSI6(2)^2);
G_C6 = mean(L_C6);
end