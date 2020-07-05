function [Particle_Category, neighboridx, L_C6, G_PSI6, G_C6, L_PSI6,L_theta6] = CATEGORIZATION(CurrFrame)

%% calculate order paramters
Particle_Category = zeros(size(CurrFrame,1),1);
CurrFrame(:,2) = CurrFrame(:,2) - mean(CurrFrame(:,2));
CurrFrame(:,3) = CurrFrame(:,3) - mean(CurrFrame(:,3));
[L_PSI6, L_C6, neighboridx, G_PSI6, G_C6,L_theta6] = ORDER_PARAMETER(CurrFrame);


%% first-pass categorization
% 1. outer rim: particles if less than 6 neighbors (including self) exist. 
% 2. GB:        particles if low local psi_6 or less than 7 neighbors
% 3. crystal:   particles otherwise
% 4. inner rim: if a particle is not an outer rim particle but has a outer
%               rim particle as neighbor
for i = 1:size(CurrFrame,1) 
    if  (size(neighboridx{i},2) < 5)
        Particle_Category(i) = -1;
    elseif L_PSI6(i) < 0.95  || L_C6(i) < 0.95 %grain boundary
%     elseif L_PSI6(i) < 0.95  && size(neighboridx{i},2) > 4 %grain boundary
        Particle_Category(i) = 0;
    else
        Particle_Category(i) = 1;
    end
end
for i = 1:size(CurrFrame,1) 
    if any(Particle_Category(neighboridx{i}) == -1) ...
            && Particle_Category(i) ~= -1
        Particle_Category(i) = -2;
    end
end


%% trim inner and outer rims
% group outer and inner particles together and divide them by their
% connectivity. If any minor divison exist, re-define them as grain
% boundary particle
%     tem_list = CurrFrame(Particle_Category==-1 | Particle_Category==-2,1);
%     if ~isempty(tem_list)
%         rim = CONNECTIVITY(tem_list, neighboridx);
%         for domain = 2:size(rim,2)
%             Particle_Category(rim{domain}) = 0;
%         end
%     end

%% trim outer rim
% there is a chance that a particle inside the ensemble has a local fluid
% state and therefore be defined as outer rim particle.
% this kind of particle should be isolated from other outer rim particles,
% so to distinguish them one can measure its nearest 3 outer rim neighbors
% for i = find(Particle_Category == -1) 
%     dist = sort((CurrFrame(Particle_Category == -1,2) - CurrFrame(i,2)).^2 + ...
%        (CurrFrame(Particle_Category == -1,3) - CurrFrame(i,3)).^2);
%    if dist(3) > 4
%        Particle_Category(i) = 0;
%    end
% end
% for i = 1:size(CurrFrame,1) 
%    if Particle_Category(i) == -1 && ...
%             size(find(Particle_Category(neighboridx{i}) == -2),1) >= 3 &&...
%             size(find(Particle_Category(neighboridx{i}) == -1),1) <= 1
%         Particle_Category(i) = 0;
%         
%    end
% end

%% trim inner rim based on updated rim particles
% Particle_Category(Particle_Category == -2) = 0;
for i = 1:size(CurrFrame,1) 
    if Particle_Category(i) ~= -1
        if any(Particle_Category(neighboridx{i}) == -1)
            Particle_Category(i) = -2;
        end
    end
end

end