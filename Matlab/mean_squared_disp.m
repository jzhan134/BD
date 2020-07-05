clear
clc
pnum = 300;
a = 1.435;
MSD = [];
%% Accumulate MSD(s) for same depletion and case
for i = 0:3
    for j = 0:2
        data = load(['../sample/',num2str(i),'xyz',num2str(j),'.dat']);
        data(1:100*300,:) = [];
        data(:,[2,3]) = data(:,[2,3])*a;
        for t = 1:99
            for t_start = 1:100-t
                data_start = data((t_start-1)*pnum+1 : t_start*pnum,[2,3]);
                data_end = data((t_start+t-1)*pnum+1 : (t_start+t)*pnum,[2,3]);
                MSD_temp(t) = mean((data_end(:,1) - data_start(:,1)).^2 + (data_end(:,2) - data_start(:,2)).^2) ;
            end
        end
        MSD = cat(1,MSD, MSD_temp);
    end
end
MSD = mean(MSD);
figure(11)
plot(1:99,MSD)
figure(21)
D = 0.25*1000 * diff(MSD);
plot(1:98,D)
% %% Calculate average MSD(s) by number of data at each second and save
% for i = 2:101
%     for type = 1:3
%         if Count{1,type}(1,i) ~= 0
%             MSD{1,type}(1,i) = MSD{1,type}(i)/Count{1,type}(i)*1.435^2;
%         end
%     end
% end
% figure(1)
% hold on
% Color = {'r','g','b'};
% for type = 1:3
%     plot(1:99,MSD{1,type},Color{type},'linewidth',2);
% end
% legend('Grain','Inner','Outer')
% xlabel('Time (ms)')
% ylabel('MSD (\mum^2/s)')
% save('./MSD_0kT_Short.mat','MSD');
% 
% %% DSL
% figure(2)
% hold on
% for type = 1:3
%     DSL(type,:) = 0.25*1000*(MSD{1,type}(2:end) - MSD{1,type}(1:end-1));
%     plot(1:100,DSL(type,:),Color{type},'linewidth',2);
% end
% % axis([ 0 100 0 0.002])
% xlabel('Time (s)')
% ylabel('D_S^L (\mum^2/s)')
% for type = 1:3
%     DSL_Val(type) = DSL(type,1);
% %     DSL_Val(type) = mean(DSL(type,20:50));
% end
% save('3kT-DSS.dat','DSL_Val','-ascii')