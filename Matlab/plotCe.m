ce = [0.4067, 0.512, 0.6965, 0.982;...
    0.4267, 0.505, 0.7043, 0.992;
    0.4367, 0.522, 0.7243, 0.992];
cb = {'r','g','b'};
figure(1)
clf
sb = subplot(1,1,1);
set(sb,'position',[0 0 1 1])
hold on
set(gcf,'position',[80,1100, 153, 81])
for i = 1:3
    plot(.4:.2:1,ce(i,:),'o','markeredgecolor',cb{i},'markerfacecolor','none','markersize',4);
end
axis([0 1.05 0 1.05])
xticks(0:.25:1)
yticks(0:.25:1)
box on
