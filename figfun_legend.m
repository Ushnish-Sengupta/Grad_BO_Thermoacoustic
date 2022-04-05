function [] = figfun_legend(colmap,labels,figname)
% figfun_legend
%
% Create a figure legend with 6 lines (6 is hardwired)
% 
% INPUTS 
% colmap    matrix containing the colormap
% labels    cell array containing the six labels 

FontSize = 15;

figure(2); clf
subplot(5,1,1); hold on
axis off
plot([0 0.8],[3 3],'-','Color',colmap(1,:),'LineWidth',6)
text(0.9,3,labels{1},'FontName','Courier New','FontSize',FontSize)
plot([2 2.8],[3 3],'-','Color',colmap(2,:),'LineWidth',5)
text(2.9,3,labels{2},'FontName','Courier New','FontSize',FontSize)
plot([0 0.8],[2 2],'-','Color',colmap(3,:),'LineWidth',4)
text(0.9,2,labels{3},'FontName','Courier New','FontSize',FontSize)
plot([2 2.8],[2 2],'-','Color',colmap(4,:),'LineWidth',3)
text(2.9,2,labels{4},'FontName','Courier New','FontSize',FontSize)
plot([0 0.8],[1 1],'-','Color',colmap(5,:),'LineWidth',2)
text(0.9,1,labels{5},'FontName','Courier New','FontSize',FontSize)
plot([2 2.8],[1 1],'-','Color',colmap(6,:),'LineWidth',1)
text(2.9,1,labels{6},'FontName','Courier New','FontSize',FontSize)
plot([0 0.8],[0 0],'-','Color',colmap(7,:),'LineWidth',2)
text(0.9,0,labels{7},'FontName','Courier New','FontSize',FontSize)
plot([2 2.8],[0 0],'-','Color',colmap(8,:),'LineWidth',1)
text(2.9,0,labels{8},'FontName','Courier New','FontSize',FontSize)
axis([0 4 0 3])
eval(['print(''-depsc2'',''figures/',figname,'_legend.eps'')'])


end
