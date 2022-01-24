function [fig12]=figure_cpt(y,cpt_fit,CPT)

fig12=figure
x0=10; y0=10; width=900; height=800;
set(gcf,'position',[x0,y0,width,height])
ax=subplot(1,1,1);
plot(y,cpt_fit(:,34),'Color',[0 0 0],'LineWidth',2); hold on;
plot(y,cpt_fit(:,33),'Color',[0 0 1],'LineWidth',2)
plot(y,cpt_fit(:,28),'Color',[0 0.5 0],'LineWidth',2)
plot(y,cpt_fit(:,7),'Color',[1 0 0],'LineWidth',2)
plot(y,CPT(:,34),'^','Color',[0 0 0],'MarkerSize',7,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0,0,0]);
plot(y,CPT(:,33),'s','Color',[0 0 1],'MarkerSize',7,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0,0,1]);
plot(y,CPT(:,28),'d','Color',[0 0.5 0],'MarkerSize',7,'MarkerEdgeColor',[0 0.5 0],'MarkerFaceColor',[0,0.5,0]);
plot(y,CPT(:,7),'o','Color',[1 0 0],'MarkerSize',7,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1,0,0]);
ylim([0 0.04]);

set(ax,'XGrid','on'); set(ax,'YGrid','on'); 
set(ax,'FontName','Times New Roman','Fontsize',22);
set(ax,'Color',[1 1 1],'Box','on'); set(ax,'YTick',[0 0.01 0.02 0.03 0.04],'YTicklabel',{0 0.01 0.02 0.03 0.04});
xlabel('y/S','FontName','Times New Roman','FontSize',28,'FontAngle','Italic','verticalalignment','middle')
ylabel('C_p_t(y)','FontName','Times New Roman','FontSize',28,'FontAngle','Italic','rotation',0,'Position',[-0.0925,0.02175,-1])
legend('(a): Re=120 f^+=0.64 φ=0.5', '(b): Re=120 f^+=0.61 φ=0.7', '(c): Re=120 f^+=1.29 φ=0.5', '(d): Re=75 f^+=0.64 φ=0.5','location','northwest','Fontsize',19)
