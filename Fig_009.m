function [] = Fig_009()
% Fig_009
%
% Diagram of the Rijke tube

figure(1); clf; subplot(2,1,1); cla; hold on

% radius of tube
r = 0.05;

% Define LineWidths
thick = 2;
thin = 1;

% Define colours
col_vol = [0.8 0.8 0.8];

FontSize = 12;

% Volume under consideration
han_vol = patch([0.25 0.30 0.30 0.25],[-r -r +r +r],'w');
set(han_vol,'FaceColor',col_vol,'LineWidth',thin)

% X
plot([0 0.25],[-2*r,-2*r],'k','LineWidth',thin)
plot([0 0],[-1.8*r -2.2*r],'k','LineWidth',thin)
plot([0.25 0.25],[-1.8*r -2.2*r],'k','LineWidth',thin)
text(0.125,-1.6*r,'$X$','Interpreter','Latex','FontSize',FontSize,'HorizontalAlignment','center')

% delta X
plot([0.25 0.30],[-2*r,-2*r],'k','LineWidth',thin)
plot([0.25 0.25],[-1.8*r -2.2*r],'k','LineWidth',thin)
plot([0.30 0.30],[-1.8*r -2.2*r],'k','LineWidth',thin)
text(0.275,-1.6*r,'$\delta X$','Interpreter','Latex','FontSize',FontSize,'HorizontalAlignment','center')

% Perimeter
theta = -pi/4:pi/100:2*pi-pi/4;
x = 0.25*r*cos(theta);
y = r*sin(theta);
plot(-0.00+x,y,'k','LineWidth',thin)
plot(+0.01+x(1:50),y(1:50),'k','LineWidth',thin)
text(+0.03,0,'$\Gamma_c$','Interpreter','Latex','FontSize',FontSize,'HorizontalAlignment','left')

% cross-sectional area
patch(+1.00+x,y,[0.8 0.8 0.8])
plot([1.00,1.03],[0,0],'k','LineWidth',thin)
text(1.05,0,'$S_c$','Interpreter','Latex','FontSize',FontSize,'HorizontalAlignment','center')

% Rijke tube
plot([0 1],[-r -r],'k','LineWidth',thick)
plot([0 1],[+r +r],'k','LineWidth',thick)


axis off
axis equal

print('-depsc2','figures/Fig_009.eps')

