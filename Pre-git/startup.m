set(0,'defaultFigureColor','w');
set(0,'defaultAxesPosition',[0.1621 0.1833 0.7750 0.7470]);
set(0,'defaultTextFontSize',16);
set(0,'defaultTextColor','k');
fontname = 'times new roman';
set(0,'DefaultFigurePaperPosition',[0.25 2.5 6.5 5])
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultLineLineWidth',2.)
set(0,'DefaultAxesFontName',fontname)
set(0,'DefaultFigurePosition',[10 277 560 420])
set(0,'defaulttextfontname',fontname)
set(0,'defaultaxesposition',[0.162143 0.143333 0.775 0.787])

h=figure(x)
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% h.Children(1) is the legend
h.Children(1).Location='northwest';
h.Children(1).FontSize=22;
% h.Children(2) is the axes
h.Children(2).LineWidth=2;
h.Children(2).FontSize=18;
set(findall(gca, 'Type', 'Line'),'LineWidth',2);

