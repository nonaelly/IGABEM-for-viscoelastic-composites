function setPlot(xtitle, ytitle, lgd)
% Set the configuration of picture.

SizeT = 16;
% Set legend.
lgd.FontSize = SizeT;
lgd.Location = 'east';
lgd.FontWeight = 'bold';
lgd.FontName = 'Times New Roman';
lgd.Interpreter = 'latex';

xlabel(xtitle,'Interpreter','latex')
ylabel(ytitle,'Interpreter','latex')

set(gca,'LineWidth',1.5,'FontSize',SizeT,'FontWeight','bold','FontName','Times New Roman')
set(gcf,'Position',[100 100 750 600])

ax = gca;
ax.TickLabelInterpreter = 'latex';

box on
grid on
end