pathname=dir('*.txt');
pf=0;EoS=[];figure();hold on;grid on;
for pf=1:length(pathname)
EoS=importdata(pathname(pf).name);
plot(EoS(:,1),EoS(:,2),'.');
end
xlabel('\bf log Energy (kg/m^3)','FontSize',18);ylabel('\bf log Pressure (Pa)','FontSize',18);
set(gca,'Box','on','xminortick','on','yminortick','on','TickDir','in','TickLength',[.02 0]);set(gca,'LineWidth',3,'fontsize',18,'fontweight','bold');
clear all
printf("\n\tDone, type [exit] or [quit] to leave.\n");
