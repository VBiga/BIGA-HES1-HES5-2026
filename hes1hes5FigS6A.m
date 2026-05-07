clear, close all, clc
Tend = 900;
tau = 29;
%HES1 parameters
par.ap1=1; % protein production
par.dp1=log(2)/22; % protein degradation
par.am1=1; %mrna production
par.dm1=log(2)/25; %mrna degradation
% hill function params
par.n1=7;
par.p0_1=390;
%HES5 parameters
par.ap5=1; % protein production
par.dp5=log(2)/80; % protein degradation
par.am5=1; %mrna production
par.dm5=log(2)/25; %mrna degradation
% hill function params
par.n5=7;
par.p0_5=390;
% HES1-HES5 parameters
par.n15=5;
par.p0_15=390;
par.n51=5;
par.p0_51=390;
% change the production and degradation rates
%% 5 onto 1
opt=ddeset('RelTol',1e-5);
thresh=[50:50:300]; 
figure
for i=1:numel(thresh)
    par.p0_51=thresh(i);
    sol = dde23(@ddePM_coupled,tau,[0 1 0 1],[0, Tend],opt,par);
    plot(sol.x'/60,sol.y(1,:));
    hold on
end 
legend('50','100','150','200','250','300')
xlabel('Time (h)'); ylabel('HES1 protein');
title('Repression threshold')
set(gca,'YLim',[0,1400]);
%% 5 onto 1
par.p0_51=390;
opt=ddeset('RelTol',1e-5);
values=[2:5]; 
figure
for i=1:numel(values)
    par.ap5=values(i);
    par.am5=values(i);
    sol = dde23(@ddePM_coupled,tau,[0 1 0 1],[0, Tend],opt,par);
    plot(sol.x'/60,sol.y(1,:));
    hold on
end 
legend('2','3','4','5')
xlabel('Time (h)'); ylabel('HES1 protein');
title('Production rates')
set(gca,'YLim',[0,1400]);
%%
par.ap5=1;
par.am5=1;
opt=ddeset('RelTol',1e-5);
thresh=[50:50:300]; 
figure
for i=1:numel(thresh)
    par.p0_15=thresh(i);
    sol = dde23(@ddePM_coupled,tau,[0 1 0 1],[0, Tend],opt,par);
    plot(sol.x'/60,sol.y(3,:));
    hold on
end 
legend('50','100','150','200','250','300')
xlabel('Time (h)'); ylabel('HES5 protein');
title('Repression threshold')
set(gca,'YLim',[0,1400]);
%% 1 onto 5
par.p0_51=390;
par.p0_15=390;
opt=ddeset('RelTol',1e-5);
values=[20:10:50]; 
figure
for i=1:numel(values)
    par.ap1=values(i);
    par.am1=values(i);
    sol = dde23(@ddePM_coupled,tau,[0 1 0 1],[0, Tend],opt,par);
    plot(sol.x'/60,sol.y(3,:));
    hold on
end 
legend('20','30','40','50')
xlabel('Time (h)'); ylabel('HES5 protein');
title('Production rates')
set(gca,'YLim',[0,1400]);
%% functions go here
function dydt=ddePM_coupled(t,y,Z,par)
% HES1
ylag1 = Z(1); % delayed protein-HES1
ylag51= Z(3);
inv_hill1=1+(ylag1/par.p0_1).^par.n1;
inv_hill51=1+(ylag51/par.p0_51).^par.n51;
inv_hill_tot1=inv_hill1*inv_hill51;
dydt(1,1) = par.ap1*y(2)-par.dp1*y(1);% protein variation
dydt(2,1)= par.am1/inv_hill_tot1-par.dm1*y(2);%mrna variation
% HES5
ylag5 = Z(3); % delayed protein-HES5
ylag15=Z(1);
inv_hill5=1+(ylag5/par.p0_5).^par.n5;
inv_hill15=1+(ylag15/par.p0_15).^par.n15;
inv_hill_tot2=inv_hill5*inv_hill15;
dydt(3,1) = par.ap5*y(4)-par.dp5*y(3);% protein variation
dydt(4,1)= par.am5/inv_hill_tot2-par.dm5*y(4);%mrna variation
end
