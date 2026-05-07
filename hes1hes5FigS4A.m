clear, close all, clc
Tend = 900;
tau = 29;
%HES1 parameters
par.ap1=1; % protein production
par.dp1=log(2)/22; % protein degradation
par.am1=1; %mrna production
par.dm1=log(2)/25; %mrna degradation
% hill function params
par.n1=5;
par.p0_1=390;
%HES5 parameters
par.ap5=1; % protein production
par.dp5=log(2)/80; % protein degradation
par.am5=1; %mrna production
par.dm5=log(2)/25; %mrna degradation
% hill function params
par.n5=8;
par.p0_5=390;
% HES1-HES5 parameters
par.n15=5;
par.p0_15=390000;
par.n51=5;
par.p0_51=390000;
%% solve the DDEs
par.n5=5;
par.n1=5;
opt=ddeset('RelTol',1e-5);
sol = dde23(@ddePM_coupled,tau,[0 1 0 1],[0, Tend],opt,par);
t = sol.x'/60;
plot(t,sol.y(1,:));
hold on
Perh1(1)=getPeriod(t,sol.y(1,:));
% 
par.n1=6;
sol = dde23(@ddePM_coupled,tau,[0 1 0 1],[0, Tend],opt,par);
t = sol.x'/60;
plot(t,sol.y(1,:));
Perh1(2)=getPeriod(t,sol.y(1,:));
%
par.n1=7;
sol = dde23(@ddePM_coupled,tau,[0 1 0 1],[0, Tend],opt,par);
t = sol.x'/60;
plot(t,sol.y(1,:));
Perh1(3)=getPeriod(t,sol.y(1,:));
xlabel ('Time (h)');
legend('n=5', 'n=6', 'n=7')
set(gca,'YLim',[0,1400])
title ('HES1 at different Hill coefficients')
%%
par.n1=5;
par.n5=5;
opt=ddeset('RelTol',1e-5);
sol = dde23(@ddePM_coupled,tau,[0 1 0 1],[0, Tend],opt,par);
t = sol.x'/60;
figure,plot(t,sol.y(3,:));
hold on
Per5(1)=getPeriod(t,sol.y(3,:));
% 
par.n5=6;
sol = dde23(@ddePM_coupled,tau,[0 1 0 1],[0, Tend],opt,par);
t = sol.x'/60;
plot(t,sol.y(3,:));
Per5(2)=getPeriod(t,sol.y(3,:));
%
par.n5=7;
sol = dde23(@ddePM_coupled,tau,[0 1 0 1],[0, Tend],opt,par);
t = sol.x'/60;
plot(t,sol.y(3,:));
Per5(3)=getPeriod(t,sol.y(3,:));
xlabel ('Time (h)');
legend('n=5', 'n=6', 'n=7')
title ('HES5 at different Hill coefficients')
set(gca,'YLim',[0,1400])
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

function y=getPeriod(t,hes)
[pks,locs]=findpeaks(hes);
y=mean(t(locs(2:end))-t(locs(1:end-1)))
end