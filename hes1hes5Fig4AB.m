clear, close all, clc
Tend = 900;
% define the dde
% HES1 self-repression + HES5 repressing Hes1
% dp1/dt=ap1*m1(t)-dp1*p1(t);
% dm1/dt=am1*HIll1(p1(t-tau1),p0_1,n1)*Hill51(p5(t-tau51),po_51,n51)-dm1*m1(t);
% HES5 self-repression + HES1 repressing Hes5
% dp5/dt=ap5*m5(t)-dp5*p5(t);
% dm5/dt=am5*Hill(p5(t-tau5),p0_5,n5)*Hill15(p15(t-tau15),p0_15,n15)-dm5*m5(t);
tau = 29;
h1h5_tau=29;
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
par.p0_15=390000;
par.n51=5;
par.p0_51=390000;
%% figure 4A- uncoupled example; cross-rep set to a very high thresh
% solve the DDEs
opt=ddeset('RelTol',1e-5);
sol = dde23(@ddePM_coupled,tau,[0 1 0 1],[0, Tend],opt,par);
t = sol.x'/60;
figure,plot(t,sol.y(1,:),'color','r');
hold on
plot(t,sol.y(3,:),'color','b');
legend('HES1','HES5');
xlabel('Time (h)');
ylabel('Protein Level');
title ('Uncoupled (free-running)');
%% figure 4B- coupled at balanced values of cross-repression
par.p0_15=390;
par.p0_51=390;
sol = dde23(@ddePM_coupled,tau,[0 1 0 1],[0, Tend],opt,par);
t = sol.x'/60;
figure,plot(t,sol.y(1,:),'color','r');
hold on
plot(t,sol.y(3,:),'color','b');
legend('HES1','HES5');
xlabel('Time (h)');
ylabel('Protein Level');
title('Coupled HES1-HES5 (balanced)')
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

