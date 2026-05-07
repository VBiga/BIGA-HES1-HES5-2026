clear, close all, clc
Tend = 1200;
tau = 29;
%HES1 parameters
par.ap1=1; % protein production
par.dp1=log(2)/22; % protein degradation
par.am1=1; %mrna production
par.dm1=log(2)/25; %mrna degradation
% hill function params- both Hill set to 7
par.n1=7;
par.p0_1=390;
%HES5 parameters
par.ap5=1; % protein production
par.dp5=log(2)/80; % protein degradation
par.am5=1; %mrna production
par.dm5=log(2)/25; %mrna degradation
par.n5=7;
par.p0_5=390;
% HES1-HES5 parameters
par.n15=5;
par.p0_15=390;
par.n51=5;
par.p0_51=390;
% set colors
col=parula(800);
%%  solve dde
opt=ddeset('RelTol',1e-5);
sol = dde23(@ddePM_coupled,tau,[0 1 0 1],[0, Tend],opt,par);
t=sol.x/60;
figure, subplot(3,2,1)
plot(t,sol.y(1,:),'color','r');
hold on
plot(t,sol.y(3,:),'color','b');
legend('HES1','HES5')
set(gca,'XLim',[0,15]);
% phase plot
subplot(3,2,2)
ph1=getHilbert(sol.y(1,:));
ph5=getHilbert(sol.y(3,:)); 
% cut to 15h
idx=find(t>15,1,'first');
ph1(idx:end)=[];
ph5(idx:end)=[];
for k=1:numel(ph1)
    plot(ph1(k),ph5(k),'.','color',col(k,:)); hold on
end
set(gca,'XLim',[-pi,pi],'YLim',[-pi,pi]);
set(gca,'YTick',[-2,0,2]);
axis square
xlabel('HES1 phase');
ylabel('HES5 phase');
%plot([-pi,pi],[-pi, pi],'k')
plot([-pi/2,pi],[-pi, pi/2],'k--');
plot([-pi,pi/2],[-pi/2, pi],'k--');
plot([pi/2,pi],[-pi,-pi/2],'k--')
plot([-pi,-pi/2],[pi/2,pi],'k--')
hold on
plot(ph1(1),ph5(1),'k*');
plot(ph1(end),ph5(end),'ko');
colormap parula
colorbar
%% initial conditions
sol = dde23(@ddePM_coupled,tau,[250 17 0 1],[0, Tend],opt,par); 
t=sol.x/60;
subplot(3,2,3),plot(t,sol.y(1,:),'color','r');
hold on
plot(t,sol.y(3,:),'color','b');
set(gca,'XLim',[0,15]);
legend('HES1','HES5')
ph1=getHilbert(sol.y(1,:));
ph5=getHilbert(sol.y(3,:)); 
% cut to 15h
idx=find(t>15,1,'first');
ph1(idx:end)=[];
ph5(idx:end)=[];
subplot(3,2,4)
for k=1:numel(ph1)
    plot(ph1(k),ph5(k),'.','color',col(k,:)); hold on
end
set(gca,'XLim',[-pi,pi],'YLim',[-pi,pi]);
set(gca,'YTick',[-2,0,2]);
axis square
xlabel('HES1 phase');
ylabel('HES5 phase');
%plot([-pi,pi],[-pi, pi],'k')
plot([-pi/2,pi],[-pi, pi/2],'k--');
plot([-pi,pi/2],[-pi/2, pi],'k--');
plot([pi/2,pi],[-pi,-pi/2],'k--')
plot([-pi,-pi/2],[pi/2,pi],'k--')
plot(ph1(1),ph5(1),'k*');
plot(ph1(end),ph5(end),'ko');
colorbar
%% initial conditions
sol = dde23(@ddePM_coupled,tau,[500,17, 0, 1],[0, Tend],opt,par);
t=sol.x/60;
subplot(3,2,5),plot(t,sol.y(1,:),'color','r');
hold on
plot(t,sol.y(3,:),'color','b');
set(gca,'XLim',[0,15]);
legend('HES1','HES5')
subplot(3,2,6)
ph1=getHilbert(sol.y(1,:));
ph5=getHilbert(sol.y(3,:)); 
% cut to 15h
idx=find(t>15,1,'first');
ph1(idx:end)=[];
ph5(idx:end)=[];
for k=1:numel(ph1)
    plot(ph1(k),ph5(k),'.','color',col(k,:)); hold on
end
set(gca,'XLim',[-pi,pi],'YLim',[-pi,pi]);
set(gca,'YTick',[-2,0,2]);
axis square
xlabel('HES1 phase');
ylabel('HES5 phase');
%plot([-pi,pi],[-pi, pi],'k')
plot([-pi/2,pi],[-pi, pi/2],'k--');
plot([-pi,pi/2],[-pi/2, pi],'k--');
plot([pi/2,pi],[-pi,-pi/2],'k--')
plot([-pi,-pi/2],[pi/2,pi],'k--')
plot(ph1(1),ph5(1),'k*');
plot(ph1(end),ph5(end),'ko');
colorbar
%% HES5 vary initial cond
sol = dde23(@ddePM_coupled,tau,[0 1 500 27],[0, Tend],opt,par); 
t=sol.x/60;
figure,subplot(3,2,1),plot(t,sol.y(1,:),'color','r');
hold on
plot(t,sol.y(3,:),'color','b');
set(gca,'XLim',[0,15]);
legend('HES1','HES5')
subplot(3,2,2)
ph1=getHilbert(sol.y(1,:));
ph5=getHilbert(sol.y(3,:)); 
% cut to 15h
idx=find(t>15,1,'first');
ph1(idx:end)=[];
ph5(idx:end)=[];
n=min(numel(ph1),numel(ph5));
ph1=ph1(1:n);
ph5=ph5(1:n);
for k=1:numel(ph1)
    plot(ph1(k),ph5(k),'.','color',col(k,:)); hold on
end
set(gca,'XLim',[-pi,pi],'YLim',[-pi,pi]);
set(gca,'YTick',[-2,0,2]);
axis square
xlabel('HES1 phase');
ylabel('HES5 phase');
%plot([-pi,pi],[-pi, pi],'k')
plot([-pi/2,pi],[-pi, pi/2],'k--');
plot([-pi,pi/2],[-pi/2, pi],'k--');
plot([pi/2,pi],[-pi,-pi/2],'k--')
plot([-pi,-pi/2],[pi/2,pi],'k--')
plot(ph1(1),ph5(1),'k*');
plot(ph1(end),ph5(end),'ko');
colormap(parula)
colorbar
%%
sol = dde23(@ddePM_coupled,tau,[250 17 500 27],[0, Tend],opt,par); 
t=sol.x/60;
subplot(3,2,3),plot(t,sol.y(1,:),'color','r');
hold on
plot(t,sol.y(3,:),'color','b');
set(gca,'XLim',[0,15]);
legend('HES1','HES5')
subplot(3,2,4)
ph1=getHilbert(sol.y(1,:));
ph5=getHilbert(sol.y(3,:)); 
idx=find(t>15,1,'first');
n=min([numel(ph1),numel(ph5),idx]);
ph1=ph1(1:n);
ph5=ph5(1:n);
for k=1:numel(ph1)
    plot(ph1(k),ph5(k),'.','color',col(k,:)); hold on
end
set(gca,'XLim',[-pi,pi],'YLim',[-pi,pi]);
set(gca,'YTick',[-2,0,2]);
axis square
xlabel('HES1 phase');
ylabel('HES5 phase');
%plot([-pi,pi],[-pi, pi],'k')
plot([-pi/2,pi],[-pi, pi/2],'k--');
plot([-pi,pi/2],[-pi/2, pi],'k--');
plot([pi/2,pi],[-pi,-pi/2],'k--')
plot([-pi,-pi/2],[pi/2,pi],'k--')
plot(ph1(1),ph5(1),'k*');
plot(ph1(end),ph5(end),'ko');
colormap(parula)
colorbar
%%
sol = dde23(@ddePM_coupled,tau,[1000 17 500 27],[0, Tend],opt,par); 
t=sol.x/60;
subplot(3,2,5),plot(t,sol.y(1,:),'color','r');
hold on
plot(t,sol.y(3,:),'color','b');
set(gca,'XLim',[0,15]);
legend('HES1','HES5')
subplot(3,2,6)
ph1=getHilbert(sol.y(1,:));
ph5=getHilbert(sol.y(3,:)); 
idx=find(t>15,1,'first');
n=min([numel(ph1),numel(ph5),idx]);
ph1=ph1(1:n);
ph5=ph5(1:n);
for k=1:numel(ph1)
    plot(ph1(k),ph5(k),'.','color',col(k,:)); hold on
end
set(gca,'XLim',[-pi,pi],'YLim',[-pi,pi]);
set(gca,'YTick',[-2,0,2]);
axis square
xlabel('HES1 phase');
ylabel('HES5 phase');
%plot([-pi,pi],[-pi, pi],'k')
plot([-pi/2,pi],[-pi, pi/2],'k--');
plot([-pi,pi/2],[-pi/2, pi],'k--');
plot([pi/2,pi],[-pi,-pi/2],'k--')
plot([-pi,-pi/2],[pi/2,pi],'k--')
plot(ph1(1),ph5(1),'k*');
plot(ph1(end),ph5(end),'ko');
colormap(parula)
colorbar
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
function ph=getHilbert(sig)
sig=sig-mean(sig);
hb=hilbert(sig);
ph=angle(hb);
% throw away last cycle-incomplete
idx=find(ph>0,1,'last');
ph(idx:end)=[];
end