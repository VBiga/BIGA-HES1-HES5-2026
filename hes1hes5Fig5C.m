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
%% solve the DDEs when production of HES5 is increasing
opt=ddeset('RelTol',1e-5);
values=[0:0.25:5];
thresh=[50:50:500];
for i=1:numel(values)
      par.ap5=values(i);
      par.am5=values(i);
      for j=1:numel(thresh)
          par.p0_51=thresh(j);
          sol = dde23(@ddePM_coupled,tau,[0 1 0 1],[0, Tend],opt,par);
          auc80(i,j)=sum(sol.y(1,:));
      end
end
%%
opt=ddeset('RelTol',1e-5);
values=[0:0.25:5];
thresh=[50:50:500];
par.dp5=log(2)/40; % protein degradation
for i=1:numel(values)
      par.ap5=values(i);
      par.am5=values(i);
      for j=1:numel(thresh)
          par.p0_51=thresh(j);
          sol = dde23(@ddePM_coupled,tau,[0 1 0 1],[0, Tend],opt,par);
          auc40(i,j)=sum(sol.y(1,:));
      end
end
%%
opt=ddeset('RelTol',1e-5);
values=[0:0.25:5];
thresh=[50:50:500];
par.dp5=log(2)/20; % protein degradation
for i=1:numel(values)
      par.ap5=values(i);
      par.am5=values(i);
      for j=1:numel(thresh)
          par.p0_51=thresh(j);
          sol = dde23(@ddePM_coupled,tau,[0 1 0 1],[0, Tend],opt,par);
          auc20(i,j)=sum(sol.y(1,:));
      end
end
%% normalise to max area under the curve
m1=max(max(auc20));
m2=max(max(auc40));
m3=max(max(auc80));
m=max([m1,m2,m3]);
auc20=auc20/m;
auc40=auc40/m;
auc80=auc80/m;
figure,subplot(1,3,1),imagesc(flipud(auc20))
set(gca,'YTick',[1:numel(values)]);
set(gca,'YTickLabel',[5:-0.25:0]);
set(gca,'XTick',[1:10]);
set(gca,'XTickLabel',thresh);%[500:-50:5]);
ylabel('Production rate');
xlabel('Cross-repression threshold')
title('Half-life 20min')
colorbar, clim([0,1]);
subplot(1,3,2),imagesc(flipud(auc40))
set(gca,'YTick',[1:numel(values)]);
set(gca,'YTickLabel',[5:-0.25:0]);
set(gca,'XTick',[1:10]);
set(gca,'XTickLabel',thresh);%[500:-50:5]);
ylabel('Production rate');
xlabel('Cross-repression threshold')
title('Half-life 40min')
colorbar, clim([0,1]);
subplot(1,3,3),imagesc(flipud(auc80))
set(gca,'YTick',[1:numel(values)]);
set(gca,'YTickLabel',[5:-0.25:0]);
set(gca,'XTick',[1:10]);
set(gca,'XTickLabel',thresh);%[500:-50:5]);
ylabel('Production rate');
xlabel('Cross-repression threshold')
title('Half-life 80min')
colorbar, clim([0,1]);
colormap bone
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
