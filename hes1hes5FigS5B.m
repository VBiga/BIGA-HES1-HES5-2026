clear, close all, clc
Tend = 1200;
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
par.p0_15=390;
par.n51=5;
par.p0_51=390;
%% solve the DDEs when HES5 onto HES1 is changing rep threshold
opt=ddeset('RelTol',1e-5);
values=[200:20:1000];
for i=1:numel(values)
    for j=1:numel(values)
      par.p0_15=values(i);
      par.p0_51=values(j);
      sol = dde23(@ddePM_coupled,tau,[0 1 0 1],[0, Tend],opt,par);
     % measure periodicity
      t = sol.x'/60;
      hes1=sol.y(1,:);
      hes5=sol.y(3,:);
      per1(i,j)=getPeriod(t,hes1);
      per2(i,j)=getPeriod(t,hes5);
      levmax1(i,j)=max(hes1);
      levmax2(i,j)=max(hes5);
      idx=find(t>7,1,'first'); % ignore first peak
      lev1(i,j)=mean(hes1(idx:end));
      lev2(i,j)=mean(hes5(idx:end));
    end
end
%%
len=numel(values);
uData = unique([per1(:);per2(:)]);
cmap = parula(numel(uData));
figure,imagesc(flipud(per1)),colormap(cmap)
colorbar, clim([2,5]);
xlabel('P5 to 1');
ylabel('P1 to 5');
set(gca,'XTick',1:5:len);
set(gca,'XTickLabel',values(1:5:len));
set(gca,'YTick',1:5:len);
set(gca,'YTickLabel',values(len:-5:1));
title('HES1 Period (h)')
%%
figure,imagesc(flipud(per2)),colormap(cmap)
colorbar, clim([2,5]);
xlabel('P5 to 1');
ylabel('P1 to 5');
set(gca,'XTick',1:5:len);
set(gca,'XTickLabel',values(1:5:len));
set(gca,'YTick',1:5:len);
set(gca,'YTickLabel',values(len:-5:1));
title('HES5 Period (h)')
%% difference between the period values
perdiff=abs(per1-per2);
figure,imagesc(flipud(perdiff))
colormap("bone")
colorbar, clim([0,2.5])
xlabel('P5 to 1');
ylabel('P1 to 5');
set(gca,'XTick',1:5:len);
set(gca,'XTickLabel',values(1:5:len));
set(gca,'YTick',1:5:len);
set(gca,'YTickLabel',values(len:-5:1));
title('Difference (absolute values)')
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
locs(1)=[]; % remove first peak
y=mean(t(locs(2:end))-t(locs(1:end-1)));
end
