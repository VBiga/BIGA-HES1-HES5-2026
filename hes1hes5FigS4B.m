clear, close all, clc
% Deterministic model single cell (Fig3B)
%DefaultParams_nonDim
tau = 29;
Tend = 800;
% define the dde
% dp/dt=ap*m-dp*p;
% dm/dt=am*HIll(p,p0,n)-dm*m;
par.ap=1; % protein production
par.dp=log(2)/80; % protein degradation
par.am=1; %mrna production
par.dm=log(2)/30; %mrna degradation % was 20
% hill function params
par.n=5;
par.p0=390;
% define the parameters of the model
opt=ddeset('RelTol',1e-5);
sol = dde23(@ddePM,tau,[0 1],[0, Tend],opt,par);
%% explore combinations and identify ones that achieve perfect adaptation
Mat=[];
mp=10:1:80;
pp=10:1:80;
for i=1:numel(mp)
    %
    for j=1:numel(pp)
        par.dm=log(2)/mp(i);
        par.dp=log(2)/pp(j);
        sol = dde23(@ddePM,tau,[0 1],[0, Tend],opt,par);
        Mat{i,j}=sol.y(1,:);
        % compute the frequency
        vect=sol.y(1,:);
        t=sol.x/60;
        [pks,locs]=findpeaks(vect);
        if numel(locs)>2
            locs=locs(1:3);
            per=mean(t(locs(2:end))-t(locs(1:end-1)));
        else
            per=t(locs(2:end))-t(locs(1:end-1));
        end
        Per(i,j)=per;
    end
end
%% display as colormap
figure,imagesc(flipud(Per))
colormap hsv
set(gca,'XTick',[1:5:71]);
set(gca,'XTickLabel',[10:5:80]);
set(gca,'YTick',[1:5:71]);
set(gca,'YTickLabel',[80:-5:1]);
xlabel('mRNA half-life');
ylabel('protein half-life');
colorbar
%% display as contour plot
figure,contour(Per)
set(gca,'XTick',[1:5:71]);
set(gca,'YTick',[1:5:71]);
set(gca,'XTickLabel',[10:5:80]);
set(gca,'YTickLabel',[10:5:80]);
xlabel('mRNA half-life');
ylabel('protein half-life');
%% functions go here
function dydt=ddePM(t,y,Z,par)
    ylag1 = Z(1); % delayed protein
    inv_hill=1+(ylag1/par.p0).^par.n;
    dydt(1,1) = par.ap*y(2)-par.dp*y(1);% protein variation
    dydt(2,1)= par.am/inv_hill-par.dm*y(2);%mrna variation
end
