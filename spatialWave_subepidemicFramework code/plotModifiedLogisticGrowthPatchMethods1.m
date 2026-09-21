function objfunction=plotsir(z)

global flag1 method1 timevect ydata I0 npatches onset_thr typedecline1 yfit


global invasions
global timeinvasions
global Cinvasions

global npatches_fixed

global onset_fixed


r=z(1);
p=z(2);
a=z(3);
K=z(4);
q=z(5);

alpha=z(end-1);
d=z(end);

%[r p a K q]

%if onset_thr>K
    
%    npatches=1;

%else


if npatches_fixed<0
    
    
    npatches2=numberSubepidemics(K,onset_thr,q,typedecline1);
    
    npatches2=max(npatches2,1);
    
    npatchesX=min(npatches,npatches2);
    
else
    
    npatchesX=npatches;
    
end

%end


if onset_thr>K
    
   npatchesX=1;
   
end

invasions=zeros(npatchesX,1);
timeinvasions=zeros(npatchesX,1);
Cinvasions=zeros(npatchesX,1);

invasions(1)=1;
timeinvasions(1)=0;
Cinvasions(1)=0;



IC=zeros(npatchesX,1);

IC(1,1)=I0;

if npatchesX>1
    IC(2:end,1)=1;
end

[t,x]=ode15s(@modifiedLogisticGrowthPatch,timevect,IC,[],r,p,a,K,npatchesX,onset_thr,q,flag1,typedecline1);


y=sum(x,2);


totinc=[y(1,1);diff(y(:,1))];

if onset_fixed==0
    totinc(1)=totinc(1)-(npatchesX-1);
end

yfit=totinc;


[objfunction,yfit]=getSpatialWaveObjective(ydata,yfit,method1,alpha,d);
end
