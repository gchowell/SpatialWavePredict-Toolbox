% Plot model fits and derive performance metrics during the calibration period for the best fitting models

clear
clear global

close all

% <============================================================================>
% <=================== Declare global variables =======================================>
% <============================================================================>

global invasions
global timeinvasions
global Cinvasions
global npatches_fixed
global onset_fixed

global method1 dist1 factor1 smoothfactor1 calibrationperiod1


% <============================================================================>
% <=================== Load parameter values supplied by user =================>
% <============================================================================>

% options.m
[outbreakx_INP, caddate1_INP, cadregion_INP, caddisease_INP,datatype_INP, DT_INP,datafilename1_INP, datevecfirst1_INP, datevecend1_INP, numstartpoints_INP,topmodelsx_INP, M_INP,flag1_INP,typedecline2_INP]=options

% <============================================================================>
% <================================ Dataset ======================================>
% <============================================================================>

outbreakx=outbreakx_INP;
caddate1=caddate1_INP;

cadregion=cadregion_INP; % string indicating the region of the time series (USA, Chile, Mexico, Nepal, etc)

caddisease=caddisease_INP;

datatype=datatype_INP;

DT=DT_INP; % temporal resolution in days (1=daily data, 7=weekly data).

if DT==1
    cadtemporal='daily';
elseif DT==7
    cadtemporal='weekly';
end

cadfilename2=strcat(cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-state-',num2str(outbreakx),'-',caddate1);


% <============================================================================>
% <============================Adjustments to data =================================>
% <============================================================================>

%smoothfactor1=1; % <smoothfactor1>-day rolling average of the case series

%calibrationperiod1=90; % calibrates model using the most recent <calibrationperiod1> days  where calibrationperiod1<length(data1)

% <=============================================================================>
% <=========================== Statistical method ====================================>
% <=============================================================================>

%method1=0; %Type of estimation method: 0 = nonlinear least squares

%dist1=0; % Normnal distribution to model error structure


% <==============================================================================>
% <========================= Growth model ==========================================>
% <==============================================================================>

npatchess2=npatches_fixed;  % maximum number of subepidemics considered in epidemic trajectory fit

%GGM=0;  % 0 = GGM
%GLM=1;  % 1 = GLM
%GRM=2;  % 2 = GRM
%LM=3;   % 3 = LM
%RICH=4; % 4 = Richards

flagss2=flag1_INP; % Growth model considered in epidemic trajectory


% <==============================================================================>
% <======== Number of best fitting models used to generate ensemble model ========================>
% <==============================================================================>

topmodels1=1:topmodelsx_INP;

MAES=[];
MSES=[];
PIS=[];
MISS=[];

epiwavesizes=[];
param_K0s=[];
param_qs=[];
numsubepidemicss=[];

cc2=1;

AICc_bests=[];

for rank1=1
    
    cc2=1;
    
    npatches_fixed=npatchess2
    
    
    flag1=flagss2;
    
    npatchess2
    
    % <========================================================================================>
    % <================================ Load model results ====================================>
    % <========================================================================================>
    
    load (strcat('./output/modifiedLogisticPatch-original-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-0-smoothing-',num2str(smoothfactor1),'-',cadfilename2,'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'-rank-',num2str(rank1),'.mat'))
    
    rank1
    
    npatches
    
    AICc_bests=[AICc_bests;AICc_best];
    
    timevect=(data1(:,1))*DT;
    
    % <========================================================================================>
    % <================================ Parameter estimates =========================================>
    % <========================================================================================>
    
    
    % Parameter values
    param_r=[mean(Phatss(:,1)) quantile(Phatss(:,1),0.025) quantile(Phatss(:,1),0.975)];
    
    param_p=[mean(Phatss(:,2)) quantile(Phatss(:,2),0.025) quantile(Phatss(:,2),0.975)];
    
    param_a=[mean(Phatss(:,3)) quantile(Phatss(:,3),0.025) quantile(Phatss(:,3),0.975)];
    
    param_K=[mean(Phatss(:,4)) quantile(Phatss(:,4),0.025) quantile(Phatss(:,4),0.975)];
    
    param_q=[mean(Phatss(:,5)) quantile(Phatss(:,5),0.025) quantile(Phatss(:,5),0.975)];
    
    param_alpha=[mean(Phatss(:,6)) quantile(Phatss(:,6),0.025) quantile(Phatss(:,6),0.975)];
    
    param_d=[mean(Phatss(:,7)) quantile(Phatss(:,7),0.025) quantile(Phatss(:,7),0.975)];
    
    
    
    cad1=strcat('r=',num2str(param_r(end,1),2),'(95%CI:',num2str(param_r(end,2),2),',',num2str(param_r(end,3),2),')');
    cad2=strcat('p=',num2str(param_p(end,1),2),'(95%CI:',num2str(param_p(end,2),2),',',num2str(param_p(end,3),2),')')
    cad3=strcat('a=',num2str(param_a(end,1),2),'(95%CI:',num2str(param_a(end,2),2),',',num2str(param_a(end,3),2),')')
    cad4=strcat('K=',num2str(param_K(end,1),3),'(95%CI:',num2str(param_K(end,2),3),',',num2str(param_K(end,3),3),')')
    cad5=strcat('q=',num2str(param_q(end,1),3),'(95%CI:',num2str(param_q(end,2),3),',',num2str(param_q(end,3),3),')')
    
    cad6=strcat('alpha=',num2str(param_alpha(end,1),3),'(95%CI:',num2str(param_alpha(end,2),3),',',num2str(param_alpha(end,3),3),')')
    
    cad7=strcat('d=',num2str(param_d(end,1),3),'(95%CI:',num2str(param_d(end,2),3),',',num2str(param_d(end,3),3),')')
    
    
  

    close all
    figure(100)
    timevect=(data1(:,1))*DT;
    
    subplot(2,4,1)
    hist(Phatss(:,1))
    xlabel('r')
    ylabel('Frequency')
    title(cad1)
    
    hold on
    
    line2=[param_r(1,2) 10;param_r(1,3) 10];
    
    line1=plot(line2(:,1),line2(:,2),'r--')
    set(line1,'LineWidth',2)
    
    axis([param_r(1,2)-0.05 param_r(1,3)+0.05 0 120])
    
    set(gca,'FontSize', 16);
    set(gcf,'color','white')
    
    
    subplot(2,4,2)
    hist(Phatss(:,2))
    xlabel('p')
    hold on
    title(cad2)
    
    line2=[param_p(1,2) 10;param_p(1,3) 10];
    
    line1=plot(line2(:,1),line2(:,2),'r--')
    set(line1,'LineWidth',2)
    
    axis([param_p(1,2)-0.025 param_p(1,3)+0.025 0 120])
    
    set(gca,'FontSize', 16);
    set(gcf,'color','white')
    
    
    subplot(2,4,3)
    hist(Phatss(:,4))
    xlabel('K')
    hold on
    title(cad4)
    
    line2=[param_K(1,2) 10;param_K(1,3) 10];
    
    line1=plot(line2(:,1),line2(:,2),'r--')
    set(line1,'LineWidth',2)
    
    axis([param_K(1,2)-50 param_K(1,3)+50 0 120])
    
    set(gca,'FontSize', 16);
    set(gcf,'color','white')
    
    subplot(2,4,4)
    hist(Phatss(:,5))
    xlabel('q')
    hold on
    title(cad5)
    
    line2=[param_q(1,2) 10;param_q(1,3) 10];
    
    line1=plot(line2(:,1),line2(:,2),'r--')
    set(line1,'LineWidth',2)
    
    axis([max(param_q(1,2)-0.1,0) param_q(1,3)+0.1 0 200])
    
    set(gca,'FontSize', 16);
    set(gcf,'color','white')
    
    
    % plot model fit and uncertainty
    
    
    subplot(2,4,5)
    %plot(timevect,curves,'c-')
    %hold on
    
    LB1=quantile(curves',0.025);
    UB1=quantile(curves',0.975);
    
    
    h=area(timevect',[LB1' UB1'-LB1'])
    hold on
    
    h(1).FaceColor = [1 1 1];
    h(2).FaceColor = [0.8 0.8 0.8];
    
    
    
    %line1=plot(timevect,mean(curves,2),'k--')  % asymptotic mean
    line1=plot(timevect,bestfit,'r-')  %LSQ fit
    
    set(line1,'LineWidth',2)
    
    line1=plot(timevect,quantile(curves',0.025),'r--')
    set(line1,'LineWidth',2)
    
    line1=plot(timevect,quantile(curves',0.975),'r--')
    set(line1,'LineWidth',2)
    
    line1=plot(timevect,data1(:,2),'ko')
    set(line1,'LineWidth',2)
    
    hold on
    
    
    % 95% prediction coverage
    coverage1=find(data1(:,2)>=quantile(curves',0.025)' & data1(:,2)<=quantile(curves',0.975)');
    
    
    xlabel('Time (days)')
    ylabel(strcat(caddisease,{' '},datatype))
    
    axis([timevect(1) timevect(end)+1 0 max(data(:,2))*1.3])
    
    set(gca,'FontSize', 16);
    set(gcf,'color','white')
    
    %goodness of fit
    
    MAE=mean(abs(bestfit-data1(:,2)))
    
    %SSE=sum((mean(curves,2)-data1(:,2)).^2);
    MSE=mean((bestfit-data1(:,2)).^2)
    
    RMSE=sqrt(MSE)
    
    Lt=quantile(curves',0.025)';
    Ut=quantile(curves',0.975)';
    
    coverage1=find(data1(:,2)>=Lt & data1(:,2)<=Ut);
    
    'prediction coverage (%)='
    PI=100*(length(coverage1)./length(data1(:,2)))
    
    MIS=mean([(Ut-Lt)+(2/0.05).*(Lt-data1(:,2)).*(data1(:,2)<Lt)+(2/0.05)*(data1(:,2)-Ut).*(data1(:,2)>Ut)])
    
    
    % residuals
    
    subplot(2,4,[7 8])
    
    resid1=bestfit-data1(:,2);
    
    stem(timevect,resid1,'b')
    hold on
    
    %anscomberesid=(3/2)*(data1(:,2).^(2/3)-bestfit.^(2/3))./(bestfit.^(1/6))
    %stem(timevect,anscomberesid,'r')
    %hold on
    
    %plot(timevect,zeros(length(timevect),1)+1.96,'k--')
    %plot(timevect,zeros(length(timevect),1)-1.96,'k--')
    
    xlabel('Time (days)')
    ylabel('Residuals')
    
    axis([timevect(1) timevect(end)+1 min(resid1)-1 max(resid1)+1])
    
    set(gca,'FontSize', 16);
    set(gcf,'color','white')
    
    
    % histogram of the predicted distribution of the number of subepidemics and epidemic
    % wave size
    
    numsubepidemics1=[];
    
    totepisize1=[];
    
    for i=1:M
        
        [Ks,Ktot,n]=getsubepidemicsizesFinite(Phatss(i,4),onset_thr,Phatss(i,5),typedecline1);
        
        %n=Phatss(i,8);
        
        totepisize1=[totepisize1;[sum(Ks)]];
        
        
    end
    
    %
    %     if typedecline1==2
    %
    %         index1=find(totepisize1(:,2)>0);
    %
    %         totepisize1=totepisize1(index1,:);
    %
    %     end
    
    
    
    
    %totepisize1
    
    %'num subepidemics'
    %param1=[mean(Phatss(:,8)) quantile(Phatss(:,8),0.025) quantile(Phatss(:,8),0.975)];
    
    
    'epidemic wave size'
    param2=[mean(totepisize1(:,1)) quantile(totepisize1(:,1),0.025) quantile(totepisize1(:,1),0.975)];
    
    
    
    figure(101)
    
%     subplot(1,2,1)
%     
%     hist(Phatss(:,6))
%     hold on
%     
%     line2=[param1(1,2) 20;param1(1,3) 20];
%     
%     line1=plot(line2(:,1),line2(:,2),'r--')
%     set(line1,'LineWidth',2)
%     
%     line3=[param1(1,1) 0;param1(1,1) 250];
%     
%     line1=plot(line3(:,1),line3(:,2),'r--')
%     set(line1,'LineWidth',2)
%     
%     
%     
%     xlabel('Number of sub-epidemics')
%     ylabel('Frequency')
%     
%     set(gca,'FontSize', 16);
%     set(gcf,'color','white')
%     
%     
    subplot(1,2,2)
    
    hist(totepisize1(:,1))
    
    hold on
    
    line2=[param2(1,2) 20;param2(1,3) 20];
    
    line1=plot(line2(:,1),line2(:,2),'r--')
    set(line1,'LineWidth',2)
    
    line3=[param2(1,1) 0;param2(1,1) 250];
    
    line1=plot(line3(:,1),line3(:,2),'r--')
    set(line1,'LineWidth',2)
    
    lineEpiSize=[sum(data(:,2)) 0;sum(data(:,2)) 250];
    
    line1=plot(lineEpiSize(:,1),lineEpiSize(:,2),'r-+')
    set(line1,'LineWidth',2)
    
    cad1=strcat('Epidemic wave size=',num2str(param2(end,1),2),'(95%CI:',num2str(param2(end,2),2),',',num2str(param2(end,3),2),')');
    
    title(cad1)
    
    xlabel('Total epidemic wave size')
    ylabel('Frequency')
    
    set(gca,'FontSize', 16);
    set(gcf,'color','white')
    

    'epidemic wave size'
    param2
    
    'actual epidemic size'
    sum(data(:,2))
    
    MAES=[MAES;MAE];
    MSES=[MSES;MSE];
    PIS=[PIS;PI];
    MISS=[MISS;MIS];
    
    
    epiwavesizes=[epiwavesizes;param2];
    param_K0s=[param_K0s;param_K];
    param_qs=[param_qs;param_q];
    
end



figure
subplot(2,2,1)
line1=plot(MAES,'b-o')
set(line1,'linewidth',2)
xlabel('Number of sub-epidemics')
ylabel('MAE')

set(gca,'FontSize', 16);
set(gcf,'color','white')

subplot(2,2,2)
line1=plot(MSES,'b-o')
set(line1,'linewidth',2)
xlabel('Number of sub-epidemics')
ylabel('MSE')

set(gca,'FontSize', 16);
set(gcf,'color','white')

subplot(2,2,3)
line1=plot(PIS,'b-o')
set(line1,'linewidth',2)
xlabel('Number of sub-epidemics')
ylabel('PI')

set(gca,'FontSize', 16);
set(gcf,'color','white')

subplot(2,2,4)
line1=plot(MISS,'b-o')
set(line1,'linewidth',2)
xlabel('Number of sub-epidemics')
ylabel('MIS')

set(gca,'FontSize', 16);
set(gcf,'color','white')


% <============================================================================>
% <========================plot sub-epidemic profile ==========================>
% <============================================================================>

figure(100)

color1=['r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';];


% generate forecast curves from each bootstrap realization
for realization=1:M
    
    
    r_hat=Phatss(realization,1);
    p_hat=Phatss(realization,2);
    a_hat=Phatss(realization,3);
    K_hat=Phatss(realization,4);
    q_hat=Phatss(realization,5);
    
    
    alpha_hat=Phatss(realization,6);
    d_hat=Phatss(realization,7);
    
    
    if method1==3
        
        dist1=3; % VAR=mean+alpha*mean;
        
        factor1=alpha_hat;
        
    elseif method1==4
        
        dist1=4; % VAR=mean+alpha*mean^2;
        
        factor1=alpha_hat;
        
    elseif method1==5
        
        dist1=5; % VAR=mean+alpha*mean^2;
        
        factor1=alpha_hat;
        
    end


    %npatches=Phatss(realization,8);
    
    IC=zeros(npatches,1);
    
    IC(1,1)=data1(1,2);
    IC(2:end,1)=1;
    
    invasions=zeros(npatches,1);
    timeinvasions=zeros(npatches,1);
    Cinvasions=zeros(npatches,1);
    
    invasions(1)=1;
    timeinvasions(1)=0;
    Cinvasions(1)=0;
    
    
    [~,x]=ode15s(@modifiedLogisticGrowthPatch,timevect,IC,[],r_hat,p_hat,a_hat,K_hat,npatches,onset_thr,q_hat,flag1,typedecline1);

    %[~,x]=ode23s(@modifiedLogisticGrowthPatch,timevect,IC,[],rs_hat,ps_hat,as_hat,Ks_hat,npatches,onset_thr,flag1,typefit1);
    
    subplot(2,4,6)
    
    for j=1:npatches
        
        incidence1=[x(1,j);diff(x(:,j))];
        
        plot(timevect,incidence1,color1(j,:))
        
        hold on
        
    end
    
    y=sum(x,2);
    
    totinc=[y(1,1);diff(y(:,1))];
    
    totinc(1)=totinc(1)-(npatches-1);
    
    bestfit=totinc;
    
    gray1=gray(10);
    
    plot(timevect,totinc,'color',gray1(7,:))
    %hold on
    plot(timevect,data1(:,2),'ko')
    
end

xlabel('\fontsize{16}Time (days)');
ylabel(strcat(caddisease,{' '},datatype))

line1=plot(data(:,1)*DT,data(:,2),'ko')
set(line1,'LineWidth',2)

axis([timevect(1) timevect(end)+1 0 max(data(:,2))*1.3])

set(gca,'FontSize',16)
set(gcf,'color','white')

title('Sub-epidemic profile')


[quantilesc,quantilesf]=computeQuantiles(data1,curves,0);



