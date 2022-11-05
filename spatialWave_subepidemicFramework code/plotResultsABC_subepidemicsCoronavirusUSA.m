
%% Plot model fits for the best fitting models

clear
clear global

close all

% <============================================================================>
% <=================== Declare global variables ===============================>
% <============================================================================>

global invasions
global timeinvasions
global Cinvasions
global npatches_fixed
global onset_fixed

global method1 dist1 factor1

global smoothfactor1

% <============================================================================>
% <================================ Dataset ====================================>
% <============================================================================>

outbreakx=52;
caddate1='05-11-20';

cadregion='USA';

caddisease='coronavirus';

datatype='cases';

DT=1; % temporal resolution in days (1=daily data, 7=weekly data).

if DT==1
    cadtemporal='daily';
elseif DT==7
    cadtemporal='weekly';
end

cadfilename1=strcat(cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-state-',num2str(outbreakx),'-',caddate1);


% <============================================================================>
% <============================Adjustments to data ============================>
% <============================================================================>

smoothfactor1=3; % <smoothfactor1>-day smoothing

calibrationperiod1=60; % calibrates model using the most recent <calibrationperiod1> days  where calibrationperiod1<length(data1)

% <=============================================================================>
% <=========================== Statistical method ==============================>
% <=============================================================================>

method1=0;

dist1=2;

% <==============================================================================>
% <========================= Mathematical model =================================>
% <==============================================================================>

npatches_fixed=3;

flag1=1;

% <==============================================================================>
% <======== Number of best fitting models used to generate ensemble model =======>
% <==============================================================================>

topmodels1=1:4;


caddatex=[2020 04 20];


AICc_bests=[];

for run_id=-1
%for run_id=0:1:58
    
    cc1=1;
    
    close all
    
    %i=(run_id)*30+1;
    %ARIMA_mean1=ARIMAforecasts(i:1:i+forecastingperiod-1,1);
    %ARIMA_lb1=ARIMAforecasts(i:1:i+forecastingperiod-1,11);
    %ARIMA_ub1=ARIMAforecasts(i:1:i+forecastingperiod-1,end-1);
    
    
    run_id
    
    if run_id==-1
        outbreakx=52;
        
        %caddate1='04-20-20';
        %caddate1='06-01-20';
        %caddate1='05-11-20';
        %caddate1='03-22-21';
        
        cadfilename2=strcat(cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-state-',num2str(outbreakx),'-',caddate1);
        
    else
        
        %
        
        datenum1=datenum(caddatex);
        
        %state_id = rem(run_id,52)+1;
        state_id = 52;
        
        nm=fix(run_id/52);
        %date_id = datetime(caddatex) + caldays(nm);
        date_id = datetime(caddatex) + run_id*7;
        
        date=datestr(date_id,'mm-dd-yy');
        
        outbreakx=state_id;
        
        caddate1=date;
        
        cadfilename2=strcat(cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-state-',num2str(outbreakx),'-',caddate1);
        
        
        %
    end
    
    
    load (strcat('./output/ABC-original-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-0-smoothing-',num2str(smoothfactor1),'-',cadfilename2,'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'.mat'))
    
    % remove repeated rows
    [RMSES,index1]=unique(RMSES,'rows','stable');
    PS=PS(index1,:);
    
    
    figure(100)
    
    color1=['r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';];
    
    
    for index1=topmodels1
        
        
        subplot(2,2,index1)
        
        npatches=RMSES(index1,1);
        
        onset_thr=RMSES(index1,2);
        
        typedecline1=RMSES(index1,3);
        
        
        AICc_best=RMSES(index1,4);
        
        
        % <=============================================================================================>
        % <============================ Get the best fit results =======================================>
        % <=============================================================================================>
        
        
        P=PS(index1,:);
        
        r_hat=PS(index1,1)
        p_hat=PS(index1,2)
        a_hat=PS(index1,3)
        K_hat=PS(index1,4)
        q_hat=PS(index1,5)
        
        alpha_hat=PS(index1,end-1);
        d_hat=PS(index1,end);
        
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
        
        
        IC=zeros(npatches,1);
        
        IC(1,1)=I0;
        IC(2:end,1)=1;
        
        
        invasions=zeros(npatches,1);
        timeinvasions=zeros(npatches,1);
        Cinvasions=zeros(npatches,1);
        
        invasions(1)=1;
        timeinvasions(1)=0;
        Cinvasions(1)=0;
        
        
        timevect
        
        npatches
        
        onset_thr
        
        flag1
        
        typedecline1
        
        
        [~,x]=ode15s(@modifiedLogisticGrowthPatch,timevect,IC,[],r_hat,p_hat,a_hat,K_hat,npatches,onset_thr,q_hat,flag1,typedecline1);
        %x=ode5(@modifiedLogisticGrowthPatch,timevect,IC,r_hat,p_hat,a_hat,K_hat,npatches,onset_thr,q_hat,flag1);
        
        
        invasions
        
        %     if sum(invasions)<npatches
        %
        %         invasions
        %
        %
        %         npatches=sum(invasions);
        %
        %         P=[rs_hat(1:npatches) ps_hat(1:npatches) as_hat(1:npatches) Ks_hat(1:npatches) alpha_hat d_hat];
        %
        %         PS(index1,:)=0;
        %
        %         PS(index1,1:length(P))=P;
        %
        %         RMSES(index1,1)=npatches;
        %
        %     end
        %
        
        for j=1:npatches
            
            incidence1=[x(1,j);diff(x(:,j))];
            
            line1=plot(timevect,incidence1,color1(j,:))
            
            set(line1,'Linewidth',2)
            hold on
            
        end
        
        y=sum(x,2);
        
        totinc=[y(1,1);diff(y(:,1))];
        
        if onset_thr>0
            totinc(1)=totinc(1)-(npatches-1);
        end
        
        bestfit=totinc;
        
        
        hold on
        line1=plot(timevect,data,'bo')
        set(line1,'Linewidth',2,'markerSize',8)
        
        
        line1=plot(timevect,totinc,'k')
        set(line1,'Linewidth',3)
        
        
        ylabel('Cases')
        
        set(gca,'FontSize',24)
        set(gcf,'color','white')
        
        
        title(strcat(num2ordinal(index1),{' '},' Ranked Model; AICc=',num2str(AICc_best,6)))
        
        legend(strcat('Sub-epidemics=',num2str(npatches),'; C_{thr}=',num2str(onset_thr)))
        
        %legend(strcat('Num. Sub-epidemics=',num2str(npatches)))
        
        
        [npatches onset_thr]
        
        
        if (method1==0 & dist1==2)  % calculate the overdispersion factor
            
            %     [coef,ns]=getMeanVarLinear(data,totinc,6);
            %
            %     if coef>0
            %         factor1=coef;
            %     else
            
            
            % estimate dispersion in data
            binsize1=7;
            
            %[ratios,~]=getMeanVarianceRatio(smooth(data,smoothfactor1),binsize1,2);
            [ratios,~]=getMeanVarianceRatio(data,binsize1,2);
            
            index1=find(ratios(:,1)>0);
            
            factor1=mean(ratios(index1,1));
            
            factor1
            
        end
        
        
        'ABC estimates:'
        '[npatches onset_thr q]'
        
        [npatches onset_thr]
        
        factor1
        
        
        
    end
end


subplot(2,2,3)
xlabel('Time (days)');

subplot(2,2,4)
xlabel('Time (days)');

save(strcat('./output/AICc_bests-modifiedLogisticPatch-original-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-0-smoothing-',num2str(smoothfactor1),'-',caddisease,'-',datatype,'-',cadregion,'-state-',num2str(outbreakx),'-dateini-',datestr(datenum(caddatex),'mm-dd-yy'),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'-topmodels-',num2str(topmodels1(end)),'.mat'),'AICc_bests','-mat')
