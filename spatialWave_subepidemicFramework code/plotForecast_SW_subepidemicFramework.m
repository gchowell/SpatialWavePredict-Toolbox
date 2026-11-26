% 
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [performanceTop, performanceEns, performanceBaseline, scoresByHorizonOut, identifiabilityOut, forecastQuantilesOut]=plotForecast_SW_subepidemicFramework(outbreakx_pass,caddate1_pass,forecastingperiod_pass,weight_type1_pass)

% Revision SW-3. Third output: the flat baseline at the full horizon,
% [MAE MSE coverage WIS]. Fourth output: the scores of every horizon, the
% table also written to performance-forecasting-byHorizon-...csv. Columns:
% model (0 = baseline, k = ranked k, 100+k = Ensemble(k)), horizon, MAE, MSE,
% coverage of the 95% PI, WIS, WIS relative to the baseline. The score of
% horizon h is the average over the first h time units of the forecast.
% Used by spatialWaveRollingForecast.m.
performanceTop=[]; performanceEns=[]; performanceBaseline=[]; scoresByHorizonOut=[];
% Revision SW-6. Fifth output: one row per ranked model, [rank, number of
% estimated parameters, number of them identified] (spatialWaveIdentifiability.m).
% The full report is written to identifiability-...csv.
identifiabilityOut=zeros(0,3); piiRanked={}; piiRanks=[];
% Revision SW-8. Sixth output: the quantile forecasts that were scored, one
% row per model and forecast step: [model, step, observation, 23 quantiles]
% (levels as in spatialWaveQuantileScores.m; model codes as in the fourth
% output). spatialWaveRollingForecast.m keeps them so that intervals can be
% widened by past forecast errors (spatialWaveModelErrorIntervals.m).
forecastQuantilesOut=zeros(0,26);
scoresByHorizon=[];

% Generate short-term forecasts - incidence

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
% <================== Load the parameter values ===============================>
% <============================================================================>

% options.m
[cumulative1_INP, outbreakx_INP, caddate1_INP, cadregion_INP, caddisease_INP, datatype_INP, DT_INP, datevecfirst1_INP, datevecend1_INP, numstartpoints_INP, topmodelsx_INP, M_INP, flag1_INP,typedecline2_INP]=options

% Revision SW-1: the error structure that goes with the estimation method,
% as in the fit driver, so the file names below match those of the fit.
dist1=spatialWaveErrorStructure(method1,dist1);

% Revision SW-2: random seed (setting rngSeed), so the same forecast call on
% the same fit files gives the same numbers.
forecastRngInfo=spatialWaveSeed('Forecast');

% options_forecast.m
[getperformance_INP, deletetempfiles_INP, forecastingperiod_INP, weight_type1_INP]=options_forecast

% <============================================================================>
% <================================ Dataset ======================================>
% <============================================================================>

if exist('outbreakx_pass','var')==1 & isempty(outbreakx_pass)==0

    outbreakx=outbreakx_pass;

else
    outbreakx=outbreakx_INP;

end

if exist('caddate1_pass','var')==1 & isempty(caddate1_pass)==0

    caddate1=caddate1_pass;
else
    caddate1=caddate1_INP;
end


cadregion=cadregion_INP;

caddisease=caddisease_INP;

datatype=datatype_INP;

datevecfirst1=datevecfirst1_INP;

datevecend1=datevecend1_INP;

cumulative1=cumulative1_INP;

DT=DT_INP; % temporal resolution in days (1=daily data, 7=weekly data).

if DT==1
    cadtemporal='daily';
elseif DT==7
    cadtemporal='weekly';
elseif DT==365
    cadtemporal='yearly';
end


cadfilename2=strcat(cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-state-',num2str(outbreakx),'-',caddate1)


% <============================================================================>
% <============================Adjustments to data ============================>
% <============================================================================>

%smoothfactor1=1; % <smoothfactor1>-day smoothing

%calibrationperiod1=90; % calibrates model using the most recent <calibrationperiod1> days  where calibrationperiod1<length(data1)

% <=============================================================================>
% <=========================== Statistical method ==============================>
% <=============================================================================>

%method1=0; % Type of estimation method: 0 = LSQ

%dist1=0;  % dist1=2 neg.binomial; dist1=0 -- >Normal distribution to model error structure

M=M_INP; %number of bootstrap realizations to generate uncertainty estimates

% <==============================================================================>
% <========================= Growth model ==========================================>
% <==============================================================================>

%npatches_fixed=4; % maximum number of subepidemics considered in epidemic trajectory fit

GGM=0;  % 0 = GGM
GLM=1;  % 1 = GLM
GRM=2;  % 2 = GRM
LM=3;   % 3 = LM
RICH=4; % 4 = Richards

flag1=flag1_INP;

typedecline2=typedecline2_INP; % 1=exponential decline in subepidemic size; 2=power-law decline in subepidemic size

% <==============================================================================>
% <================= Number of best fitting models used to generate ensemble model ===============>
% <==============================================================================>

finalRanking=loadSpatialWaveFinalRanking(strcat('./output/ABC-original-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-',num2str(onset_fixed),'-typedecline-',num2str(sum(typedecline2)),'-smoothing-',num2str(smoothfactor1),'-',cadfilename2,'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'.mat'));
% Only the selected/refitted models have finalized rank files.
topmodels1=1:min(topmodelsx_INP,finalRanking.numModels);

if npatches_fixed==1
    topmodels1=1;
end

% Revision SW-1: a grid that always holds every model (the grid built from
% the prime factors was too small for 16, 24, 32, ... models).
[rows,cols]=spatialWavePanelGrid(length(topmodels1));

% Revision SW-10: Ensemble(2) to Ensemble(n) take one panel each of figure 500,
% three per row. The grid was fixed at 1 by 3 inside getensemblesubepidemics.m,
% so the ensemble of five models asked for a fourth panel and the run stopped.
% Up to four models the grid is 1 by 3, as before.
ensemblePanelGrid=[max(1,ceil((length(topmodels1)-1)/3)) 3];

% <==============================================================================>
% <========================== Forecasting parameters ============================>
% <==============================================================================>


getperformance=getperformance_INP; % flag or indicator variable (1/0) to calculate forecasting performance or not

deletetempfiles=deletetempfiles_INP; %flag or indicator variable (1/0) to delete Forecast..mat files after use

if exist('forecastingperiod_pass','var')==1 & isempty(forecastingperiod_pass)==0
    forecastingperiod=forecastingperiod_pass; %forecast horizon (number of data points ahead)
else
    forecastingperiod=forecastingperiod_INP; %forecast horizon (number of data points ahead)
end

printscreen1=1; % print plots with the results

% <==============================================================================>
% <====================== weighting scheme for ensemble model ============================>
% <==============================================================================>

if exist('weight_type1_pass','var')==1 & isempty(weight_type1_pass)==0
    weight_type1=weight_type1_pass; % -1= equally weighted from the top models, 0=same as 1 (it was 1/AICc, which is not a weight), 1= based on relative likelihood (Akaike weights), 2=based on WISC during calibration, 3=based on WISF during forecasting performance at previous time period (week)

else
    weight_type1=weight_type1_INP; % -1= equally weighted from the top models, 0=same as 1 (it was 1/AICc, which is not a weight), 1= based on relative likelihood (Akaike weights), 2=based on WISC during calibration, 3=based on WISF during forecasting performance at previous time period (week)
end


WISC_hash=zeros(length(topmodels1),1); % vector that saves the WISC based on calibration to be used with weight_type1=2

WISF_hash=zeros(length(topmodels1),200); % vector that saves the WISF based on prior forecasting performance to be used with weight_type1=2


% <=======================================================================================>
% <========== Initialize variables to store forecast metrics across models ===================================>
% <=======================================================================================>

RMSECSS=[];
MSECSS=[];
MAECSS=[];
PICSS=[];
MISCSS=[];
RMSEFSS=[];
MSEFSS=[];
MAEFSS=[];
PIFSS=[];
MISFSS=[];

WISCSS=[];
WISFSS=[];

quantilescs=[];
quantilesfs=[];

% <===========================================================================================>
% <========== plot short-term forecasts for best models and the ensemble model ==================================>
% <===========================================================================================>

forecasts_best=[];
forecasts_ENS2=[];
forecasts_ENS3=[];
forecasts_ENS4=[];

npatches=npatches_fixed;

cc1=1;

for run_id=-1

    cc1=1;

    close all

    run_id=0;

    cadfilename2=strcat(cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-state-',num2str(outbreakx),'-',caddate1);

    AICc_rank1=[];
    relativelik_rank1=[];

    for rankx=topmodels1

        % Resetting the quantile array
        combinedQuantiles = []; 

        % <========================================================================================>
        % <================================ Load model results =========================================>
        % <========================================================================================>

        modelIdentity=load(strcat('./output/modifiedLogisticPatch-original-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-',num2str(onset_fixed),'-typedecline-',num2str(sum(typedecline2)),'-smoothing-',num2str(smoothfactor1),'-',cadfilename2,'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'-rank-',num2str(rankx),'.mat'),'fitRunId','fitResult');
        assertSpatialWaveFitIdentity(modelIdentity,rankx,finalRanking);
        fitSettings=[]; % revision SW-5: read from the fit file below (absent in fits made before SW-2)
        bootstrapDiagnostics=[]; % revision SW-6: likewise
        load (strcat('./output/modifiedLogisticPatch-original-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-',num2str(onset_fixed),'-typedecline-',num2str(sum(typedecline2)),'-smoothing-',num2str(smoothfactor1),'-',cadfilename2,'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'-rank-',num2str(rankx),'.mat'))

        % Revision SW-5: day-of-week factors of this fit, estimated again from
        % its calibration data with the fit settings saved in the fit file
        % (not those in force now). weekdayCurve is [] when the fit was made
        % without the adjustment.
        weekday=spatialWaveWeekdayFactors(data1(:,2),fitSettings);

        % Revision SW-6: identifiability report of this ranked model, from the
        % refits and search bounds saved with the fit (no fit is repeated).
        % The file is written again after every ranked model, so it holds
        % all models read so far.
        piiOne=spatialWaveSavedIdentifiability(Phatss,P0,bootstrapDiagnostics);
        if ~isempty(piiOne)
            piiRanked{end+1}=piiOne; %#ok<AGROW>
            piiRanks(end+1)=rankx; %#ok<AGROW>
            identifiabilityOut=[identifiabilityOut;[rankx piiOne.nFree piiOne.nIdentified]]; %#ok<AGROW>
            spatialWaveWritePII(strcat('./output/identifiability-',spatialWaveOutputTag(),'onsetfixed-',num2str(onset_fixed),'-typedecline-',num2str(sum(typedecline2)),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'),piiRanked,piiRanks);
            notNames=piiOne.name(strcmp(piiOne.verdict,'not identified'))';
            if isempty(notNames), notText='none'; else, notText=strjoin(notNames,', '); end
            if isempty(piiOne.tradeoffs), tradeText='none'; else, tradeText=strjoin(piiOne.tradeoffs',', '); end
            fprintf('Identifiability of ranked model %d: %d of %d estimated parameter(s) identified; not identified: %s; trade-offs: %s\n', ...
                rankx,piiOne.nIdentified,piiOne.nFree,notText,tradeText);
        end

        rankx
        AICc_rank1=[AICc_rank1;[rank1 AICc_best]];
        relativelik_rank1=[relativelik_rank1;[rank1 relativelik_i(rankx)]];

        %         if(relativelik_i(rankx))<0.05
        %             relativelik_i(rankx)
        %             strcat('drop model rank=',num2str(rankx))
        %             topmodels1=1:1:rank1-1;
        %             break
        %         end

        % Revision SW-1. The forecast simulates the model that was fitted: the
        % number of sub-epidemics it was fitted with (patchCap), which can be
        % larger than the number active so far. Beyond that the wave may
        % travel on by spatialWaveForecastExtraPatches() further
        % sub-epidemics, which can start only AFTER the calibration period,
        % so every forecast curve equals its fitted curve over the
        % calibration period. Before, "npatches+10" was applied from the
        % first day (changing the curve inside the calibration period when
        % the limit npatches_fixed was binding), and not at all to a wave
        % model with one active sub-epidemic, which could then never start
        % a second one.
        patchCap=modelIdentity.fitResult.patchCap;
        waveExtension=struct('extraPatches',spatialWaveForecastExtraPatches(), ...
            'unlockTime',data1(end,1));
        nSkippedRefits=0;
        firstSkipMessage='';
        nImmediateAtUnlock=0;

        d_hat=1;

        cc3=1;

        % <========================================================================================>
        % <================================ Compute short-term forecast ==================================>
        % <========================================================================================>

        timevect=(data1(:,1));

        timevect2=(0:t_window(end)-1+forecastingperiod);

        % vector to store forecast mean curves
        curvesforecasts1=[];

        % vector to store forecast prediction curves
        curvesforecasts2=[];

        color1=['r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';];


        if printscreen1
            figure(10)
            %subplot(1,length(topmodels1),cc1)
            subplot(rows,cols,cc1)

        end

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

            try
                weekdayCurve=spatialWaveWeekdayCurve(weekday,length(timevect2));
                if isempty(weekdayCurve)
                    I0forecast=data1(1,2);
                else
                    I0forecast=data1(1,2)/weekdayCurve(1);
                end
                [~,x,totinc,patchIncidence,activation]=simulateSpatialWave(timevect2,I0forecast,patchCap, ...
                    onset_fixed,onset_thr,flag1,typedecline1,r_hat,p_hat,a_hat,K_hat,q_hat,[],waveExtension);
            catch simulationError
                % A refit that cannot be simulated over the forecast period
                % is left out and counted (spatialWaveSkippedRefits.m).
                if ~strncmp(simulationError.identifier,'SpatialWave:',12)
                    rethrow(simulationError);
                end
                nSkippedRefits=nSkippedRefits+1;
                if isempty(firstSkipMessage)
                    firstSkipMessage=simulationError.message;
                end
                continue
            end
            nImmediateAtUnlock=nImmediateAtUnlock+double(activation.immediateAtUnlock);

            if ~isempty(weekdayCurve) % revision SW-5: the curve on the scale of the data
                totinc=totinc(:).*weekdayCurve;
                patchIncidence=bsxfun(@times,patchIncidence,weekdayCurve);
            end

            bestfit=totinc;

            for j=1:min(size(patchIncidence,2),size(color1,1))

                incidence1=patchIncidence(:,j);

                if printscreen1
                    plot(timevect2,incidence1,color1(j,:))
                end

                hold on

            end

            gray1=gray(10);

            if printscreen1
                plot(timevect2,totinc,'color',gray1(7,:))
                hold on
            end

            curvesforecasts1=[curvesforecasts1 totinc];

            forecasts2=AddPoissonError(cumsum(totinc),20,dist1,factor1,d_hat);

            curvesforecasts2=[curvesforecasts2 forecasts2];

        end

        % Revision SW-1: one policy for refits that could not be simulated.
        spatialWaveSkippedRefits('Forecast',rankx,nSkippedRefits,M,firstSkipMessage);
        nForecastCurves=size(curvesforecasts1,2);
        if nImmediateAtUnlock>0
            warning('SpatialWave:PatchLimitBinding', ...
                ['Ranked model %d: in %d of %d forecast curves a further sub-epidemic starts on the first ', ...
                 'forecast day, because the limit npatches_fixed = %d held it back during calibration. ', ...
                 'Consider fitting with a larger npatches_fixed.'], ...
                rankx,nImmediateAtUnlock,nForecastCurves,npatches_fixed);
        end

        [quantilesc,quantilesf]=computeQuantiles(data1,curvesforecasts2,forecastingperiod);

        quantilescs=[quantilescs;quantilesc];

        quantilesfs=[quantilesfs;quantilesf];

        % Exporting the quantile forecasts
        quantNamesRanked = {'Q_0.010', 'Q_0.025', 'Q_0.050', 'Q_0.100', 'Q_0.150', 'Q_0.200', 'Q_0.250', 'Q_0.300', 'Q_0.350', 'Q_0.400', 'Q_0.450', 'Q_0.500', 'Q_0.550', 'Q_0.600', 'Q_0.650', 'Q_0.700', 'Q_0.750', 'Q_0.800', 'Q_0.850', 'Q_0.900', 'Q_0.950', 'Q_0.975', 'Q_0.990'};
       
        % Creating the array with the quantile forecasts 
        combinedQuantiles = [quantilesc; quantilesf];
        combinedQuantilesTable = array2table(combinedQuantiles, 'VariableNames', quantNamesRanked);
        
        % Saving the results 
         writetable(combinedQuantilesTable,strcat('./output/quantile-ranked(', num2str(rank1),')-',spatialWaveOutputTag(),'onsetfixed-',num2str(onset_fixed),'-typedecline-',num2str(sum(typedecline2)),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-horizon-',num2str(forecastingperiod),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'))

        % compute doubling times

        % Revision SW-1: meandoublingtime starts empty. It used to start as M zeros
        % and the means were appended after them, so half of the values behind the
        % mean doubling time and its interval were zeros. Curves that never double
        % (NaN) are left out.
        % The loops run over the curves that exist (refits can be skipped).
        meandoublingtime=[];

        doublingtimess=zeros(30,nForecastCurves)+NaN;

        maxd=1;

        for j=1:nForecastCurves

            [tds,C0data,curve,doublingtimes]=getDoublingTimeCurve(max(curvesforecasts1(:,j),0),DT,0);

            doublingtimess(1:length(doublingtimes),j)=doublingtimes;

            if maxd<length(doublingtimes)
                maxd=length(doublingtimes);
            end

            meandoublingtime=[meandoublingtime;mean(doublingtimes)];
        end

        doublingtimess=doublingtimess(1:maxd,1:nForecastCurves);

        seq_doublingtimes=[];

        for j=1:maxd

            index1=find(~isnan(doublingtimess(j,:)));

            seq_doublingtimes=[seq_doublingtimes;[j mean(doublingtimess(j,index1)) quantile(doublingtimess(j,index1),0.025) quantile(doublingtimess(j,index1),0.975) length(index1)./nForecastCurves]];

        end

        seq_doublingtimes % [ith doubling, mean, 95%CI LB, 95%CI UB, prob. i_th doubling]

        % Mean doubling times
        meandoublingtime=meandoublingtime(~isnan(meandoublingtime));
        if isempty(meandoublingtime) % no curve doubled
            dmean=NaN; dLB=NaN; dUB=NaN;
        else
            dmean=mean(meandoublingtime);
            dLB=quantile(meandoublingtime,0.025);
            dUB=quantile(meandoublingtime,0.975);
        end

        param_doubling=[dmean dLB dUB]

        % <=============================================================================================>
        % <============================== Save file with doubling time estimates =======================>
        % <=============================================================================================>

        T = array2table(seq_doublingtimes);
        T.Properties.VariableNames(1:5) = {'i_th doubling','db mean','db 95%CI LB','db 95% CI UB','prob. i_th doubling'};

        writetable(T,strcat('./output/doublingTimes-ranked(', num2str(rank1),')-',spatialWaveOutputTag(),'onsetfixed-',num2str(onset_fixed),'-typedecline-',num2str(sum(typedecline2)),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-horizon-',num2str(forecastingperiod),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'))


        % <==============================================================================================>
        % <================================ Plot short-term forecast ============================================>
        % <==============================================================================================>

        if printscreen1

            title(strcat(num2ordinal(rank1),' Selected Refit'))

            line1=plot(data1(:,1),data1(:,2),'ko');
            set(line1,'LineWidth',2)


            axis([0 length(timevect2)-1 0 max(data1(:,2))*2])


            line2=[timevect(end) 0;timevect(end) max(data1(:,2))*2];


            line1=plot(line2(:,1),line2(:,2),'k--');
            set(line1,'LineWidth',2)


            datenum1=datenum([str2num(caddate1(7:10)) str2num(caddate1(1:2)) str2num(caddate1(4:5))]);

            datevec1=datevec(datenum1+forecastingperiod*DT);

            wave=[datevecfirst1 datevec1(1:3)];


            % plot dates in x axis
            'day='
            datenum1=datenum(wave(1:3))+timelags*DT; % start of fall wave (reference date)
            datestr(datenum1)

            datenumIni=datenum1;
            datenumEnd=datenum(wave(4:6))

            dates1=datestr(datenumIni:DT:datenumEnd,'mm-dd');

            if DT==1

                set(gca, 'XTick', 0:3:length(dates1(:,1))-1);
                set(gca, 'XTickLabel', strcat('\fontsize{14}',dates1(1:3:end,:)));

                xticklabel_rotate;

            elseif DT==7

                set(gca, 'XTick', 0:2:length(dates1(:,1))-1);
                set(gca, 'XTickLabel', strcat('\fontsize{14}',dates1(1:2:end,:)));

                xticklabel_rotate;

            elseif DT==365

                years1=wave(1)+timelags:wave(1)+timelags+length(dates1(:,1))-1;

                set(gca,'XTick',0:1:length(years1)-1);
                set(gca, 'XTickLabel', strcat('\fontsize{14}',num2str(years1')));

            end


            line1=plot(line2(:,1),line2(:,2),'k--');
            set(line1,'LineWidth',2)

            ylabel(strcat(caddisease,{' '},datatype))

            %title(strcat('Sub-epidemic Model Forecast',{' '},getUSstateName(outbreakx),{' '},'- Reported by',{' '},caddate1))

            set(gca,'FontSize',24)
            set(gcf,'color','white')


        end


        LB1=quantile(curvesforecasts2',0.025);
        LB1=(LB1>=0).*LB1;

        UB1=quantile(curvesforecasts2',0.975);
        UB1=(UB1>=0).*UB1;


        if rank1==1 %top model
            %store forecast for best model
            forecasts_best=[forecasts_best;[median(curvesforecasts2(end-forecastingperiod+1:end,:),2) LB1(end-forecastingperiod+1:end)' UB1(end-forecastingperiod+1:end)']];
        end



        if printscreen1

            figure(11)
            %subplot(1,length(topmodels1),cc1)
            subplot(rows,cols,cc1)


            hold on

            h=area(timevect2',[LB1' UB1'-LB1'])
            hold on

            h(1).FaceColor = [1 1 1];
            h(2).FaceColor = [0.8 0.8 0.8];

            %line1=plot(timevect2,quantile(curvesforecasts2',0.5),'r-')

            line1=plot(timevect2,median(curvesforecasts2,2),'r-');

            set(line1,'LineWidth',2)

            if  1
                line1=plot(timevect2,quantile(curvesforecasts2',0.025),'k--');
                set(line1,'LineWidth',2)

                line1=plot(timevect2,quantile(curvesforecasts2',0.975),'k--');
                set(line1,'LineWidth',2)
            end



            gray1=gray(10);

            % plot time series datalatest
            line1=plot(data1(:,1),data1(:,2),'ko');
            set(line1,'LineWidth',2)


            axis([0 length(timevect2)-1 0 max(quantile(curvesforecasts2',0.975))*1.2])

            line2=[timevect(end) 0;timevect(end) max(quantile(curvesforecasts2',0.975))*1.20];

            box on

            line1=plot(line2(:,1),line2(:,2),'k--');
            set(line1,'LineWidth',2)


            % plot dates in x axis
            'day='
            datenum1=datenum(wave(1:3))+timelags*DT; % start of fall wave (reference date)
            datestr(datenum1)

            datenumIni=datenum1;
            datenumEnd=datenum(wave(4:6))

            dates1=datestr(datenumIni:DT:datenumEnd,'mm-dd');


            if DT==1

                set(gca, 'XTick', 0:3:length(dates1(:,1))-1);
                set(gca, 'XTickLabel', strcat('\fontsize{14}',dates1(1:3:end,:)));

                xticklabel_rotate;

            elseif DT==7

                set(gca, 'XTick', 0:2:length(dates1(:,1))-1);
                set(gca, 'XTickLabel', strcat('\fontsize{14}',dates1(1:2:end,:)));

                xticklabel_rotate;

            elseif DT==365

                years1=wave(1)+timelags:wave(1)+timelags+length(dates1(:,1))-1;

                set(gca,'XTick',0:1:length(years1)-1);
                set(gca, 'XTickLabel', strcat('\fontsize{14}',num2str(years1')));

            end

            line1=plot(line2(:,1),line2(:,2),'k--');
            set(line1,'LineWidth',2)

            ylabel(strcat(caddisease,{' '},datatype))

            title(strcat(num2ordinal(rank1),' Selected Refit'))

            %title(strcat('Sub-epidemic Model Forecast-',{' '},getUSstateName(outbreakx),{' '},'- Reported by',{' '},caddate1))

            set(gca,'FontSize',24)
            set(gcf,'color','white')

        end



        % <=========================================================================================>
        % <================================ Save short-term forecast results ==================================>
        % <=========================================================================================>

        save(strcat('./output/Forecast-modifiedLogisticPatch-original-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-',num2str(onset_fixed),'-typedecline-',num2str(sum(typedecline2)),'-smoothing-',num2str(smoothfactor1),'-',cadfilename2,'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'-forecastingperiod-',num2str(forecastingperiod),'-rank-',num2str(rankx),'.mat'),'curvesforecasts1','curvesforecasts2','datevecfirst1','datevecend1','timevect','timevect2','timelags','cadtemporal','data1','fitRunId','fitResult')


        % <=============================================================================================>
        % <=================== Plot data for the forecast period (if getperformance=1) ==================================>
        % <=============================================================================================>

        % Use the same verification-date rule as the ensemble path.
        datenum1=getSpatialWaveVerificationDate(caddate1,DT);

        if getperformance && forecastingperiod>0

            data2=getData(cumulative1,cadtemporal,caddisease,datatype,cadregion,DT,datevecfirst1,datevecend1,datevec(datenum1),outbreak1,forecastingperiod);

            timevect2=(data1(end,1)+1:(data1(end,1)+1+forecastingperiod-1));

            if printscreen1

                line2=plot(timevect2,data2,'ro')
                set(line2,'LineWidth',2)

            end


            % <=============================================================================================>
            % <============================== Save file with forecast ======================================>
            % <=============================================================================================>

            forecastdata=[str2num(datestr((datenumIni:DT:datenumEnd)','yyyy')) str2num(datestr((datenumIni:DT:datenumEnd)','mm')) str2num(datestr((datenumIni:DT:datenumEnd)','dd')) [data1(:,2);data2] median(curvesforecasts2,2) LB1' UB1'];

            T = array2table(forecastdata);
            T.Properties.VariableNames(1:7) = {'year','month','day','data','median','LB','UB'};
            writetable(T,strcat('./output/ranked(', num2str(rank1),')-',spatialWaveOutputTag(),'onsetfixed-',num2str(onset_fixed),'-typedecline-',num2str(sum(typedecline2)),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-horizon-',num2str(forecastingperiod),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'))


        else

            timevect2=[];

            data2=[];

            % <=============================================================================================>
            % <============================== Save file with forecast ======================================>
            % <=============================================================================================>

            forecastdata=[str2num(datestr((datenumIni:DT:datenumEnd)','yyyy')) str2num(datestr((datenumIni:DT:datenumEnd)','mm')) str2num(datestr((datenumIni:DT:datenumEnd)','dd')) [data1(:,2);zeros(forecastingperiod,1)+NaN] median(curvesforecasts2,2) LB1' UB1'];

            T = array2table(forecastdata);
            T.Properties.VariableNames(1:7) = {'year','month','day','data','median','LB','UB'};
            writetable(T,strcat('./output/ranked(', num2str(rank1),')-',spatialWaveOutputTag(),'onsetfixed-',num2str(onset_fixed),'-typedecline-',num2str(sum(typedecline2)),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-horizon-',num2str(forecastingperiod),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'))


        end

        datalatest2=[data1;[timevect2' data2]];

        % <==================================================================================================>
        % <========== Get forecast performance metrics for the model (if getperformance=1) =====================================>
        % <==================================================================================================>

        if getperformance

            [RMSECS_model1 MSECS_model1 MAECS_model1  PICS_model1 MISCS_model1 RMSEFS_model1 MSEFS_model1 MAEFS_model1 PIFS_model1 MISFS_model1]=computeforecastperformance(data1,datalatest2,curvesforecasts1,curvesforecasts2,forecastingperiod);

            [WISC,WISFS,~,quantilesfModel]=computeWIS(data1,datalatest2,curvesforecasts2,forecastingperiod);
            if forecastingperiod>0 % revision SW-8
                forecastQuantilesOut=[forecastQuantilesOut;[rankx+zeros(forecastingperiod,1) (1:forecastingperiod)' datalatest2(size(data1,1)+(1:forecastingperiod)',2) quantilesfModel]]; %#ok<AGROW>
            end

            WISC_hash(rankx,1)=WISC;  %save WISC for the model to be used to weight the models in the ensemble in function getensemblesubepidemics

            WISF_hash(rankx,run_id+1)=WISFS(end,end);  %save WISF for the model to be used to weight the models in the ensemble in function getensemblesubepidemics

            % Revision SW-3: the forecast scores of every horizon
            if forecastingperiod>0
                scoresByHorizon=[scoresByHorizon;[rankx+zeros(forecastingperiod,1) MAEFS_model1(:,1) MAEFS_model1(:,2) MSEFS_model1(:,2) PIFS_model1(:,2) WISFS(:,2)]]; %#ok<AGROW>
            end


            % store metrics for calibration
            RMSECSS=[RMSECSS;[rankx outbreakx datenum(caddate1) RMSECS_model1(end,end)]];
            MSECSS=[MSECSS;[rankx outbreakx datenum(caddate1) MSECS_model1(end,end)]];
            MAECSS=[MAECSS;[rankx outbreakx datenum(caddate1) MAECS_model1(end,end)]];
            PICSS=[PICSS;[rankx outbreakx datenum(caddate1) PICS_model1(end,end)]];
            MISCSS=[MISCSS;[rankx outbreakx datenum(caddate1) MISCS_model1(end,end)]];

            WISCSS=[WISCSS;[rankx outbreakx datenum(caddate1) WISC(end,end)]];

            % store metrics for short-term forecasts
            if forecastingperiod>0
                RMSEFSS=[RMSEFSS;[rankx outbreakx datenum(caddate1) RMSEFS_model1(end,end)]];
                MSEFSS=[MSEFSS;[rankx outbreakx datenum(caddate1) MSEFS_model1(end,end)]];
                MAEFSS=[MAEFSS;[rankx outbreakx datenum(caddate1) MAEFS_model1(end,end)]];
                PIFSS=[PIFSS;[rankx outbreakx datenum(caddate1) PIFS_model1(end,end)]];
                MISFSS=[MISFSS;[rankx outbreakx datenum(caddate1) MISFS_model1(end,end)]];

                WISFSS=[WISFSS;[rankx outbreakx datenum(caddate1) WISFS(end,end)]];

            end

        end

        cc1=cc1+1;

        %end

        % <=======================================================================================>
        % <=================== Get performance metrics for ensemble model(Ensemble(K)) ========================>
        % <=======================================================================================>

        if rankx>1

            %[RMSECS_model1 MSECS_model1 MAECS_model1  PICS_model1 MISCS_model1 WISC RMSEFS_model1 MSEFS_model1 MAEFS_model1 PIFS_model1 MISFS_model1 WISFS forecast1]=getensemblesubepidemics(cadfilename2,datevecfirst1,npatches_fixed,onset_fixed,smoothfactor1,outbreakx,caddate1,flag1,method1,dist1,calibrationperiod1,1:rankx,forecastingperiod,getperformance,printscreen1);

            if run_id==0
                [RMSECS_model1 MSECS_model1 MAECS_model1  PICS_model1 MISCS_model1 WISC RMSEFS_model1 MSEFS_model1 MAEFS_model1 PIFS_model1 MISFS_model1 WISFS forecast1 quantilesc quantilesf]=getensemblesubepidemics(cumulative1,cadfilename2,datevecfirst1,npatches_fixed,onset_fixed,typedecline2,smoothfactor1,outbreakx,cadregion,caddate1,caddisease,datatype,flag1,method1,dist1,calibrationperiod1,1:rankx,forecastingperiod,getperformance,weight_type1,WISC_hash,WISC_hash,printscreen1,ensemblePanelGrid);
            else
                [RMSECS_model1 MSECS_model1 MAECS_model1  PICS_model1 MISCS_model1 WISC RMSEFS_model1 MSEFS_model1 MAEFS_model1 PIFS_model1 MISFS_model1 WISFS forecast1 quantilesc quantilesf]=getensemblesubepidemics(cumulative1,cadfilename2,datevecfirst1,npatches_fixed,onset_fixed,typedecline2,smoothfactor1,outbreakx,cadregion,caddate1,caddisease,datatype, flag1,method1,dist1,calibrationperiod1,1:rankx,forecastingperiod,getperformance,weight_type1,WISC_hash,WISF_hash(:,run_id),printscreen1,ensemblePanelGrid);
            end

            quantilescs=[quantilescs;quantilesc];

            quantilesfs=[quantilesfs;quantilesf];


            %store calibration performance metrics
            RMSECSS=[RMSECSS;[100+rankx outbreakx datenum(caddate1) RMSECS_model1(end,end)]];
            MSECSS=[MSECSS;[100+rankx outbreakx datenum(caddate1) MSECS_model1(end,end)]];
            MAECSS=[MAECSS;[100+rankx outbreakx datenum(caddate1) MAECS_model1(end,end)]];
            PICSS=[PICSS;[100+rankx outbreakx datenum(caddate1) PICS_model1(end,end)]];
            MISCSS=[MISCSS;[100+rankx outbreakx datenum(caddate1) MISCS_model1(end,end)]];

            WISCSS=[WISCSS;[100+rankx outbreakx datenum(caddate1) WISC(end,end)]];

            if forecastingperiod>0

                %store metrics for short-term forecasts
                RMSEFSS=[RMSEFSS;[100+rankx outbreakx datenum(caddate1) RMSEFS_model1(end,end)]];
                MSEFSS=[MSEFSS;[100+rankx outbreakx datenum(caddate1) MSEFS_model1(end,end)]];
                MAEFSS=[MAEFSS;[100+rankx outbreakx datenum(caddate1) MAEFS_model1(end,end)]];
                PIFSS=[PIFSS;[100+rankx outbreakx datenum(caddate1) PIFS_model1(end,end)]];
                MISFSS=[MISFSS;[100+rankx outbreakx datenum(caddate1) MISFS_model1(end,end)]];

                WISFSS=[WISFSS;[100+rankx outbreakx datenum(caddate1) WISFS(end,end)]];
            end

            % Revision SW-3: the forecast scores of every horizon
            if getperformance && forecastingperiod>0 && size(MAEFS_model1,1)==forecastingperiod
                scoresByHorizon=[scoresByHorizon;[100+rankx+zeros(forecastingperiod,1) MAEFS_model1(:,1) MAEFS_model1(:,2) MSEFS_model1(:,2) PIFS_model1(:,2) WISFS(:,2)]]; %#ok<AGROW>
                if isequal(size(quantilesf),[forecastingperiod 23]) % revision SW-8
                    forecastQuantilesOut=[forecastQuantilesOut;[100+rankx+zeros(forecastingperiod,1) (1:forecastingperiod)' datalatest2(size(data1,1)+(1:forecastingperiod)',2) quantilesf]]; %#ok<AGROW>
                end
            end

            % store the first 3 ensemble forecasts including the 95% Prediction interval.
            switch rankx
                case 2
                    forecasts_ENS2=[forecasts_ENS2;[forecast1(end-forecastingperiod+1:end,:)]];

                case 3
                    forecasts_ENS3=[forecasts_ENS3;[forecast1(end-forecastingperiod+1:end,:)]];

                case 4
                    forecasts_ENS4=[forecasts_ENS4;[forecast1(end-forecastingperiod+1:end,:)]];

            end

        end

    end

    if deletetempfiles %flag or indicator variable (1/0) to delete Forecast..mat files after use

        for j=topmodels1

            delete(strcat('./output/Forecast-modifiedLogisticPatch-original-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-',num2str(onset_fixed),'-typedecline-',num2str(sum(typedecline2)),'-smoothing-',num2str(smoothfactor1),'-',cadfilename2,'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'-forecastingperiod-',num2str(forecastingperiod),'-rank-',num2str(j),'.mat'))
        end
    end

end

if getperformance

    % <============================================================================>
    % <=================plot forecasting performance metrics for the top-ranked models ==============>
    % <============================================================================>

    index1=find(MAEFSS(:,1)>=100);

    index2=setdiff(1:length(MAEFSS(:,1)),index1);

    figure(400)

    subplot(2,2,1)
    line1=plot(MAEFSS(index2,1),MAEFSS(index2,4),'k-o');
    set(line1,'linewidth',2)
    xlabel('Rank among selected refits')
    ylabel('MAE')

    set(gca,'FontSize', 16);
    set(gcf,'color','white')

    subplot(2,2,2)
    line1=plot(MSEFSS(index2,1),MSEFSS(index2,4),'k-o');
    set(line1,'linewidth',2)
    xlabel('Rank among selected refits')
    ylabel('MSE')

    set(gca,'FontSize', 16);
    set(gcf,'color','white')

    subplot(2,2,3)
    line1=plot(PIFSS(index2,1),PIFSS(index2,4),'k-o');
    set(line1,'linewidth',2)
    xlabel('Rank among selected refits')
    ylabel('Coverage of the 95% PI')

    set(gca,'FontSize', 16);
    set(gcf,'color','white')

    subplot(2,2,4)

    line1=plot(WISFSS(index2,1),WISFSS(index2,4),'k-o');
    set(line1,'linewidth',2)
    xlabel('Rank among selected refits')
    ylabel('WIS')

    set(gca,'FontSize', 16);
    set(gcf,'color','white')

    % <=============================================================================================>
    % <========================== Save file with top-ranked models' forecasting performance metrics ===================>
    % <=============================================================================================>

    performanceTop=[topmodels1' MAEFSS(index2,4) MSEFSS(index2,4) PIFSS(index2,4) WISFSS(index2,4) AICc_rank1(:,2) relativelik_rank1(:,2)];

    T = array2table(performanceTop);
    T.Properties.VariableNames(1:7) = {'i_th-ranked model','MAE','MSE','Coverage 95%PI','WIS','AICc','RelativeLikelihood'};
    writetable(T,strcat('./output/performance-forecasting-topRanked-',spatialWaveOutputTag(),'onsetfixed-',num2str(onset_fixed),'-typedecline-',num2str(sum(typedecline2)),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-horizon-',num2str(forecastingperiod),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'))

    % <============================================================================>
    % <=================plot forecasting performance metrics of the ensemble models ==============>
    % <============================================================================>

    if isempty(index1)~=1

        figure(401)

        subplot(2,2,1)
        line1=plot(1+(1:length(index1)),MAEFSS(index1,4),'k-o');
        set(line1,'linewidth',2)
        xlabel('Ensemble(i) model')
        ylabel('MAE')

        set(gca,'FontSize', 16);
        set(gcf,'color','white')

        subplot(2,2,2)
        line1=plot(1+(1:length(index1)),MSEFSS(index1,4),'k-o');
        set(line1,'linewidth',2)
        xlabel('Ensemble(i) model')
        ylabel('MSE')

        set(gca,'FontSize', 16);
        set(gcf,'color','white')

        subplot(2,2,3)
        line1=plot(1+(1:length(index1)),PIFSS(index1,4),'k-o');
        set(line1,'linewidth',2)
        xlabel('Ensemble(i) model')
        ylabel('Coverage of the 95% PI')

        set(gca,'FontSize', 16);
        set(gcf,'color','white')

        subplot(2,2,4)

        line1=plot(1+(1:length(index1)),WISFSS(index1,4),'k-o');
        set(line1,'linewidth',2)
        xlabel('Ensemble(i) model')
        ylabel('WIS')

        set(gca,'FontSize', 16);
        set(gcf,'color','white')

    end

    % <=============================================================================================>
    % <=========================== Save file with ensemble forecasting performance metrics =========================>
    % <=============================================================================================>

    performanceEns=[(1+(1:length(index1)))' MAEFSS(index1,4) MSEFSS(index1,4) PIFSS(index1,4) WISFSS(index1,4)];

    T = array2table(performanceEns);
    T.Properties.VariableNames(1:5) = {'Ensemble(i) model','MAE','MSE','Coverage 95%PI','WIS'};
    writetable(T,strcat('./output/performance-forecasting-Ensemble-',spatialWaveOutputTag(),'onsetfixed-',num2str(onset_fixed),'-typedecline-',num2str(sum(typedecline2)),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-horizon-',num2str(forecastingperiod),'-weight_type-',num2str(weight_type1),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'))

    % <=============================================================================================>
    % <======= Revision SW-3: flat baseline, and the forecast scores of every horizon ================>
    % <=============================================================================================>
    % The baseline (spatialWaveBaselineForecast.m) repeats the mean of the last
    % observations, with uncertainty from the errors the same rule made earlier
    % in the calibration data. It is scored with the functions that score the
    % models, so its rows can be compared with theirs. A model is worth using at
    % a horizon when its WIS there is below the baseline's.
    if forecastingperiod>0

        baselineSettings=spatialWaveFitSettings();
        baselineWindow=baselineSettings.baselineWindow;
        if isempty(baselineWindow)
            if DT==1
                baselineWindow=7; % daily data: one week, to average the weekly reporting cycle
            else
                baselineWindow=1; % weekly or yearly data: the last observation
            end
        end

        baseline=spatialWaveBaselineForecast(data1(:,2),forecastingperiod,baselineWindow);

        % curves in the layout of the model curves: calibration rows hold the data
        % (they are not scored), forecast rows hold the baseline
        baselineCurves1=[data1(:,2);baseline.point];
        baselineCurves2=[repmat(data1(:,2),1,size(baseline.samples,2));baseline.samples];

        [~,~,~,~,~,~,MSEFS_base,MAEFS_base,PIFS_base]=computeforecastperformance(data1,datalatest2,baselineCurves1,baselineCurves2,forecastingperiod);

        if all(isfinite(baseline.samples(:)))
            [~,WISFS_base,~,quantilesfBase]=computeWIS(data1,datalatest2,baselineCurves2,forecastingperiod);
            forecastQuantilesOut=[[zeros(forecastingperiod,1) (1:forecastingperiod)' datalatest2(size(data1,1)+(1:forecastingperiod)',2) quantilesfBase];forecastQuantilesOut];
        else
            % too few past errors for an interval: point scores only
            WISFS_base=[(1:forecastingperiod)' NaN(forecastingperiod,1)];
            PIFS_base(:,2)=NaN;
            warning('SpatialWave:BaselineWithoutInterval', ...
                ['The calibration series is too short to give the flat baseline an interval at every horizon ', ...
                 '(%d observations, window %d, horizon %d). Its MAE and MSE are reported; its WIS and coverage are not.'], ...
                size(data1,1),baseline.window,forecastingperiod);
        end

        % model code 0 = baseline; then the ranked models, then the ensembles (100 + k)
        scoresByHorizon=[[zeros(forecastingperiod,1) MAEFS_base(:,1) MAEFS_base(:,2) MSEFS_base(:,2) PIFS_base(:,2) WISFS_base(:,2)];scoresByHorizon];
        scoresByHorizon=sortrows(scoresByHorizon,[1 2]);

        % WIS of each row relative to the baseline at the same horizon (below 1 = better than the baseline)
        relativeWIS=scoresByHorizon(:,6)./WISFS_base(scoresByHorizon(:,2),2);

        T = array2table([scoresByHorizon relativeWIS]);
        T.Properties.VariableNames(1:7) = {'model','horizon','MAE','MSE','Coverage 95%PI','WIS','WIS relative to baseline'};
        writetable(T,strcat('./output/performance-forecasting-byHorizon-',spatialWaveOutputTag(),'onsetfixed-',num2str(onset_fixed),'-typedecline-',num2str(sum(typedecline2)),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-horizon-',num2str(forecastingperiod),'-weight_type-',num2str(weight_type1),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'))

        % fourth output of the function: the table just written
        scoresByHorizonOut=[scoresByHorizon relativeWIS];

        % third output of the function: the baseline at the full horizon
        performanceBaseline=[MAEFS_base(end,2) MSEFS_base(end,2) PIFS_base(end,2) WISFS_base(end,2)];

        % the full horizon on the screen: baseline first, then the models
        atEnd=scoresByHorizon(:,2)==forecastingperiod;
        rowsEnd=[scoresByHorizon(atEnd,:) relativeWIS(atEnd)];
        fprintf('\nForecast scores over %d time unit(s). Baseline: flat at %.6g, the mean of the last %d observation(s).\n', ...
            forecastingperiod,baseline.level,baseline.window);
        fprintf('  %-14s %10s %10s %10s %14s\n','model','MAE','coverage','WIS','WIS/baseline');
        for rowEnd=1:size(rowsEnd,1)
            if rowsEnd(rowEnd,1)==0
                modelLabel='baseline';
            elseif rowsEnd(rowEnd,1)>100
                modelLabel=sprintf('Ensemble(%d)',rowsEnd(rowEnd,1)-100);
            else
                modelLabel=sprintf('ranked %d',rowsEnd(rowEnd,1));
            end
            fprintf('  %-14s %10.2f %10.1f %10.2f %14.2f\n',modelLabel,rowsEnd(rowEnd,3),rowsEnd(rowEnd,5),rowsEnd(rowEnd,6),rowsEnd(rowEnd,7));
        end
        fprintf('  (scores of every horizon from 1 to %d are in the file performance-forecasting-byHorizon-... in ./output)\n',forecastingperiod);

    end

end
