function summary = plotReproductionNumber_SW(outbreakx_pass,caddate1_pass,forecastingperiod_pass)
%PLOTREPRODUCTIONNUMBER_SW Effective reproduction number of the ranked models.
%   plotReproductionNumber_SW
%   summary = plotReproductionNumber_SW(outbreakx,caddate1,forecastingperiod)
%
%   Reads the fit that Run_SW_subepidemicFramework.m saved for the settings
%   in options.m and computes, for every ranked model, the effective
%   reproduction number Rt over the calibration period and
%   forecastingperiod time units beyond it, with a 95% interval from the
%   bootstrap refits (spatialWaveReproductionNumber.m). The generation
%   interval is set in options_Rt_sw.m.
%
%   outbreakx, caddate1   area and date ('mm-dd-yyyy'); [] or omitted: options.m
%   forecastingperiod     [] or omitted: options_forecast.m
%
%   Written to ./output:
%     Rt-ranked(k)-...csv   one row per time unit: time, whether it is in
%                           the calibration period, Rt of the ranked fit,
%                           median and 95% interval over the refits, and the
%                           share of refits with Rt above 1
%     Rt-summary-...csv     one row per ranked model: Rt at the end of the
%                           calibration period and at the end of the forecast
%   and one figure with a panel per ranked model.
%
%   summary  one row per ranked model:
%     [rank, time of the calibration end, estimate, median, LB, UB, share above 1,
%      last time, estimate, median, LB, UB, share above 1]
%
%   No fit is repeated. Rt comes from the model curve, so it is only as good
%   as the fit, and beyond the calibration period only as good as the
%   forecast; the interval reflects parameter uncertainty given the model.
%
%   Revision SW-7 (after plotReproductionNumber.m of the SubEpiPredict toolbox).
%
%   See also spatialWaveReproductionNumber, spatialWaveRt, options_Rt_sw

global method1 npatches_fixed onset_fixed dist1 smoothfactor1 calibrationperiod1 %#ok<GVMIS>

[~,outbreakx_INP,caddate1_INP,cadregion,caddisease,datatype,DT,~,~,~,topmodelsx,~,flag1,typedecline2]=options;
[~,~,forecastingperiod_INP]=options_forecast;
[typeGI,meanGI,varGI,minCoverage]=options_Rt_sw;
rt=struct('type',typeGI,'mean',meanGI,'variance',varGI,'minCoverage',minCoverage);

if nargin>=1 && ~isempty(outbreakx_pass), outbreakx=outbreakx_pass; else, outbreakx=outbreakx_INP; end
if nargin>=2 && ~isempty(caddate1_pass), caddate1=caddate1_pass; else, caddate1=caddate1_INP; end
if nargin>=3 && ~isempty(forecastingperiod_pass), forecastingperiod=forecastingperiod_pass; else, forecastingperiod=forecastingperiod_INP; end

switch DT
    case 1
        cadtemporal='daily';
    case 7
        cadtemporal='weekly';
    case 365
        cadtemporal='yearly';
    otherwise
        error('SpatialWave:InvalidRtSettings','DT = %s in options.m; the toolbox reads daily (1), weekly (7) and yearly (365) data.',mat2str(DT));
end
distFit=spatialWaveErrorStructure(method1,dist1);
fitStem=strcat('./output/modifiedLogisticPatch-original-npatchesfixed-',num2str(npatches_fixed), ...
    '-onsetfixed-',num2str(onset_fixed),'-typedecline-',num2str(sum(typedecline2)),'-smoothing-',num2str(smoothfactor1),'-', ...
    cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-state-',num2str(outbreakx),'-',caddate1, ...
    '-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(distFit), ...
    '-calibrationperiod-',num2str(calibrationperiod1),'-rank-');
if exist(strcat(fitStem,'1.mat'),'file')~=2
    error('SpatialWave:RtFitNotFound', ...
        ['No fit was found for the settings in options.m and the date %s (file %s). ', ...
         'Run Run_SW_subepidemicFramework first, or check npatches_fixed, typedecline2, smoothfactor1, method1 and calibrationperiod1.'], ...
        caddate1,strcat(fitStem,'1.mat'));
end
fileTag=strcat(spatialWaveOutputTag(npatches_fixed,smoothfactor1,calibrationperiod1),'onsetfixed-',num2str(onset_fixed),'-typedecline-',num2str(sum(typedecline2)),'-flag1-',num2str(flag1(1)), ...
    '-method-',num2str(method1),'-dist-',num2str(distFit),'-horizon-',num2str(forecastingperiod),'-', ...
    cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1);

fprintf('\nEffective reproduction number, area %d, data to %s, %d time unit(s) beyond the calibration period.\n',outbreakx,caddate1,forecastingperiod);
generationNames={'gamma','exponential','fixed'};
if rt.type==1
    fprintf('Generation interval (options_Rt_sw.m): %s, mean %g, variance %g, in the time unit of the data.\n',generationNames{rt.type},rt.mean,rt.variance);
else
    fprintf('Generation interval (options_Rt_sw.m): %s, mean %g, in the time unit of the data.\n',generationNames{rt.type},rt.mean);
end

summary=zeros(0,13);
results={};
for rank1=1:topmodelsx
    fitName=strcat(fitStem,num2str(rank1),'.mat');
    if exist(fitName,'file')~=2
        break   % fewer ranked models were saved than topmodelsx
    end
    S=load(fitName);
    out=spatialWaveReproductionNumber(S,forecastingperiod,rt);
    results{end+1}=out; %#ok<AGROW>
    summary=[summary;[rank1 out.atCalibrationEnd out.atEnd]]; %#ok<AGROW>

    T=array2table([out.time double(out.calibration) out.estimate out.median out.LB out.UB out.aboveOne]);
    T.Properties.VariableNames(1:7)={'time','calibration period (1/0)','Rt of the ranked fit','Rt median','Rt 95% LB','Rt 95% UB','share of refits with Rt above 1'};
    writetable(T,strcat('./output/Rt-ranked(',num2str(rank1),')-',fileTag,'.csv'));

    if isnan(out.firstTime)
        warning('SpatialWave:RtNotEstimable', ...
            ['Ranked model %d: the curve is shorter than the generation interval needs (coverage below %g everywhere), ', ...
             'so no Rt value is reported. Check the units and values in options_Rt_sw.m.'],rank1,rt.minCoverage);
    else
        fprintf('Ranked model %d: Rt from time %g on (earlier values are blank: the curve does not yet cover %g of the generation interval).\n', ...
            rank1,out.firstTime,rt.minCoverage);
    end
    if out.nSkipped>0
        warning('SpatialWave:RtSkippedRefits','Ranked model %d: %d of %d refits could not be simulated and are left out. First reason: %s', ...
            rank1,out.nSkipped,out.nSkipped+out.nRefits,out.firstSkipMessage);
    end
end

T=array2table(summary);
T.Properties.VariableNames(1:13)={'ranked model','calibration end (time)','Rt fit','Rt median','Rt 95% LB','Rt 95% UB','share above 1', ...
    'last time','Rt fit (last)','Rt median (last)','Rt 95% LB (last)','Rt 95% UB (last)','share above 1 (last)'};
writetable(T,strcat('./output/Rt-summary-',fileTag,'.csv'));

fprintf('\nRt at the end of the calibration period (time %g):\n',summary(1,2));
fprintf('  %-8s %10s %10s %22s %16s\n','model','fit','median','95% interval','share above 1');
for k=1:size(summary,1)
    fprintf('  %-8s %10.3f %10.3f %22s %16.2f\n',sprintf('ranked %d',summary(k,1)),summary(k,3),summary(k,4), ...
        sprintf('%.3f to %.3f',summary(k,5),summary(k,6)),summary(k,7));
end
if forecastingperiod>0
    fprintf('Rt at the end of the forecast (time %g):\n',summary(1,8));
    for k=1:size(summary,1)
        fprintf('  %-8s %10.3f %10.3f %22s %16.2f\n',sprintf('ranked %d',summary(k,1)),summary(k,9),summary(k,10), ...
            sprintf('%.3f to %.3f',summary(k,11),summary(k,12)),summary(k,13));
    end
end
fprintf(['Rt is computed from the fitted model curve under the generation interval above. The interval reflects the uncertainty ', ...
    'of the parameter estimates given the model; beyond the calibration period Rt is only as good as the forecast.\n']);

% ---- figure: one panel per ranked model
figure(600)
[rows,cols]=spatialWavePanelGrid(numel(results));
for k=1:numel(results)
    out=results{k};
    subplot(rows,cols,k)
    shown=isfinite(out.median);
    if any(shown)
        plot(out.time(shown),out.UB(shown),'r--','LineWidth',1); hold on
        plot(out.time(shown),out.LB(shown),'r--','LineWidth',1);
        plot(out.time(shown),out.median(shown),'r-','LineWidth',2);
        plot(out.time(shown),out.estimate(shown),'k-');
        line([out.time(1) out.time(end)],[1 1],'Color','k','LineStyle',':');
        line(out.atCalibrationEnd(1)*[1 1],[0 max([out.UB(shown);1])],'Color','k','LineStyle','--');
    end
    xlabel('Time'); ylabel('R_t');
    title(sprintf('Ranked model %d',k));
    set(gca,'FontSize',14);
end
end
