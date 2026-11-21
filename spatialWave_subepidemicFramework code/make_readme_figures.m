function files = make_readme_figures(outbreakx,caddate1,forecastingperiod,weight_type1,folder)
%MAKE_README_FIGURES Save the example figures shown in the README.
%   make_readme_figures
%   files = make_readme_figures(outbreakx,caddate1,forecastingperiod,weight_type1,folder)
%
%   Runs the plotting functions of the toolbox on a fit that is already in
%   ./output and saves their figures as PNG files with the names the README
%   refers to. Nothing is fitted: run Run_SW_subepidemicFramework first, with
%   the settings of the README's quick start (or your own).
%
%   outbreakx, caddate1   area and date ('mm-dd-yyyy'); [] or omitted: options.m
%   forecastingperiod     [] or omitted: options_forecast.m
%   weight_type1          [] or omitted: options_forecast.m
%   folder                where the files go; default ../docs/images (the
%                         docs folder at the repository root, when MATLAB's
%                         current folder is the code directory)
%
%   Files written (those whose figure the settings produce):
%     relativelikelihood.png      AICc, relative likelihood and evidence ratio of the ranked models
%     rankings.png                the ranked models with their sub-epidemic profiles
%     model_fit.png               ranked model 1: parameter distributions, fit, sub-epidemics, residuals
%     performance_calibration.png calibration performance of the ranked models
%     parameters.png              parameter estimates across the ranked models
%     forecasts.png               forecasts of the ranked models
%     ensembles.png               ensemble forecasts
%     performance_forecast.png    forecast performance of the ranked models (needs getperformance = 1)
%     performance_ensemble.png    forecast performance of the ensembles (needs getperformance = 1)
%     reproductionnumber.png      effective reproduction number of the ranked models
%
%   The reproduction number uses the generation interval in options_Rt_sw.m:
%   set it for your pathogen before publishing that figure.
%
%   files   cell array of the files written.

if nargin<1, outbreakx=[]; end
if nargin<2, caddate1=[]; end
if nargin<3, forecastingperiod=[]; end
if nargin<4, weight_type1=[]; end
if nargin<5 || isempty(folder)
    folder=fullfile('..','docs','images');
end
if ~exist(folder,'dir')
    mkdir(folder);
end
files={};

close all
plotRankings_SW_subepidemicFramework(outbreakx,caddate1);
files=localSave(files,99,folder,'relativelikelihood.png',[100 100 1500 450]);
files=localSave(files,100,folder,'rankings.png',[100 100 1500 900]);

close all
plotFit_SW_subepidemicFramework(outbreakx,caddate1);
files=localSave(files,101,folder,'model_fit.png',[100 100 1700 850]);
files=localSave(files,300,folder,'performance_calibration.png',[100 100 1200 800]);
files=localSave(files,400,folder,'parameters.png',[100 100 1200 900]);

close all
plotForecast_SW_subepidemicFramework(outbreakx,caddate1,forecastingperiod,weight_type1);
files=localSave(files,10,folder,'forecasts.png',[100 100 1500 900]);
files=localSave(files,11,folder,'ensembles.png',[100 100 1500 900]);
files=localSave(files,400,folder,'performance_forecast.png',[100 100 1200 800]);
files=localSave(files,401,folder,'performance_ensemble.png',[100 100 1200 800]);

close all
plotReproductionNumber_SW(outbreakx,caddate1,forecastingperiod);
files=localSave(files,600,folder,'reproductionnumber.png',[100 100 1500 900]);

fprintf('\n%d figure file(s) written to %s:\n',numel(files),folder);
fprintf('  %s\n',files{:});
end

% -----------------------------------------------------------------------------
function files = localSave(files,number,folder,name,position)
% Save figure NUMBER if the plotting function made it.
if isempty(findobj(0,'Type','figure','Number',number))
    fprintf('Figure %d was not produced with these settings; %s is not written.\n',number,name);
    return
end
handle=figure(number);
set(handle,'Color','w','Position',position);
drawnow
target=fullfile(folder,name);
if exist('exportgraphics','file')
    exportgraphics(handle,target,'Resolution',150);
else
    set(handle,'PaperPositionMode','auto');
    print(handle,target,'-dpng','-r150');
end
files{end+1,1}=target;
end
