% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [cumulative1, outbreakx, caddate1, cadregion, caddisease,datatype, DT, datevecfirst1, datevecend1, numstartpoints,topmodelsx, M,flag1,typedecline2]=options

% parameter values for the Spatial wave sub-epidemic framework.

% <============================================================================>
% <=================== Declare global variables =======================================>
% <============================================================================>

global method1 %LSQ=0, MLE Poisson=1, Pearson chi-squared=2, MLE (Neg Binomial)=3,MLE (Neg Binomial)=4, MLE (Neg Binomial)=5
global npatches_fixed
global onset_fixed
global dist1
global factor1
global smoothfactor1
global calibrationperiod1

% <============================================================================>
% <================================ Datasets properties =======================>
% <============================================================================>
% Located in the input folder, the time series data file is a text file with the extension *.txt. The data file can contain one or more incidence curves (one per
% column in the file). Each column corresponds to the number of new cases over time for each epidemic corresponding to a different area/group.
% For instance, each column could correspond to different states in
% the U.S or countries in the world. In the options.m file, a specific data column in the file can be accessed using the parameter <outbreakx> (see below).

% if the time series file contains cumulative incidence count data, the name of the time series data file starts with "cumulative" with the
% following format:

% 'cumulative-<cadtemporal>-<caddisease>-<datatype>-<cadregion>-<caddate1>.txt');
%  For example: 'cumulative-daily-coronavirus-deaths-USA-05-11-2020.txt'

% Otherwise, if the time series file contains incidence data, the name of the data file follows the format:

% <cadtemporal>-<caddisease>-<datatype>-<cadregion>-<caddate1>.txt');
%  For example: 'daily-coronavirus-deaths-USA-05-11-2020.txt'

cumulative1=1; % flag to indicate if the data file contains cumulative incidence counts (cumulative1=1) or not (cumulative1=0)

outbreakx=52;  % identifier for the spatial area of interest

caddate1='05-11-2020';  % data file time stamp in format: mm-dd-yyyy

cadregion='USA'; % string indicating the geographic region of the time series contained in the file (Georgia, USA, World, Asia, Africa, etc.)

caddisease='coronavirus'; % string indicating the name of the disease related to the time series data

datatype='cases'; % string indicating the nature of the data (cases, deaths, hospitalizations, etc)

DT=1; % temporal resolution in days (1=daily data, 7=weekly data).

if DT==1
    cadtemporal='daily'; % string indicating the temporal resolution of the data
elseif DT==7
    cadtemporal='weekly';
end

datevecfirst1=[2020 02 27]; % date corresponding to the first data point in time series data in format [year_number month_number day_number]

datevecend1=[2021 05 31]; % date of the most recent data file in format [year_number month_number day_number]. This data file is accessed to assess forecast performance

% <============================================================================>
% <============================Adjustments to data =================================>
% <============================================================================>

smoothfactor1=1; % The span of the moving average smoothing of the case series (smoothfactor1=1 indicates no smoothing)

calibrationperiod1=90; % calibrates model using the most recent <calibrationperiod1> data points where <calibrationperiod> does not exceed the length of the time series data otherwise it will use the maximum length of the data

% <=============================================================================>
% <=========================== Parameter estimation and bootstrapping=====================>
% <=============================================================================>

method1=0; % Type of estimation method. See below:

% Nonlinear least squares (LSQ)=0,
% MLE Poisson=1,
% MLE (Neg Binomial)=3, with VAR=mean+alpha*mean;
% MLE (Neg Binomial)=4, with VAR=mean+alpha*mean^2;
% MLE (Neg Binomial)=5, with VAR=mean+alpha*mean^d;

dist1=0; % Define dist1 which is the type of error structure. See below:

%dist1=0; % Normal distribution to model error structure (method1=0)
%dist1=1; % Poisson error structure (method1=0 OR method1=1)
%dist1=2; % Neg. binomial error structure where var = factor1*mean where
                  % factor1 is empirically estimated from the time series
                  % data (method1=0)
%dist1=3; % MLE (Neg Binomial) with VAR=mean+alpha*mean  (method1=3)
%dist1=4; % MLE (Neg Binomial) with VAR=mean+alpha*mean^2 (method1=4)
%dist1=5; % MLE (Neg Binomial)with VAR=mean+alpha*mean^d (method1=5)

numstartpoints=10; % Number of initial guesses for optimization procedure using MultiStart

M=300; % number of bootstrap realizations to characterize parameter uncertainty

% <==============================================================================>
% <========================= Spatial wave sub-epidemic model ============================>
% <==============================================================================>

npatches_fixed=4; % maximum number of subepidemics considered in epidemic wave model fit

topmodelsx=4; % number of best fitting models (based on AICc) that will be generated to derive ensemble models

if npatches_fixed==1  % if one sub-epidemic is employed, then there is only one model
    topmodelsx=1;
end

flag1=1; % Type of growth model used to model a subepidemic

% 0 = GGM
% 1 = GLM
% 2 = GRM
% 3 = LM
% 4 = Richards

onset_fixed=0; % flag to indicate if the onset timing of subepidemics fixed at time 0 (onset_fixed=1) or not (onset_fixed=0).

typedecline2=[2]; % Type of functional declines that will be considered for the sequential sub-epidemic sizes where 1=exponential decline in subepidemic size; 2=power-law decline in subepidemic size

