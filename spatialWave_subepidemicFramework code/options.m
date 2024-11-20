% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [cumulative1, outbreakx, caddate1, cadregion, caddisease, datatype, DT, datevecfirst1, datevecend1, numstartpoints, topmodelsx, B, flag1, typedecline2] = options

% Parameter values for the Spatial Wave Sub-Epidemic Framework.

% <============================================================================>
% <======================== Declare Global Variables ==========================>
% <============================================================================>
% Define global variables used for parameter estimation, model fitting, and data adjustments.

global method1           % Parameter estimation method:
                         % 0 = LSQ (Least Squares)
                         % 1 = MLE Poisson
                         % 2 = Pearson chi-squared
                         % 3, 4, 5 = MLE (Negative Binomial)
                         % 6 = SAD (Sum of Absolute Deviations)
global npatches_fixed    % Maximum number of sub-epidemics considered in the wave model
global onset_fixed       % Indicator for fixed onset timing of subepidemics (1=fixed, 0=not fixed)
global dist1             % Error structure assumption for the data
global smoothfactor1     % Smoothing factor for the time series data
global calibrationperiod1 % Calibration period for parameter estimation

% <============================================================================>
% <======================== Dataset Properties ================================>
% <============================================================================>
% Characteristics of the time series data used for inference:
% - File format: *.txt
% - The file may contain one or more incidence curves (one per column), 
%   where each column corresponds to new cases over time for a specific area/group.
% - Examples:
%   - "cumulative-daily-coronavirus-deaths-USA-05-11-2020.txt" (cumulative data)
%   - "daily-coronavirus-deaths-USA-05-11-2020.txt" (incidence data)

cumulative1 = 1;          % Boolean: Indicates if the data file contains cumulative counts (1) or not (0).
outbreakx = 52;           % Identifier for the spatial area or group of interest.
caddate1 = '05-11-2020';  % Date stamp of the data file (format: mm-dd-yyyy).
cadregion = 'USA';        % Geographic region (e.g., "Georgia", "USA", "World", etc.).
caddisease = 'coronavirus'; % Name of the disease associated with the data.
datatype = 'cases';       % Type of data (e.g., "cases", "deaths", "hospitalizations").

DT = 1;                   % Temporal resolution in days:
                          % 1 = daily, 7 = weekly, 365 = yearly.

% Assign temporal resolution description
if DT == 1
    cadtemporal = 'daily';
elseif DT == 7
    cadtemporal = 'weekly';
elseif DT == 365
    cadtemporal = 'yearly';
end

datevecfirst1 = [2020, 02, 27]; % Date of the first data point in the series [yyyy, mm, dd].
datevecend1 = [2021, 05, 31];   % Date of the most recent data file [yyyy, mm, dd].

% <============================================================================>
% <======================== Data Adjustments ==================================>
% <============================================================================>
% Smoothing and calibration settings for the time series data.

smoothfactor1 = 7;        % Span of the moving average for smoothing:
                          % smoothfactor1=1 applies no smoothing.

calibrationperiod1 = 90;  % Number of most recent data points used for calibration:
                          % If this value exceeds the series length, the full data length is used.

% <============================================================================>
% <=============== Parameter Estimation and Bootstrapping =====================>
% <============================================================================>
% Configuration for parameter estimation methods and error structure assumptions.

method1 = 0;              % Parameter estimation method:
                          % 0 = LSQ (Least Squares)
                          % 1 = MLE Poisson
                          % 3 = MLE Negative Binomial (VAR = mean + alpha * mean)
                          % 4 = MLE Negative Binomial (VAR = mean + alpha * mean^2)
                          % 5 = MLE Negative Binomial (VAR = mean + alpha * mean^d)
                          % 6 = SAD (Laplace distribution).

dist1 = 0;                % Error structure assumption:
                          % 0 = Normal distribution (method1 = 0)
                          % 1 = Poisson error structure (method1 = 0 or 1)
                          % 2 = Negative Binomial (VAR = factor1 * mean)
                          % 3 = MLE Negative Binomial (VAR = mean + alpha * mean)
                          % 4 = MLE Negative Binomial (VAR = mean + alpha * mean^2)
                          % 5 = MLE Negative Binomial (VAR = mean + alpha * mean^d)
                          % 6 = SAD (Laplace distribution).

numstartpoints = 10;      % Number of initial guesses for parameter optimization (MultiStart).
B = 300;                  % Number of bootstrap realizations for uncertainty characterization.

% <============================================================================>
% <========== Spatial Wave Sub-Epidemic Growth Model Configuration ============>
% <============================================================================>
% Configuration for sub-epidemic wave models.

npatches_fixed = 3;       % Maximum number of sub-epidemics considered in the model fit.

topmodelsx = 4;           % Number of best-fitting models (based on AICc) used to generate ensemble models:
                          % If npatches_fixed=1, topmodelsx is set to 1.

if npatches_fixed == 1
    topmodelsx = 1;
end

flag1 = 1;                % Growth model type:
                          % 0 = GGM (Generalized Growth Model)
                          % 1 = GLM (Logistic Model)
                          % 2 = GRM (Generalized Richards Model)
                          % 3 = LM (Linear Model)
                          % 4 = Richards Model.

onset_fixed = 0;          % Boolean: Fix onset timing of sub-epidemics:
                          % 1 = Fixed onset at time 0, 0 = Flexible onset.

typedecline2 = [2];       % Type of functional decline for sub-epidemic sizes:
                          % 1 = Exponential decline.
                          % 2 = Power-law decline.