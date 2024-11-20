function [getperformance, deletetempfiles, forecastingperiod, weight_type1] = options_forecast

% <==============================================================================>
% <=========================== Forecasting Parameters ==========================>
% <==============================================================================>
% This section defines parameters related to the forecasting configuration.

getperformance = 1; % Boolean variable to enable/disable forecasting performance metrics:
                    % 1 = Calculate performance metrics (e.g., error, accuracy).
                    % 0 = Skip performance metric calculations.

deletetempfiles = 0; % Boolean variable to determine whether temporary forecast files should be deleted:
                     % 1 = Delete "Forecast..mat" files from the output folder after use.
                     % 0 = Retain "Forecast..mat" files for further analysis.

forecastingperiod = 30; % Integer variable specifying the forecast horizon:
                        % Number of time units to predict ahead.

% <==============================================================================>
% <================== Weighting Scheme for Ensemble Model ======================>
% <==============================================================================>
% This section defines the weighting scheme used to construct an ensemble model 
% from the top-ranked individual models.

weight_type1 = 1; % Integer variable specifying the weighting scheme for the ensemble:
                  % -1 = Equally weighted ensemble from the top models.
                  %  0 = Weighted ensemble based on AICc (Akaike Information Criterion corrected).
                  %  1 = Weighted ensemble based on relative likelihood (Akaike weights).
                  %  2 = Weighted ensemble based on the Weighted Interval Score during calibration (WISC).
