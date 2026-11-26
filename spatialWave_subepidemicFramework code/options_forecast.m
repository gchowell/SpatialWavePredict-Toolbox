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
                  %  0 = Same as 1. It used to be 1/AICc, which is not a weight; the code now warns and uses Akaike weights.
                  %  1 = Weighted ensemble based on relative likelihood (Akaike weights).
                  %  2 = Weighted ensemble based on the Weighted Interval Score during calibration (WISC).
                  %  3 = Weighted ensemble based on the Weighted Interval Score of earlier forecasts (WISF);
                  %      in a single run the driver supplies WISC, so it equals 2.
                  % Codes 2 and 3 need getperformance = 1 and finite, positive scores.
