function [type_GId1, mean_GI1, var_GI1, minGIcoverage] = options_Rt_sw
% OPTIONS_RT_SW  Generation interval for the effective reproduction number (Rt)
%
%   [type_GId1, mean_GI1, var_GI1, minGIcoverage] = options_Rt_sw;
%
%   Read by plotReproductionNumber_SW.m. The file name ends in _sw so that it
%   does not clash with options_Rt.m of the SubEpiPredict toolbox when both
%   are on the MATLAB path.
%
%   type_GId1      generation-interval distribution:
%                  1 = gamma, 2 = exponential, 3 = delta (fixed interval)
%   mean_GI1       mean of the generation interval, in the time unit of the
%                  data (days for daily data, WEEKS for weekly data)
%   var_GI1        its variance, in that unit squared (gamma only)
%   minGIcoverage  Rt is blank (NaN) until the curve covers this share of the
%                  generation-interval distribution (0 to 1). Before that the
%                  renewal equation misses infectors from before the first
%                  data point and overstates Rt (see spatialWaveRt.m).
%                  0 = no blanking.
%
%   SET THESE FOR YOUR PATHOGEN. The values below are the defaults of
%   options_Rt.m in the SubEpiPredict toolbox, not a recommendation for any
%   particular disease. Rt depends directly on them: a longer generation
%   interval moves every Rt value away from 1.
%
%   Revision SW-7.

type_GId1 = 1;        % 1 = gamma, 2 = exponential, 3 = delta
mean_GI1 = 11.7;      % mean of the generation interval
var_GI1 = 3.7;        % variance of the generation interval
minGIcoverage = 0.99; % share of the generation interval the curve must cover
end
