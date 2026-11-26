function checkSpatialWaveSettings(S)
%CHECKSPATIALWAVESETTINGS Check the settings of a fit before any fitting starts.
%   checkSpatialWaveSettings(S) returns silently when the settings can run
%   and otherwise stops with error SpatialWave:InvalidSettings, which lists
%   EVERY problem found, each with the name of the setting in options.m.
%   It follows subepidemicCheckSettings.m of the SubEpiPredict toolbox.
%
%   S is a struct with the fields
%     method1, dist1, flag1, typedecline2, npatches_fixed, onset_fixed,
%     smoothfactor1, calibrationperiod1, numstartpoints, topmodelsx, B, DT,
%     cumulative1, outbreakx
%
%   Call spatialWaveErrorStructure.m first: it sets dist1 for the methods
%   that leave no choice and reports a method/distribution pair that cannot
%   run.
%
%   See also Run_SW_subepidemicFramework, spatialWaveErrorStructure, options

problems={};
isInt=@(v) isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v) && v==fix(v);

if ~(isInt(S.method1) && any(S.method1==0:6))
    problems{end+1}='method1 must be 0, 1, 2, 3, 4, 5 or 6.';
end
if ~(isInt(S.dist1) && any(S.dist1==0:6))
    problems{end+1}='dist1 must be an integer from 0 to 6.';
end
if ~(isInt(S.flag1) && any(S.flag1==0:5))
    problems{end+1}=['flag1 must be ONE integer from 0 to 5 (0 GGM, 1 generalized logistic, ', ...
        '2 generalized Richards, 3 logistic, 4 Richards, 5 Gompertz); a vector is not supported.'];
end
if ~(isnumeric(S.typedecline2) && isvector(S.typedecline2) && ~isempty(S.typedecline2) && ...
        all(ismember(S.typedecline2,[1 2])) && numel(unique(S.typedecline2))==numel(S.typedecline2))
    problems{end+1}='typedecline2 must be 1 (exponential), 2 (power law) or [1 2].';
end
if ~(isInt(S.npatches_fixed) && S.npatches_fixed>=1)
    problems{end+1}='npatches_fixed must be a positive integer.';
end
if ~(isscalar(S.onset_fixed) && (S.onset_fixed==0 || S.onset_fixed==1))
    problems{end+1}='onset_fixed must be 0 or 1.';
end
if ~(isInt(S.smoothfactor1) && S.smoothfactor1>=1)
    problems{end+1}='smoothfactor1 must be a positive integer (1 = no smoothing).';
end
if ~(isInt(S.calibrationperiod1) && S.calibrationperiod1>=1)
    problems{end+1}='calibrationperiod1 must be a positive integer.';
end
if ~(isInt(S.numstartpoints) && S.numstartpoints>=0)
    problems{end+1}='numstartpoints must be a nonnegative integer.';
end
if ~(isInt(S.topmodelsx) && S.topmodelsx>=1)
    problems{end+1}='topmodelsx must be a positive integer.';
end
if ~(isInt(S.B) && S.B>=1)
    problems{end+1}='B (bootstrap realizations) must be a positive integer.';
end
if ~(isnumeric(S.DT) && isscalar(S.DT) && any(S.DT==[1 7 365]))
    problems{end+1}='DT must be 1 (daily), 7 (weekly) or 365 (yearly).';
end
if ~(isscalar(S.cumulative1) && (S.cumulative1==0 || S.cumulative1==1))
    problems{end+1}='cumulative1 must be 0 or 1.';
end
if ~(isnumeric(S.outbreakx) && isvector(S.outbreakx) && ~isempty(S.outbreakx) && ...
        all(isfinite(S.outbreakx)) && all(S.outbreakx>=1) && all(S.outbreakx==fix(S.outbreakx)))
    problems{end+1}='outbreakx must hold positive integers (columns of the data file).';
end
if isInt(S.flag1) && S.flag1==0 && isInt(S.npatches_fixed) && S.npatches_fixed>1
    % not an error: the fit runs, but the user should know what it means
    warning('SpatialWave:GGMWave', ...
        ['flag1 = 0 (GGM) has no sub-epidemic size, so sizes cannot decline along a wave: ', ...
         'candidates with a threshold above 1 have a single sub-epidemic.']);
end

if ~isempty(problems)
    error('SpatialWave:InvalidSettings','The settings cannot run:\n  - %s', ...
        strjoin(problems,sprintf('\n  - ')));
end
end
