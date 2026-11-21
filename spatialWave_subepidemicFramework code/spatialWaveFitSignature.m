function [same,differences] = spatialWaveFitSignature(saved)
%SPATIALWAVEFITSIGNATURE Was a saved fit made the way a fit would be made now?
%   [same,differences] = spatialWaveFitSignature(saved)
%
%   saved   struct with the fields fitRevision and fitSettings read from the
%           file of a ranked model ([] or missing fields: a fit made before
%           these were saved).
%   same    true when the saved fit was made by the current revision of the
%           fitting layer (spatialWaveFitRevision.m) with the fit settings in
%           force now (spatialWaveFitSettings.m).
%   differences   cell array of short texts, one per difference.
%
%   Settings that do not change what a fit computes are not compared:
%   useParallel (where the same computations run), baselineWindow (a setting
%   of the forecast scores) and the settings of the identifiability report
%   (PIIscale, notIdentifiedRatio, weaklyIdentifiedPII, tradeoffCorrelation)
%   and of the intervals for model error (modelError...). A setting added
%   after revision SW-2 (weekdayAdjustment) is taken as switched off in a fit
%   saved before it existed, so such a fit stays usable while the setting is
%   off and is out of date once it is switched on.
%
%   Revision SW-3. Used by spatialWaveRollingForecast.m to decide whether the
%   saved fit of a date can be used again.

notCompared = {'useParallel','baselineWindow','PIIscale','notIdentifiedRatio','weaklyIdentifiedPII','tradeoffCorrelation', ...
    'modelErrorIntervals','modelErrorWindow','modelErrorMinErrors','modelErrorMaxScale'};
% Settings added after revision SW-2 that change the fit, with the value that
% reproduces a fit made before they existed.
addedLater = struct('weekdayAdjustment',false);
differences = {};
if ~(isstruct(saved) && isfield(saved,'fitRevision') && isfield(saved,'fitSettings') && isstruct(saved.fitSettings))
    same = false;
    differences = {'the saved fit does not record its revision and fit settings (it was made before revision SW-2)'};
    return
end
if ~isequal(saved.fitRevision,spatialWaveFitRevision())
    differences{end+1} = sprintf('fit revision %s, now %s',localText(saved.fitRevision),spatialWaveFitRevision());
end
current = spatialWaveFitSettings();
names = fieldnames(current);
for k = 1:numel(names)
    name = names{k};
    if any(strcmp(name,notCompared))
        continue
    end
    if isfield(saved.fitSettings,name)
        savedValue = saved.fitSettings.(name);
        recorded = true;
    elseif isfield(addedLater,name)
        % the saved fit was made before this setting existed, that is, with
        % the value that leaves the fit as it was
        savedValue = addedLater.(name);
        recorded = true;
    else
        recorded = false;
    end
    if recorded
        if ~isequal(savedValue,current.(name))
            differences{end+1} = sprintf('setting %s was %s, now %s',name, ...
                localText(savedValue),localText(current.(name))); %#ok<AGROW>
        end
    end
end
same = isempty(differences);
end

% -----------------------------------------------------------------------------
function t = localText(v)
if ischar(v)
    t = v;
elseif islogical(v) && isscalar(v)
    if v, t = 'true'; else, t = 'false'; end
elseif isnumeric(v) || islogical(v)
    t = mat2str(double(v));
elseif isstruct(v)
    t = sprintf('a struct with %d field(s)',numel(fieldnames(v)));
else
    t = class(v);
end
end
