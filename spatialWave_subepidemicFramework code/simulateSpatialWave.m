function [t,x,totinc,patchIncidence,activation]=simulateSpatialWave(timevect,I0,npatches,onset_fixed,onset_thr,flag1,typedecline1,r,p,a,K,q,odeOptions,extension)
%SIMULATESPATIALWAVE Event-driven simulation of the spatial wave model.
%   [t,x,totinc,patchIncidence,activation] = simulateSpatialWave(timevect,I0, ...
%       npatches,onset_fixed,onset_thr,flag1,typedecline1,r,p,a,K,q)
%   [...] = simulateSpatialWave(...,odeOptions,extension)
%
%   The one simulator of the toolbox: the objective, the search, the accepted
%   fit, the bootstrap, the plots and the forecast all call it, so they cannot
%   disagree about initial states, activation or incidence. It follows
%   simulateSubepidemic.m of the SubEpiPredict toolbox, adapted to parameters
%   shared by all sub-epidemics and to sizes that decline along the wave.
%
%   INPUT
%     timevect      observation times, strictly increasing, starting at the
%                   initial time of the calibration series
%     I0            first observation
%     npatches      largest number of sub-epidemics the wave may have here
%     onset_fixed   0: sub-epidemic j+1 starts when sub-epidemic j reaches
%                   onset_thr cumulative cases; 1: all start at timevect(1)
%     onset_thr     the threshold C_thr (not used when onset_fixed = 1)
%     flag1         growth model, 0 to 5 (modifiedLogisticGrowthPatch.m)
%     typedecline1  1 exponential, 2 power-law decline of the sizes
%     r, p, a       growth parameters, shared by all sub-epidemics
%     K, q          size of the first sub-epidemic and decline rate
%                   (spatialWaveSizes.m)
%     odeOptions    odeset structure; [] = spatialWaveOdeOptions()
%     extension     [] or a struct that lets the wave continue beyond the
%                   fitted sub-epidemics in a forecast:
%                     extraPatches  further sub-epidemics allowed
%                     unlockTime    they can start from this time on
%                   Until unlockTime the model is exactly the fitted one
%                   (at most npatches sub-epidemics), so a forecast curve
%                   equals the fitted curve over the calibration period.
%                   The extension is not applied when it has no meaning:
%                   onset_fixed = 1, npatches = 1 (a single sub-epidemic
%                   model), or onset_thr <= 1 (every sub-epidemic starts at
%                   once, because a seed of one case already meets the
%                   threshold).
%
%   OUTPUT
%     t               timevect as a column
%     x               cumulative state at every time, one column per
%                     sub-epidemic (npatches + applied extraPatches columns)
%     totinc          incidence per time interval (spatialWaveIncidence.m)
%     patchIncidence  the same for each sub-epidemic
%     activation      struct:
%       active            final activation mask
%       nActive           number of active sub-epidemics at the last time
%       nActiveAt         the same at every row of t
%       time              activation time of each sub-epidemic (NaN: never)
%       predecessorCount  cumulative count of the predecessor at activation
%       extraPatches      sub-epidemics added by the extension
%       immediateAtUnlock true when a sub-epidemic started at unlockTime
%                         because the fitted limit npatches had held it
%                         back: the limit was binding during calibration
%       nSegments         integration segments used
%
%   ACTIVATION. Sub-epidemics are activated in order, between integration
%   segments. Each segment is integrated with a fixed set of active
%   sub-epidemics and stops at the first upward crossing of the threshold by
%   the last active one (an ode15s event); the next one is then activated and
%   a new segment starts from the event state. States are continuous: nothing
%   is added or removed at an activation. If the threshold is already met at
%   the start of a segment, the activation happens there. A sub-epidemic is
%   activated only after its predecessor (the earlier code also allowed, for
%   growth models 2 to 5, a sub-epidemic to start from the dormant seed of a
%   predecessor that had not started; that can only occur when
%   I0 < onset_thr <= 1). When onset_thr > K no second sub-epidemic starts:
%   the first one cannot reach the threshold. The fitting objective has
%   always applied that rule (also for the GGM, where K is fixed at 1, so a
%   GGM wave has one sub-epidemic unless onset_thr <= 1).
%
%   The global variables invasions, timeinvasions and Cinvasions are written
%   after a successful call, for scripts that still read them. They are never
%   read here.
%
%   Revision SW-1 (2026-11-20).
%
%   See also modifiedLogisticGrowthPatch, spatialWaveSizes,
%            spatialWaveInitialState, spatialWaveIncidence, spatialWaveOdeOptions

global invasions timeinvasions Cinvasions
invasions=[]; timeinvasions=[]; Cinvasions=[];

if nargin<13 || isempty(odeOptions)
    odeOptions=spatialWaveOdeOptions();
end
if nargin<14
    extension=[];
end

% ---- checks -------------------------------------------------------------
if ~(isnumeric(timevect) && isreal(timevect) && isvector(timevect) && ...
        ~isempty(timevect) && all(isfinite(timevect(:))) && all(diff(timevect(:))>0))
    error('SpatialWave:InvalidTimeGrid', ...
        'timevect must be a finite, real, strictly increasing time vector.');
end
if ~(isnumeric(onset_thr) && isreal(onset_thr) && isscalar(onset_thr) && ...
        isfinite(onset_thr) && onset_thr>=0)
    error('SpatialWave:InvalidThreshold','onset_thr must be finite and nonnegative.');
end
if ~(isnumeric(flag1) && isreal(flag1) && isscalar(flag1) && any(flag1==0:5))
    error('SpatialWave:InvalidModelFlag','flag1 must be an integer from 0 to 5.');
end
growth={r,p,a,K,q};
for k=1:numel(growth)
    v=growth{k};
    if ~(isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v))
        error('SpatialWave:InvalidGrowthParameters', ...
            'r, p, a, K and q must be finite real scalars.');
    end
end
if ~(isstruct(odeOptions) && isscalar(odeOptions))
    error('SpatialWave:InvalidOdeOptions','odeOptions must be an odeset structure.');
end
if ~isempty(odeget(odeOptions,'Events',[]))
    error('SpatialWave:ReservedEventsOption', ...
        'The simulator sets Events itself for the activation of sub-epidemics.');
end

% ---- extension of the wave beyond the fitted sub-epidemics ---------------
[~,~]=spatialWaveInitialState(I0,npatches,onset_fixed); % validates I0, npatches, onset_fixed
extraPatches=0;
unlockTime=-Inf;
if ~isempty(extension)
    if ~(isstruct(extension) && isscalar(extension) && isfield(extension,'extraPatches'))
        error('SpatialWave:InvalidExtension','extension must be a struct with the field extraPatches.');
    end
    extraPatches=extension.extraPatches;
    if ~(isnumeric(extraPatches) && isscalar(extraPatches) && isfinite(extraPatches) && ...
            extraPatches>=0 && extraPatches==fix(extraPatches))
        error('SpatialWave:InvalidExtension','extension.extraPatches must be a nonnegative integer.');
    end
    if isfield(extension,'unlockTime') && ~isempty(extension.unlockTime)
        unlockTime=extension.unlockTime;
        if ~(isnumeric(unlockTime) && isreal(unlockTime) && isscalar(unlockTime) && ~isnan(unlockTime))
            error('SpatialWave:InvalidExtension','extension.unlockTime must be a real scalar.');
        end
    end
    if onset_fixed==1 || npatches==1 || onset_thr<=1
        extraPatches=0; % no meaning for these models (see the help text)
    end
end
% A threshold above the size of the first (largest) sub-epidemic can never
% start a second one. The fitting objective has always applied this rule;
% now every script does, because it is applied here.
wavePossible=onset_fixed==0 && onset_thr<=K;
if ~wavePossible
    extraPatches=0;
end
capFit=npatches;
capTotal=npatches+extraPatches;

[IC,active0]=spatialWaveInitialState(I0,capFit,onset_fixed);
IC=[IC;ones(extraPatches,1)];             % extra sub-epidemics: dormant seeds
active=[logical(active0);false(extraPatches,1)];
Ks=spatialWaveSizes(K,q,capTotal,typedecline1);

t=double(timevect(:));
tEnd=t(end);
x=nan(numel(t),capTotal);
x(1,:)=IC.';
nActiveAt=nan(numel(t),1);
activationTimes=nan(capTotal,1);
activationTimes(active)=t(1);
activationCounts=nan(capTotal,1);
immediateAtUnlock=false;
currentTime=t(1);
currentState=double(IC);
nSegments=0;
wasLocked=extraPatches>0 && currentTime<unlockTime;

while true
    locked=extraPatches>0 && currentTime<unlockTime;
    if locked
        cap=capFit;
    else
        cap=capTotal;
    end
    justUnlocked=wasLocked && ~locked;
    wasLocked=locked;
    nBefore=sum(active);
    activateAtBoundary(cap);
    if justUnlocked && sum(active)>nBefore
        immediateAtUnlock=true;
    end
    nActiveAt(t==currentTime)=sum(active);
    if currentTime>=tEnd
        break
    end

    nSegments=nSegments+1;
    if nSegments>capTotal+2
        error('SpatialWave:ActivationProgressFailure', ...
            'Too many integration segments: every restart must activate a sub-epidemic.');
    end

    segmentEnd=tEnd;
    if locked && unlockTime<tEnd
        segmentEnd=unlockTime;
    end
    inner=t(t>currentTime & t<segmentEnd);
    tspan=[currentTime;inner;segmentEnd];
    if numel(tspan)==2
        % With two points ode15s returns every solver step; a third point
        % makes it return the requested times only.
        tspan=[currentTime;(currentTime+segmentEnd)/2;segmentEnd];
    end

    segmentActive=active; % fixed for the whole segment
    lastActive=find(segmentActive,1,'last');
    rhs=@(tt,xx)modifiedLogisticGrowthPatch(tt,xx,r,p,a,Ks,flag1,segmentActive);
    watching=wavePossible && lastActive<cap;

    hasEvent=false;
    try
        if watching
            eventFcn=@(tt,xx)spatialWaveActivationEvent(tt,xx,onset_thr,lastActive);
            [tt,xx,te,xe]=ode15s(rhs,tspan,currentState,odeset(odeOptions,'Events',eventFcn));
            hasEvent=~isempty(te);
        else
            [tt,xx]=ode15s(rhs,tspan,currentState,odeOptions);
        end
    catch solverError
        % A solver that gives up with an error (not a warning) is reported
        % like one that stops early, so a fit can penalize the parameter
        % vector and go on. Errors of the toolbox itself are passed on.
        if strncmp(solverError.identifier,'SpatialWave:',12)
            rethrow(solverError);
        end
        error('SpatialWave:IncompleteIntegration','ode15s failed: %s',solverError.message);
    end

    if hasEvent
        % Cut at the EARLIEST recorded event and restart from its state;
        % rows the solver returned beyond it are discarded.
        [stopTime,firstEvent]=min(te(:));
        stopState=xe(firstEvent,:).';
    else
        stopTime=tt(end);
        stopState=xx(end,:).';
        if stopTime~=segmentEnd
            error('SpatialWave:IncompleteIntegration', ...
                'ode15s stopped at t = %g, before the requested time %g.',stopTime,segmentEnd);
        end
    end
    if ~(isfinite(stopTime) && stopTime>=currentTime && stopTime<=segmentEnd) || ...
            ~isreal(stopState) || any(~isfinite(stopState(:)))
        error('SpatialWave:InvalidSegment','Invalid end point of an integration segment.');
    end

    % Requested times strictly inside the segment come from the solver
    % output; a requested time equal to the stop time gets the stop state.
    rows=find(t>currentTime & t<stopTime);
    if ~isempty(rows)
        [found,location]=ismember(t(rows),tt);
        if ~all(found)
            error('SpatialWave:IncompleteIntegration', ...
                'ode15s did not return every requested time of a segment.');
        end
        x(rows,:)=xx(location,:);
        nActiveAt(rows)=sum(segmentActive);
    end
    x(t==stopTime,:)=repmat(stopState.',sum(t==stopTime),1);

    currentTime=stopTime;
    currentState=stopState;

    if hasEvent
        % The event is the crossing by lastActive. Activate its successor
        % even if roundoff left the event state a hair below the threshold.
        activatePatch(lastActive+1);
    end
end

if ~isreal(x) || any(~isfinite(x(:)))
    error('SpatialWave:InvalidTrajectory', ...
        'The simulation must return a finite, real state at every requested time.');
end
if nargout>2
    [totinc,patchIncidence]=spatialWaveIncidence(x,onset_fixed);
end
if nargout>4
    activation=struct('active',active,'nActive',sum(active),'nActiveAt',nActiveAt, ...
        'time',activationTimes,'predecessorCount',activationCounts, ...
        'extraPatches',extraPatches,'immediateAtUnlock',immediateAtUnlock, ...
        'nSegments',nSegments);
end

% Output-only summaries for scripts that still read the globals.
invasions=double(active);
legacyTimes=activationTimes; legacyTimes(isnan(legacyTimes))=0;
legacyCounts=activationCounts; legacyCounts(isnan(legacyCounts))=0;
timeinvasions=legacyTimes;
Cinvasions=legacyCounts;

    function activateAtBoundary(capNow)
        % Thresholds already met at the start of a segment (onset_thr <= 1
        % activates the whole wave at once, as before).
        if ~wavePossible
            return
        end
        while true
            last=find(active,1,'last');
            if last>=capNow || currentState(last)<onset_thr
                break
            end
            activatePatch(last+1);
        end
    end

    function activatePatch(index)
        active(index)=true;
        activationTimes(index)=currentTime;
        activationCounts(index)=currentState(index-1);
    end
end
