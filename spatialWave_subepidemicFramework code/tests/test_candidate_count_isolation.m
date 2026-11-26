function tests = test_candidate_count_isolation
% Candidate-count regression tests for the actual model-search function.
% The fit of one candidate (fitSpatialWaveModel.m), the simulator and
% graphics are replaced in a temporary working folder. These tests exercise
% search control flow, not fit accuracy. No test double is installed on the
% normal toolbox path.
%
% Only functions of the TOOLBOX and graphics commands are replaced; nothing
% of MATLAB's optimization toolboxes is.
%
% Revision SW-2. The search passes each candidate to the fit in an explicit
% context (number of sub-epidemics, threshold, decline family), so the
% stand-in for the fit records what it is given. Before, the candidate was
% read from global variables at the moment the optimizer ran.
tests = functiontests(localfunctions);
end

function setup(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.SourceRoot = root;
oldWarnings = warning;
testCase.addTeardown(@() warning(oldWarnings));
warning('off','MATLAB:dispatcher:nameConflict');
warning('off','MATLAB:dispatcher:ShadowedFunction');
savedGlobals = captureGlobals;
testCase.addTeardown(@() restoreGlobals(savedGlobals));

work = testCase.applyFixture( ...
    matlab.unittest.fixtures.WorkingFolderFixture);
testCase.TestData.WorkFolder = work.Folder;
mkdir(fullfile(work.Folder,'output'));
sourceFiles = {'fittingModifiedLogisticFunctionPatchABC.m','getAICc.m'};
for k = 1:numel(sourceFiles)
    copyfile(fullfile(root,sourceFiles{k}), ...
        fullfile(work.Folder,sourceFiles{k}));
end
installDoubles(work.Folder);
clear fittingModifiedLogisticFunctionPatchABC fitSpatialWaveModel simulateSpatialWave options_fit_settings_sw
rehash;
assertEqual(testCase,which('fittingModifiedLogisticFunctionPatchABC'), ...
    fullfile(work.Folder,'fittingModifiedLogisticFunctionPatchABC.m'));
for name = {'fitSpatialWaveModel','simulateSpatialWave','options_fit_settings_sw','close','figure','plot'}
    assertEqual(testCase,which(name{1}),fullfile(work.Folder,[name{1} '.m']));
end
end

function testEveryDeclineFamilyUsesConfiguredCount(testCase)
r = runSearch(3,0,[1 2]);
expected = expectedDimensions(3,0,2);
verifyEqual(testCase,r.log.fit(:,1),expected);
% In the controlled scenario, family 1 activates two components at Cthr=5;
% family 2 at the same threshold must still be fitted with three.
row = r.log.fit(:,2)==5 & r.log.fit(:,3)==2;
assertEqual(testCase,sum(row),1);
verifyEqual(testCase,r.log.fit(row,1),3);
end

function testEveryThresholdUsesConfiguredCount(testCase)
r = runSearch(3,0,1);
verifyEqual(testCase,r.log.fit(:,1),expectedDimensions(3,0,1));
% The preceding candidate at threshold 5 activates only two components.
row = r.log.fit(:,2)==10;
assertEqual(testCase,sum(row),1);
verifyEqual(testCase,r.log.fit(row,1),3);
end

function testReversingFamiliesPreservesCandidateResults(testCase)
forward = runSearch(3,0,[1 2]);
reverse = runSearch(3,0,[2 1]);
verifyEqual(testCase,canonicalResults(forward),canonicalResults(reverse), ...
    'AbsTol',1e-12);
verifyEqual(testCase,sortrows(waveRows(forward.log.fit)), ...
    sortrows(waveRows(reverse.log.fit)));
end

function testActivatedCountsRemainInRankingAndAICc(testCase)
r = runSearch(3,0,[1 2]);
fits = r.log.fit;
expected = zeros(size(fits,1),4);
for k = 1:size(fits,1)
    nActive = fits(k,4);
    expected(k,:) = [nActive fits(k,2:3) ...
        getAICc(0,1,nActive,1,1,r.log.fval(k),12,0)];
end
verifyEqual(testCase,sortrows(r.RMSES),sortrows(expected),'AbsTol',1e-12);
% Ensure the test covers both actual reductions, not just fully active fits.
verifyTrue(testCase,any(fits(:,1)==3 & fits(:,4)==2));
verifyTrue(testCase,any(fits(:,1)==3 & fits(:,4)==1));
end

function testSelectedModelAndSavedResultsRetainInterface(testCase)
r = runSearch(3,0,[1 2]);
% one single sub-epidemic candidate, and 12 thresholds for each of 2 families
verifySize(testCase,r.RMSES,[25 4]);
verifySize(testCase,r.PS,[25 7]);
verifyEqual(testCase,[r.npatches r.threshold r.family],r.RMSES(1,1:3));
verifyEqual(testCase,r.P,r.PS(1,:));
verifyEqual(testCase,r.globalCount,r.npatches);
verifyEqual(testCase,r.configuredCount,3);
verifyEqual(testCase,r.saved.RMSES,r.RMSES);
verifyEqual(testCase,r.saved.PS,r.PS);
verifyEqual(testCase,r.saved.npatches,r.npatches);
% The final simulation uses the number of sub-epidemics the best candidate
% was FITTED with, not the number that became active.
verifySize(testCase,r.log.best,[1 4]);
verifyEqual(testCase,r.log.best(1), ...
    spatialWavePatchCap(r.npatches,r.threshold,0,3));
end

function testFixedOnsetRetainsConfiguredCounts(testCase)
r = runSearch(3,1,[1 2]);
verifyEqual(testCase,r.log.fit(:,1),expectedDimensions(3,1,2));
verifyEqual(testCase,r.log.fit(:,4),r.log.fit(:,1));
verifyEqual(testCase,r.log.fit(:,2),zeros(5,1));
verifySize(testCase,r.RMSES,[5 4]);
end

function testSingleComponentConfigurationStaysSingle(testCase)
r = runSearch(1,0,[1 2]);
% One sub-epidemic is one model: it is fitted once, whatever the number of
% decline families (the family has no effect on a single sub-epidemic).
verifyEqual(testCase,r.log.fit(:,1),1);
verifyEqual(testCase,r.RMSES(:,1),1);
verifyEqual(testCase,r.npatches,1);
end

function testSingleSubEpidemicCandidateIsFittedOnce(testCase)
r = runSearch(3,0,[1 2]);
single = r.log.fit(:,1)==1;
verifyEqual(testCase,sum(single),1);
verifyEqual(testCase,r.log.fit(single,2),0);
verifyEqual(testCase,r.log.fit(single,3),1); % labelled with the first family
end

function testFullyActivatedCandidatesKeepConfiguredCounts(testCase)
r = runSearch(3,0,2);
verifyEqual(testCase,r.log.fit(:,1),expectedDimensions(3,0,1));
verifyEqual(testCase,r.log.fit(:,4),r.log.fit(:,1));
end

function testEachThresholdAlsoStartsFromThePreviousSolution(testCase)
r = runSearch(3,0,[1 2]);
fits = r.log.fit;
starts = r.log.nStarts;
% no start points besides the seed were requested (numstartpoints = 0)
verifyEqual(testCase,starts(fits(:,1)==1),1);
for family = 1:2
    rows = find(fits(:,1)==3 & fits(:,3)==family);
    verifyEqual(testCase,fits(rows,2)',0:5:55);          % thresholds in order
    verifyEqual(testCase,starts(rows(1)),1);             % first threshold: the seed
    verifyEqual(testCase,starts(rows(2:end)),2*ones(11,1)); % then seed + previous solution
end
end

function testCandidateLoopDoesNotAssignActivatedCountToGlobal(testCase)
text = fileread(fullfile(testCase.TestData.SourceRoot, ...
    'fittingModifiedLogisticFunctionPatchABC.m'));
first = strfind(text,'for npatches2 = npatchess');
last = strfind(text,'end % npatches2');
assertEqual(testCase,numel(first),1);
assertEqual(testCase,numel(last),1);
loop = text(first:last+numel('end % npatches2')-1);
assignments = regexp(loop,'(?m)^\s*npatches\s*=\s*([^;\r\n]+);','tokens');
assertEqual(testCase,numel(assignments),1);
verifyEqual(testCase,strtrim(assignments{1}{1}),'configuredNpatches');
thresholdLoop = strfind(loop,'for onset_thr = onset_thrs');
reset = strfind(loop,'npatches = configuredNpatches;');
context = strfind(loop,'''npatches'',configuredNpatches');
optimizer = strfind(loop,'[P,fval,fitInfo]=fitSpatialWaveModel(');
assertEqual(testCase,numel(context),1);
verifyGreaterThan(testCase,reset,thresholdLoop);
verifyLessThan(testCase,reset,optimizer);
verifyLessThan(testCase,context,optimizer);
end

function r = runSearch(configuredCount,onsetMode,families)
global flag1 method1 timevect ydata yfit
global I0 npatches onset_thr typedecline1
global invasions timeinvasions Cinvasions
global npatches_fixed onset_fixed dist1 factor1
global smoothfactor1 calibrationperiod1 LBe UBe
global SW_CANDIDATE_COUNT_TEST
method1 = 0;
dist1 = 1;
factor1 = 1;
smoothfactor1 = 1;
calibrationperiod1 = 12;
npatches_fixed = configuredCount;
onset_fixed = onsetMode;
% Stale initial state must not control the search.
npatches = 99;
onset_thr = -7;
typedecline1 = 9;
SW_CANDIDATE_COUNT_TEST = struct('phase','search','fit',zeros(0,4), ...
    'nStarts',zeros(0,1),'best',zeros(0,4),'fval',zeros(0,1));
data = [(0:11)' 5*ones(12,1)];
% Zero random starts makes the harness deterministic; the seed remains.
[RMSES,PS,nBest,tBest,fBest,P] = fittingModifiedLogisticFunctionPatchABC( ...
    'candidate_count_regression.txt',data,1,1:12,1,1,families,0);
r = struct('RMSES',RMSES,'PS',PS,'npatches',nBest, ...
    'threshold',tBest,'family',fBest,'P',P, ...
    'globalCount',npatches,'configuredCount',npatches_fixed, ...
    'log',SW_CANDIDATE_COUNT_TEST);
outputFile = fullfile('output',sprintf( ...
    ['ABC-original-npatchesfixed-%d-onsetfixed-%d-typedecline-%d' ...
    '-smoothing-1-candidate_count_regression-flag1-1-method-0' ...
    '-dist-1-calibrationperiod-12.mat'], ...
    configuredCount,onsetMode,sum(families)));
r.saved = load(outputFile,'RMSES','PS','npatches');
end

function expected = expectedDimensions(n,onsetMode,nFamilies)
% Configured number of sub-epidemics of every fit, in search order.
if onsetMode==1
    counts = 1:n;
elseif n==1
    counts = 1;
else
    counts = [1 n];
end
expected = zeros(0,1);
for count = counts
    nThresholds = 1;
    if onsetMode==0 && count>1
        nThresholds = 12;
    end
    familiesFitted = nFamilies;
    if count==1
        familiesFitted = 1; % one sub-epidemic is fitted once
    end
    expected = [expected; repmat(count,familiesFitted*nThresholds,1)]; %#ok<AGROW>
end
end

function rows = canonicalResults(r)
% Sort rankings and their associated parameter rows together. The single
% sub-epidemic candidate carries the label of whichever family is listed
% first, so that label is not compared.
table = [r.RMSES r.PS];
single = table(:,1)==1 & table(:,2)==0;
table(single,3) = 0;
rows = sortrows(table);
end

function rows = waveRows(fits)
% Fits of wave candidates (the single candidate's family label depends on
% the order of the families).
rows = fits(fits(:,1)>1,:);
end

function saved = captureGlobals
names = {'flag1','method1','timevect','ydata','yfit', ...
    'I0','npatches','onset_thr','typedecline1', ...
    'invasions','timeinvasions','Cinvasions', ...
    'npatches_fixed','onset_fixed','dist1','factor1', ...
    'smoothfactor1','calibrationperiod1','LBe','UBe', ...
    'SW_CANDIDATE_COUNT_TEST'};
existing = whos('global');
saved = struct('names',{names},'existed',ismember(names,{existing.name}), ...
    'values',{cell(size(names))});
for k = 1:numel(names)
    if saved.existed(k)
        eval(['global ' names{k}]);
        saved.values{k} = eval(names{k});
    end
end
end

function restoreGlobals(saved)
for k = 1:numel(saved.names)
    name = saved.names{k};
    if saved.existed(k)
        eval(['global ' name]);
        eval([name ' = saved.values{k};']);
    else
        eval(['clear global ' name]);
    end
end
end

function installDoubles(folder)
% These files exist only inside WorkingFolderFixture, never in the toolbox.

% Fixed fit settings, whatever options_fit_settings_sw.m the user has: the
% regular threshold grid alone and no extra start points.
writeDouble(folder,'options_fit_settings_sw',{
    'function U = options_fit_settings_sw()'
    'U = struct(''thresholdLowPoints'',0,''firstThresholdStartFactor'',1,''rngSeed'',[]);'
    'end'});

% Stands in for the fit of one candidate. It records the model it is given
% in the context and returns the seed. Which sub-epidemics become active
% depends only on the configured count, the threshold and the family.
writeDouble(folder,'fitSpatialWaveModel',{
    'function [P,fval,fitInfo] = fitSpatialWaveModel(context,z0,~,~,~,starts,~,~)'
    'global SW_CANDIDATE_COUNT_TEST'
    'n = context.npatches; threshold = context.onset_thr; family = context.typedecline1;'
    'active = n;'
    'if context.onset_fixed==0 && family==1 && threshold>0'
    '    active = min(n,2);'
    '    if threshold>=20, active = 1; end'
    'end'
    'P = z0;'
    'if n==1'
    '    fval = 100;'
    'else'
    '    fval = 100+6*n+family+0.01*threshold;'
    'end'
    'fitInfo = struct(''success'',true,''message'',''stand-in'',''exitflag'',1,''nActive'',active, ...'
    '    ''startAgreement'',1,''gradientInfNorm'',0,''nLocalRuns'',size(starts,1), ...'
    '    ''nConvergedRuns'',size(starts,1),''atLowerBound'',[],''atUpperBound'',[],''fallback'',[]);'
    'fitInfo.polish = struct(''run'',false,''improved'',false);'
    's = SW_CANDIDATE_COUNT_TEST;'
    's.fit(end+1,:) = [n threshold family active];'
    's.nStarts(end+1,1) = size(starts,1);'
    's.fval(end+1,1) = fval;'
    'SW_CANDIDATE_COUNT_TEST = s;'
    'end'});

% The search simulates once more, for the best model, after the loop.
writeDouble(folder,'simulateSpatialWave',{
    'function [t,x,totinc,patchIncidence,activation] = simulateSpatialWave(timevect,I0,n,onsetFixed,threshold,~,family,varargin)'
    'global SW_CANDIDATE_COUNT_TEST'
    't = timevect(:);'
    'if onsetFixed==0, IC = [I0 ones(1,n-1)]; else, IC = repmat(I0/n,1,n); end'
    'x = repmat(IC,numel(t),1);'
    'for k = 1:n'
    '    x(:,k) = IC(k)+(t-t(1))*(6*n+family+k);'
    'end'
    'patchIncidence = [x(1,:); diff(x,1,1)];'
    'if onsetFixed==0, patchIncidence(1,2:end) = patchIncidence(1,2:end)-1; end'
    'totinc = sum(patchIncidence,2);'
    'activation = struct(''nActive'',n,''active'',true(n,1));'
    'SW_CANDIDATE_COUNT_TEST.best(end+1,:) = [n threshold family n];'
    'end'});
for entry = {'close','figure','plot','hold','xlabel','ylabel','title'}
    name = entry{1};
    writeDouble(folder,name,{
        ['function varargout = ' name '(varargin)']
        'varargout = cell(1,nargout);'
        'end'});
end
end

function writeDouble(folder,name,lines)
file = fullfile(folder,[name '.m']);
fid = fopen(file,'w');
assert(fid~=-1,'CandidateCountTest:WriteFailed', ...
    'Could not create temporary test double: %s',file);
closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',lines{:});
end
