function tests = test_spatialWaveSelfCheck
% Runs spatialWaveSelfCheck.m (simulation layer and fit layer) under runtests.
tests = functiontests(localfunctions);
end

function testSimulationLayerPassesItsSelfCheck(testCase)
testsFolder = fileparts(mfilename('fullpath'));
root = fileparts(testsFolder);
checkFile = fullfile(testsFolder,'spatialWaveSelfCheck.m');
% Say plainly when the file is not there, instead of "Unrecognized function".
assert(exist(checkFile,'file')==2, 'SpatialWave:SelfCheckMissing', ...
    'The file %s is missing. Copy tests/spatialWaveSelfCheck.m from the revision archive.', checkFile);
assert(exist(fullfile(root,'simulateSpatialWave.m'),'file')==2, 'SpatialWave:SelfCheckMissing', ...
    'simulateSpatialWave.m is not in %s: this is not the revised toolbox folder.', root);
% Run it from its own folder, so it is found whatever the path holds.
addpath(root);
previousFolder = cd(testsFolder);
restoreFolder = onCleanup(@() cd(previousFolder));
rehash;
report = spatialWaveSelfCheck();
verifyGreaterThanOrEqual(testCase,numel(report),27);
end
