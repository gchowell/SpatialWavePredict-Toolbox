function tests = test_spatialWaveEnsembleCheck
% Runs spatialWaveEnsembleCheck.m (ensemble layer, revision SW-10) under runtests.
tests = functiontests(localfunctions);
end

function testEnsembleLayerPassesItsCheck(testCase)
testsFolder = fileparts(mfilename('fullpath'));
root = fileparts(testsFolder);
checkFile = fullfile(testsFolder,'spatialWaveEnsembleCheck.m');
% Say plainly when the file is not there, instead of "Unrecognized function".
assert(exist(checkFile,'file')==2, 'SpatialWave:EnsembleCheckMissing', ...
    'The file %s is missing. Copy tests/spatialWaveEnsembleCheck.m from the revision archive.', checkFile);
assert(exist(fullfile(root,'spatialWaveWeightedQuantile.m'),'file')==2, 'SpatialWave:EnsembleCheckMissing', ...
    'spatialWaveWeightedQuantile.m is not in %s: this is not the revised toolbox folder.', root);
% Run it from its own folder, so it is found whatever the path holds.
addpath(root);
previousFolder = cd(testsFolder);
restoreFolder = onCleanup(@() cd(previousFolder));
rehash;
report = spatialWaveEnsembleCheck();
verifyGreaterThanOrEqual(testCase,numel(report),38);
end
