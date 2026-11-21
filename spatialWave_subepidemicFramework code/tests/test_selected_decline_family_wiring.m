function tests = test_selected_decline_family_wiring
% Regression test for propagation of the selected decline family to the
% objective of the selected model and of its bootstrap refits.
%
% Revision SW-2. The family is no longer shared through a global variable:
% it is a field of the fit context that is handed to the objective. The
% original defect (the refit declared a global of another name, so the
% objective used whatever family the search had left behind) cannot occur
% with a context, and this test checks that the context is built from the
% selected row of the ranking and reaches the simulator.
tests = functiontests(localfunctions);
end

function testSelectedFamilyIsPassedToTheObjective(testCase)
testFile = mfilename('fullpath');
toolboxRoot = fileparts(fileparts(testFile));

refitText = fileread(fullfile(toolboxRoot, ...
    'fittingModifiedLogisticFunctionPatchMultiple.m'));
objectiveText = fileread(fullfile(toolboxRoot,'spatialWaveObjective.m'));

staleGlobal = regexp(refitText, '\<typefit1\>', 'once');
selectedAssignment = regexp(refitText, ...
    'typedecline1\s*=\s*RMSES\s*\(\s*index1\s*,\s*3\s*\)', 'once');
contextField = regexp(refitText, ...
    '''typedecline1''\s*,\s*typedecline1', 'once');
contextUse = regexp(refitText, ...
    'spatialWaveObjective\s*\(\s*P\s*,\s*context\s*\)', 'once');
bootstrapUse = regexp(refitText, ...
    'spatialWaveBootstrapRefits\s*\(\s*context\s*,', 'once');
objectiveUse = regexp(objectiveText, ...
    'simulateSpatialWave\([^;]*context\.typedecline1', 'once');

verifyEmpty(testCase, staleGlobal, ...
    'The obsolete typefit1 name must not reappear.');
verifyNotEmpty(testCase, selectedAssignment, ...
    'The selected decline family must be read from RMSES(:,3).');
verifyNotEmpty(testCase, contextField, ...
    'The fit context must carry the selected decline family.');
verifyNotEmpty(testCase, contextUse, ...
    'The ranked fit must be scored through the context.');
verifyNotEmpty(testCase, bootstrapUse, ...
    'The bootstrap refits must receive the same context.');
verifyNotEmpty(testCase, objectiveUse, ...
    'The objective must pass the family of its context to the simulator.');
verifyLessThan(testCase, selectedAssignment, contextField, ...
    'The selected decline family must be set before the context is built.');
verifyLessThan(testCase, contextField, contextUse, ...
    'The context must be built before the ranked fit is scored.');
end

function testFamilyOfTheContextChangesTheCurve(testCase)
% Behaviour, not text: the same parameters with the other family in the
% context give another objective value.
toolboxRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(toolboxRoot);
t = (0:60)';
context = struct('timevect',t,'ydata',[],'I0',5,'npatches',3,'onset_fixed',0, ...
    'onset_thr',400,'flag1',1,'typedecline1',1,'method1',0, ...
    'odeOptions',spatialWaveOdeOptions());
z = [0.3 0.9 1 1000 0.5 0 0];
[~,~,curve] = simulateSpatialWave(t,5,3,0,400,1,1,z(1),z(2),z(3),z(4),z(5));
context.ydata = curve;
[sameFamily,detail1] = spatialWaveObjective(z,context);
context.typedecline1 = 2;
[otherFamily,detail2] = spatialWaveObjective(z,context);
verifyTrue(testCase,detail1.success && detail2.success);
verifyLessThan(testCase,sameFamily,1e-12);
verifyGreaterThan(testCase,otherFamily,1);
end
