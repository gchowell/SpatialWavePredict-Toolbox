function results=spatialWaveFitDiagnostic()
%SPATIALWAVEFITDIAGNOSTIC How the fit layer converges on one small problem.
%   results = spatialWaveFitDiagnostic() fits a two-sub-epidemic model (all
%   start at the first time point) to a series made from known parameters,
%   under several settings, and prints one line per fit. Nothing is saved.
%
%   Two series: the model curve itself (no noise), and the same curve with a
%   fixed +/-15% disturbance, which is what real data look like. A fit is
%   judged by its sum of squared errors against the SSE of the TRUE
%   parameters on the same series: an optimizer that has converged does at
%   least as well as the truth (ratio <= 1). How close the parameters are to
%   the truth is shown as well, but it also depends on how well the data
%   identify them.
%
%   Columns: SSE of the fit, SSE of the true parameters, their ratio, the
%   largest relative error of r, p, K, q, converged/started local runs,
%   start agreement, largest gradient component, whether the polish ran and
%   improved, seconds.

t=(0:40)'; truth=[0.3 0.9 1 1000 0.3 0 0];
[~,~,clean]=simulateSpatialWave(t,5,2,1,0,1,1,truth(1),truth(2),truth(3),truth(4),truth(5));
k=(1:numel(t))';
noisy=clean.*(1+0.15*sin(1.7*k)); noisy(1)=clean(1);

base=spatialWaveFitSettings();
base.useParallel=false;
cases={ ...
    'A  no noise, polish off, 1 random start, 1 pass',clean,false,1,1; ...
    'B  no noise, polish on,  1 random start, 1 pass',clean,true,1,1; ...
    'C  noise,    polish off, 1 random start, 1 pass',noisy,false,1,1; ...
    'D  noise,    polish on,  1 random start, 1 pass',noisy,true,1,1; ...
    'E  noise,    polish on,  5 random starts, default passes',noisy,true,5,base.fitMaxAttempts; ...
    'F  no noise, polish on,  5 random starts, default passes',clean,true,5,base.fitMaxAttempts};

fprintf('%-58s %11s %11s %9s %8s %7s %6s %9s %7s %6s\n','case','SSE fit','SSE truth','ratio','par.err','runs','agree','gradient','polish','sec');
results=struct('name',{},'P',{},'sse',{},'sseTruth',{},'parameterError',{},'info',{});
for c=1:size(cases,1)
    y=cases{c,2};
    S=base; S.polish=cases{c,3}; S.fitMaxAttempts=cases{c,5};
    [LB,UB,isLog,z0]=spatialWaveSearchSpace(y,1,0,S);
    context=struct('timevect',t,'ydata',y,'I0',5,'npatches',2,'onset_fixed',1,'onset_thr',0, ...
        'flag1',1,'typedecline1',1,'method1',0,'odeOptions',spatialWaveOdeOptions());
    sseTruth=spatialWaveObjective(truth,context);
    rng(1);
    tic;
    [P,fval,info]=fitSpatialWaveModel(context,z0,LB,UB,isLog,[],cases{c,4},S);
    seconds=toc;
    if info.success
        parameterError=max(abs(P([1 2 4 5])-truth([1 2 4 5]))./truth([1 2 4 5]));
    else
        parameterError=NaN;
    end
    polishText='no';
    if info.polish.run && info.polish.improved, polishText='better';
    elseif info.polish.run, polishText='ran'; end
    fprintf('%-58s %11.4g %11.4g %9.3g %8.2g %3d/%-3d %6.2f %9.2g %7s %6.1f\n',cases{c,1},fval,sseTruth, ...
        fval/max(sseTruth,realmin),parameterError,info.nConvergedRuns,info.nLocalRuns,info.startAgreement, ...
        info.gradientInfNorm,polishText,seconds);
    fprintf('     [r p K q] = %s   exit flag %g\n',mat2str(P([1 2 4 5]),5),info.exitflag);
    results(end+1)=struct('name',cases{c,1},'P',P,'sse',fval,'sseTruth',sseTruth, ...
        'parameterError',parameterError,'info',info); %#ok<AGROW>
end
fprintf('True [r p K q] = %s\n',mat2str(truth([1 2 4 5]),5));
end
