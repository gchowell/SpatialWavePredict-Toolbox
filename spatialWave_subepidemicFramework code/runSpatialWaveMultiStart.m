function [P,fval,exitflag,output,allmins]=runSpatialWaveMultiStart(ms,problem,startPoints)
%RUNSPATIALWAVEMULTISTART Run MultiStart for one candidate of the search.
%   [P,fval,exitflag,output,allmins] = runSpatialWaveMultiStart(ms,problem,startPoints)
%   is run(ms,problem,startPoints) of the Global Optimization Toolbox and
%   nothing else.
%
%   It exists so that the tests of the search
%   (tests/test_candidate_count_isolation.m) can replace the optimizer with
%   a function of this name in a temporary folder. They used to replace the
%   class MultiStart itself; MATLAB kept using the toolbox class, so those
%   tests stopped inside the real MultiStart.
[P,fval,exitflag,output,allmins]=run(ms,problem,startPoints);
end
