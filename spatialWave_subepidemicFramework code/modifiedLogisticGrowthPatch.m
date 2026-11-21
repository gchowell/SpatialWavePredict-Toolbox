function dx=modifiedLogisticGrowthPatch(t,x,r,p,a,Ks,flag1,activeMask)
%MODIFIEDLOGISTICGROWTHPATCH Pure derivative of the spatial wave model for one segment.
%   dx = modifiedLogisticGrowthPatch(t,x,r,p,a,Ks,flag1,activeMask)
%
%   x           cumulative state of every sub-epidemic (column)
%   r, p, a     growth parameters shared by all sub-epidemics (scalars)
%   Ks          size of each sub-epidemic, from spatialWaveSizes.m (the decline
%               family and q are applied there, not here)
%   flag1       0 GGM, 1 generalized logistic, 2 generalized Richards,
%               3 logistic, 4 Richards, 5 Gompertz
%   activeMask  which sub-epidemics grow during this integration segment
%
%   The function reads and changes no global variable. Which sub-epidemics are
%   active is decided by simulateSpatialWave.m between segments, never here.
%   Until revision SW-1 the activation flags were global variables that this
%   function switched on while ode15s was integrating, also on trial steps the
%   solver later rejected, so an activation could depend on the step sequence.
%
%   Gompertz (flag1 = 5) is r*x*log(K_j/x), as in the SubEpiPredict toolbox.
%   It used to be r*x*exp(-a*t) with t the time since the START of the series:
%   every later sub-epidemic was damped by a clock it had no part in, and K
%   and q had no effect on the curve.
%
%   The earlier 11-argument call is rejected. Use simulateSpatialWave.m.
%
%   See also simulateSpatialWave, spatialWaveSizes

if nargin~=8
    error('SpatialWave:KernelInterface', ...
        ['modifiedLogisticGrowthPatch now takes (t,x,r,p,a,Ks,flag1,activeMask). ', ...
         'Call simulateSpatialWave for a complete trajectory.']);
end

dx=zeros(numel(x),1);
j=find(activeMask(:));
if isempty(j)
    return
end
xj=x(j);
Kj=Ks(j);
Kj=Kj(:);

switch flag1
    case 0 % GGM
        dx(j)=r*xj.^p;
    case 1 % generalized logistic
        dx(j)=r*(xj.^p).*(1-xj./Kj);
    case 2 % generalized Richards (existing expression retained)
        dx(j)=r*(xj.^p).*(1-xj./Kj).^a;
    case 3 % logistic
        dx(j)=r*xj.*(1-xj./Kj);
    case 4 % Richards
        dx(j)=r*xj.*(1-(xj./Kj).^a);
    case 5 % Gompertz
        dx(j)=r*xj.*log(Kj./xj);
    otherwise
        error('SpatialWave:InvalidModelFlag','flag1 must be an integer from 0 to 5.');
end
end
