function [Ks,Ktot,n]=getsubepidemicsizesFinite(K0,C_thr,q,typedecline1,maxPatches)
%GETSUBEPIDEMICSIZESFINITE Sub-epidemics a spatial wave produces, and their sizes.
%   [Ks,Ktot,n] = getsubepidemicsizesFinite(K0,C_thr,q,typedecline1,maxPatches)
%
%   K0, q, typedecline1   sizes K_j of the sub-epidemics (spatialWaveSizes.m)
%   C_thr                 onset threshold: sub-epidemic j+1 starts when
%                         sub-epidemic j reaches C_thr cumulative cases
%   maxPatches            largest number of sub-epidemics allowed (default Inf)
%
%   n     number of sub-epidemics the wave ends with
%   Ks    their sizes, n-by-1
%   Ktot  sum(Ks), the size of the wave
%
%   RULE (the one simulateSpatialWave.m implements). Sub-epidemic j grows
%   towards K_j, so it reaches the threshold only if K_j > C_thr. With m the
%   number of sub-epidemics that do, the wave has m+1: the last one is
%   started by its predecessor but stays below the threshold itself.
%     C_thr > K0     the first sub-epidemic never reaches it: n = 1
%     C_thr <= 1     a seed of one case already meets it: all maxPatches
%                    sub-epidemics start at once
%     q = 0          the sizes do not decline, so the wave never ends by
%                    itself: n = maxPatches
%   When the wave has no end and maxPatches is Inf, n and Ktot are Inf and
%   Ks is empty. Ks is also empty when n exceeds one million (possible only
%   without maxPatches); Ktot is still returned.
%
%   Revision SW-1 replaces three errors of the earlier version: (1) Ks held
%   100 sizes whatever n was, and Ktot was their sum, so the reported wave
%   size did not depend on the threshold; (2) n left out the last
%   sub-epidemic, the one that starts but does not reach the threshold;
%   (3) q = 0 returned n = 18 and an infinite size.
%
%   See also spatialWaveSizes, simulateSpatialWave, spatialWaveForecastExtraPatches

if nargin<5 || isempty(maxPatches)
    maxPatches=Inf;
end
if ~(isscalar(maxPatches) && maxPatches>=1 && (maxPatches==fix(maxPatches) || isinf(maxPatches)))
    error('SpatialWave:InvalidPatchCount','maxPatches must be a positive integer or Inf.');
end
if ~(isscalar(K0) && isfinite(K0) && K0>0 && isscalar(q) && isfinite(q) && q>=0 && ...
        isscalar(C_thr) && ~isnan(C_thr) && C_thr>=0)
    error('SpatialWave:InvalidGrowthParameters', ...
        'K0 must be positive, q nonnegative and C_thr nonnegative.');
end
if ~any(typedecline1==[1 2])
    error('SpatialWave:InvalidDeclineFamily', ...
        'The decline family must be 1 (exponential) or 2 (power law).');
end

if C_thr>K0
    n=1;
elseif C_thr<=1 || q==0
    if C_thr==K0 && C_thr>1
        n=1; % q = 0 and K_j = C_thr: the threshold is approached, never reached
    else
        n=maxPatches;
    end
else
    % m = number of j >= 1 with K_j > C_thr
    switch typedecline1
        case 1 % K0*exp(-q*(j-1)) > C_thr  <=>  j < 1 + log(K0/C_thr)/q
            m=ceil(log(K0/C_thr)/q);
        case 2 % K0*j^(-q) > C_thr         <=>  j < (K0/C_thr)^(1/q)
            m=ceil((K0/C_thr)^(1/q))-1;
    end
    m=max(m,0);
    n=min(m+1,maxPatches);
end

if isinf(n)
    Ks=[];
    Ktot=Inf;
    return
end
if n>1e6
    % A slow power-law decline can give an astronomically long wave (25^10
    % sub-epidemics for K0/C_thr = 25 and q = 0.1). The sizes are not listed
    % then; their sum is computed without them.
    Ks=[];
    switch typedecline1
        case 1 % geometric series
            Ktot=K0*(-expm1(-q*n))/(-expm1(-q));
        case 2 % first 1e6 terms, then the Euler-Maclaurin tail
            N0=1e6;
            head=sum((1:N0).^(-q));
            if q==1
                integral1=log(n/N0);
            else
                integral1=(n^(1-q)-N0^(1-q))/(1-q);
            end
            Ktot=K0*(head+integral1+(n^(-q)-N0^(-q))/2);
    end
    return
end
Ks=spatialWaveSizes(K0,q,n,typedecline1);
Ktot=sum(Ks);
end
