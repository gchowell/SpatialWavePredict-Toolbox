function Ks=spatialWaveSizes(K0,q,npatches,typedecline1)
%SPATIALWAVESIZES Size of each sub-epidemic of a spatial wave.
%   Ks = spatialWaveSizes(K0,q,npatches,typedecline1) returns a column with
%       typedecline1 = 1 (exponential)   K_j = K0*exp(-q*(j-1))
%       typedecline1 = 2 (power law)     K_j = K0*j^(-q)
%   for j = 1..npatches. One function holds the two laws, so the simulator,
%   the fit and the reported wave size cannot disagree about them.
%
%   See also simulateSpatialWave, getsubepidemicsizesFinite

if ~(isnumeric(npatches) && isscalar(npatches) && isfinite(npatches) && npatches>=1 && npatches==fix(npatches))
    error('SpatialWave:InvalidPatchCount','npatches must be a positive integer.');
end
if ~(isnumeric(K0) && isreal(K0) && isscalar(K0) && isfinite(K0) && ...
        isnumeric(q) && isreal(q) && isscalar(q) && isfinite(q))
    error('SpatialWave:InvalidGrowthParameters','K0 and q must be finite real scalars.');
end
j=(1:npatches)';
switch typedecline1
    case 1
        Ks=K0*exp(-q*(j-1));
    case 2
        Ks=K0*j.^(-q);
    otherwise
        error('SpatialWave:InvalidDeclineFamily', ...
            'The decline family must be 1 (exponential) or 2 (power law).');
end
end
