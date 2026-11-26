function revision = spatialWaveFitRevision()
%SPATIALWAVEFITREVISION Name of the revision of the fitting layer.
%   It is saved with every fit (variable fitRevision) and changes only when a
%   revision changes what a fit computes. Revision SW-3 added forecast scoring
%   and left the fit as it was, so the name is still 'SW-2' and fits made with
%   SW-2 remain usable.
revision = 'SW-2';
end
