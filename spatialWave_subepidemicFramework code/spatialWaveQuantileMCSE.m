function [mcseLB, mcseUB, mcsePII] = spatialWaveQuantileMCSE(Phatss, B2)
%SPATIALWAVEQUANTILEMCSE Monte Carlo error of the 95% interval and of the PII.
%   [mcseLB, mcseUB, mcsePII] = spatialWaveQuantileMCSE(Phatss, B2)
%
%   Phatss   M-by-d bootstrap refits, one column per parameter.
%   B2       number of resamples (default 500).
%
%   The M refits of each parameter are resampled with replacement B2 times.
%   mcseLB and mcseUB are the standard deviations of the 2.5% and 97.5%
%   quantiles over the resamples, mcsePII that of log10(upper/lower) over the
%   resamples in which both are positive. They say how far the interval
%   limits and the PII would move if the bootstrap were repeated: a
%   difference between two PII values below about twice mcsePII is noise.
%   NaN for a parameter whose refits are all equal (a fixed parameter) and
%   when M < 2 or B2 < 2.
%
%   The resampling uses its own small generator (below) with a fixed seed.
%   It does not touch MATLAB's random number generator, and it gives the
%   same numbers in every MATLAB release, in GNU Octave and in the Python
%   reference of the tests. (subepidemicQuantileMCSE.m of the SubEpiPredict
%   toolbox uses a RandStream instead, so its values differ from these in
%   the digits that are noise.)
%
%   Revision SW-6.
%
%   See also spatialWaveIdentifiability, spatialWavePII

if nargin < 2 || isempty(B2)
    B2 = 500;
end
Phatss = double(full(Phatss));
[M, d] = size(Phatss);
mcseLB = nan(1, d); mcseUB = nan(1, d); mcsePII = nan(1, d);
if M < 2 || B2 < 2, return, end
idx = localIndices(M, B2);                                 % one resample per column
for k = 1:d
    x = Phatss(:, k);
    if all(isnan(x)) || max(x) - min(x) == 0, continue, end   % fixed / degenerate
    xr  = reshape(x(idx), M, B2);                          % M-by-B2 resampled estimates
    q   = quantile(xr, [0.025 0.975], 1);                  % 2-by-B2
    mcseLB(k) = std(q(1, :), 0, 2);
    mcseUB(k) = std(q(2, :), 0, 2);
    ok = q(1, :) > 0 & q(2, :) > 0;
    if nnz(ok) >= 2
        mcsePII(k) = std(log10(q(2, ok)) - log10(q(1, ok)), 0, 2);
    end
end
end

% -----------------------------------------------------------------------------
function idx = localIndices(M, B2)
% M-by-B2 indices between 1 and M from the Park-Miller generator
% s <- mod(16807*s, 2^31-1), started at 20261030 and filled column by
% column. Every product is below 2^53, so the arithmetic is exact in doubles.
modulus = 2147483647;
s = 20261030;
idx = zeros(M, B2);
for k = 1:M*B2
    s = mod(16807*s, modulus);
    idx(k) = min(floor((s - 1)/(modulus - 1)*M) + 1, M);
end
end
