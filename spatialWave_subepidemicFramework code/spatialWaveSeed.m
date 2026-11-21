function info = spatialWaveSeed(label)
%SPATIALWAVESEED Apply the random seed setting and record the generator state.
%   info = spatialWaveSeed(label) is called at the start of the fit driver
%   (once per area) and of the forecast script.
%
%   Setting rngSeed (options_fit_settings_sw.m; default 1):
%     an integer   the generator is reset with rng(rngSeed), so the same call
%                  on the same files gives the same numbers again
%     empty        the generator is left as it is
%
%   info.seed   the seed used, or [] when none was set
%   info.state  the generator state at this point (the output of rng). The
%               fit driver saves it with its results.
%
%   The start points of the search, the bootstrap data and the start points
%   of the bootstrap refits are drawn on the MATLAB client, so one seed fixes
%   them with and without a parallel pool.
%
%   Revision SW-2 (after subepidemicSeed.m of the SubEpiPredict toolbox).

if nargin<1 || isempty(label)
    label = 'Spatial wave';
end
settings = spatialWaveFitSettings();
info = struct('seed',[],'state',[]);
info.seed = settings.rngSeed;
if ~isempty(info.seed)
    rng(info.seed);
    fprintf('%s: random seed %d (setting rngSeed).\n',label,info.seed);
end
info.state = rng;
end
