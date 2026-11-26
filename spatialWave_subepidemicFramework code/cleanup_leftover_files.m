function cleanup_leftover_files(doIt)
%CLEANUP_LEFTOVER_FILES Move files the toolbox no longer uses out of the code directory.
%   cleanup_leftover_files          lists what would be moved; changes nothing
%   cleanup_leftover_files(true)    moves the files
%
%   Run it once, in the code directory, after installing revision SW-9. The
%   revision's own archive already has this layout; this function brings an
%   existing folder to it. Nothing is deleted:
%
%   - Earlier configuration files (options_covid19.m, options_tutorial.m, ...)
%     go to ./option_presets. Each defines a function called "options"
%     whatever the file is called, so next to the real options.m they are a
%     source of confusion. See option_presets/README.md.
%   - Files no code refers to (an earlier copy of the search function, an
%     editor autosave) go to ./removed_files. Delete that folder when you
%     have looked at it.
%
%   options.m, options_forecast.m, options_Rt_sw.m, options_fit_settings_sw.m
%   and its example are not touched.
if nargin<1
    doIt=false;
end
presets={'options_SARS.m','options_SARS2003Singapore.m','options_covid19.m','options_fmd.m','options_mpox.m', ...
    'options_mpx.m','options_mpx2.m','options_sw.m','options_tutorial.m','options_yearly.m','optionsx.m','options_forecast_sw.m'};
unused={'fittingModifiedLogisticFunctionPatchABC_old_09042025.m','getensemblesubepidemics.asv', ...
    'cumulative-yearly-monkeypox-cases-USA-07-20-2079.txt'};
if exist('Run_SW_subepidemicFramework.m','file')~=2 || exist(fullfile(pwd,'options.m'),'file')~=2
    error('SpatialWave:Cleanup','Run this in the code directory (the folder with Run_SW_subepidemicFramework.m and options.m).');
end
nMoved=localMove(presets,'option_presets',doIt);
nMoved=nMoved+localMove(unused,'removed_files',doIt);
% any other editor autosave
autosaves=dir('*.asv');
nMoved=nMoved+localMove({autosaves.name},'removed_files',doIt);
if nMoved==0
    fprintf('Nothing to move: the folder is already tidy.\n');
elseif doIt
    fprintf('%d file(s) moved. Run "rehash; clear functions" before the next command.\n',nMoved);
else
    fprintf('%d file(s) would be moved. Run cleanup_leftover_files(true) to move them.\n',nMoved);
end
end

% -----------------------------------------------------------------------------
function n = localMove(names,folder,doIt)
n=0;
for k=1:numel(names)
    if exist(fullfile(pwd,names{k}),'file')~=2
        continue
    end
    n=n+1;
    if doIt
        if ~exist(folder,'dir')
            mkdir(folder);
        end
        movefile(fullfile(pwd,names{k}),fullfile(pwd,folder,names{k}));
        fprintf('moved   %s  ->  %s/\n',names{k},folder);
    else
        fprintf('would move   %s  ->  %s/\n',names{k},folder);
    end
end
end
