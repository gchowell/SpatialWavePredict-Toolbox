function finalRanking = finalizeSpatialWaveFits(stagedFiles,finalFiles,rankingFile)
%FINALIZESPATIALWAVEFITS Rerank complete fits and publish matched rank files.
% STAGEDFILES are selected-model MAT files saved in search-selection order.
% FINALFILES are the legacy output names for final ranks 1:N. No refit is
% performed. Bootstrap draws move with their own accepted model.
% Only selected/refitted models are ranked here. The ABC search file stays
% untouched. A manifest is published LAST; readers reject mixed run IDs.

N=numel(stagedFiles);
if N<1 || numel(finalFiles)~=N
    error('SpatialWave:InvalidFinalFitFiles','Expected equally sized, nonempty staging and destination lists.');
end
fits=cell(N,1);
rows=zeros(N,4);
for j=1:N
    s=load(stagedFiles{j},'fitResult','DT');
    if ~isfield(s,'fitResult')
        error('SpatialWave:MissingAcceptedFit','A staged file lacks its accepted-fit record.');
    end
    if ~isfield(s,'DT') || ~isscalar(s.DT) || ~isfinite(s.DT) || s.DT<=0
        error('SpatialWave:MissingFitTimeStep','A staged file lacks a valid DT.');
    end
    if j==1, DT=s.DT; elseif DT~=s.DT
        error('SpatialWave:IncomparableFinalFits','Selected fits must have the same DT.');
    end
    fits{j}=s.fitResult;
    f=fits{j};
    if ~isfinite(f.AICc) || ~isreal(f.AICc) || numel(f.modelSpec)~=3
        error('SpatialWave:InvalidAcceptedAICc','All selected accepted fits must have finite real scores.');
    end
    if j>1
        first=fits{1};
        if ~isequaln(f.data,first.data) || ~isequaln(f.target,first.target) || ...
                ~isequal([f.method f.distribution f.growthModel f.onsetFixed], ...
                         [first.method first.distribution first.growthModel first.onsetFixed])
            error('SpatialWave:IncomparableFinalFits', ...
                'Final rankings require identical calibration data, targets and scoring conventions.');
        end
    end
    % Do not deduplicate the FINAL rows: identical scores are not evidence
    % that two parameter vectors, trajectories or bootstrap samples coincide.
    rows(j,:)=[f.modelSpec f.AICc];
end
% Preserve selection order as an explicit last tie breaker.
[~,order]=sortrows([rows(:,4) rows(:,1) (1:N)'],[1 2 3]);
RMSES=rows(order,:);
fits=fits(order);
PS=zeros(N,7);
relativelik_i=exp(-0.5*(RMSES(:,4)-min(RMSES(:,4))));
evidenceRatio=1./relativelik_i;
[~,fitRunId]=fileparts(tempname);
for rank1=1:N
    fits{rank1}.rank=rank1;
    fits{rank1}.relativeLikelihood=relativelik_i(rank1);
    PS(rank1,:)=fits{rank1}.parameters;
end
finalRanking=struct('schemaVersion',1,'fitRunId',fitRunId, ...
    'rankingScope','selected-refitted-models', ...
    'RMSES',RMSES,'PS',PS,'relativelik_i',relativelik_i, ...
    'evidenceRatio',evidenceRatio,'fits',{fits}, ...
    'sourceSelectionOrder',order,'numModels',N,'DT',DT);

for rank1=1:N
    s=load(stagedFiles{order(rank1)});
    f=fits{rank1};
    % Retain the search-stage ranking/parameters separately for provenance.
    s.RMSES_search=s.RMSESx;
    s.PS_search=s.PS;
    s.searchRank=f.searchRank;
    s.searchAICc=f.searchAICc;
    s.P0_search=f.initialParameters;
    s.RMSESx=RMSES; s.RMSES=RMSES; s.PS=PS;
    s.relativelik_i=relativelik_i; s.evidenceRatio=evidenceRatio;
    s.P0=f.parameters; s.P=f.parameters; s.Ptrue=f.parameters;
    s.bestfit=f.bestfit; s.data1=f.data;
    s.AICc_best=f.AICc; s.RelLik_best=f.relativeLikelihood;
    s.fval=f.objectiveValue; s.numparams=f.numParameters;
    s.ydata=f.target; s.yfit=f.scoringCurve;
    if isfield(s,'AICmin'), s.AICmin_search=s.AICmin; end
    s.AICmin=RMSES(1,4);
    s.npatches=f.modelSpec(1); s.onset_thr=f.modelSpec(2);
    s.typedecline1=f.modelSpec(3); s.d=f.parameters(end);
    s.rank1=rank1; s.topmodelsx=N;
    s.fitResult=f; s.fitRunId=fitRunId;
    s.rankingScope=finalRanking.rankingScope;
    s.cadfilename1=finalFiles{rank1};
    % Staging paths are implementation details, not reproducible run inputs.
    drop=intersect(fieldnames(s),{'stageDir','stageFile','stagedFiles','finalFiles','finalRanking'});
    if ~isempty(drop), s=rmfield(s,drop); end
    writeMatAtomically(finalFiles{rank1},s);
end
writeMatAtomically(rankingFile,finalRanking);
end

function writeMatAtomically(filename,s)
folder=fileparts(filename);
if isempty(folder), folder=pwd; end
tmp=[tempname(folder) '.mat'];
cleanup=onCleanup(@() removeTemporary(tmp));
save(tmp,'-struct','s','-mat');
[ok,msg]=movefile(tmp,filename,'f');
if ~ok, error('SpatialWave:FinalFitSaveFailed','%s',msg); end
end

function removeTemporary(filename)
if exist(filename,'file')==2, delete(filename); end
end
