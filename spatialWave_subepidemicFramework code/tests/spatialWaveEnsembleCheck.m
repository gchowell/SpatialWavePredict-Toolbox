function report=spatialWaveEnsembleCheck()
%SPATIALWAVEENSEMBLECHECK Quick checks of the ensemble layer (revision SW-10).
%   report = spatialWaveEnsembleCheck() runs in MATLAB (and in GNU Octave)
%   with the toolbox folder on the path, stops with an error at the first
%   check that fails and otherwise prints one line per check. No files, no
%   figures and no fits are involved; it takes a second or two.
%
%   It checks the weighted quantile (spatialWaveWeightedQuantile.m), the
%   weights of the member models (spatialWaveEnsembleWeights.m) and the
%   weighted scores (computeQuantiles.m, computeWIS.m,
%   computeforecastperformance.m), which together make the ensemble of
%   getensemblesubepidemics.m: every curve of member k carries the weight
%   (weight of member k)/(number of curves of member k).
%
%   runtests('tests') runs this function through test_spatialWaveEnsembleCheck.m.

report={};
p=[0.01 0.025 0.05 0.1:0.1:0.9 0.95 0.975 0.99];

% deterministic numbers, the same in MATLAB and in Octave (no random generator)
T=9; N=40;
[tt,nn]=ndgrid(1:T,1:N);
X=100+50*sin(1.7*tt+0.37*nn.^1.3)+7*cos(0.11*tt.*nn)+0.01*nn;

% 1. equal weights give quantile()
q0=quantile(X,p,2);
check(maxdiff(spatialWaveWeightedQuantile(X,p,ones(1,N)),q0)<1e-10,'equal weights give quantile()');
check(maxdiff(spatialWaveWeightedQuantile(X,p),q0)<1e-10,'omitted weights give quantile()');
check(maxdiff(spatialWaveWeightedQuantile(X,p,[]),q0)<1e-10,'empty weights give quantile()');

% 2. the result does not depend on the scale of the weights or the order of the curves
w=1+mod(7*(1:N),11);
q1=spatialWaveWeightedQuantile(X,p,w);
check(maxdiff(spatialWaveWeightedQuantile(X,p,5*w),q1)<1e-10,'the scale of the weights does not matter');
order=[N:-2:1 1:2:N];
check(isequal(sort(order),1:N) && maxdiff(spatialWaveWeightedQuantile(X(:,order),p,w(order)),q1)<1e-10, ...
    'the order of the curves does not matter');

% 3. the weighted quantile against an independent calculation (no equal values in X)
q2=zeros(T,numel(p));
for i=1:T
    [xs,k]=sort(X(i,:)); ws=w(k)/sum(w);
    P=cumsum(ws)-ws/2;
    q2(i,:)=interp1(P,xs,min(max(p,P(1)),P(end)));
end
check(maxdiff(q1,q2)<1e-10,'weighted quantile equals the interpolated weighted distribution');

% 4. curves of weight 0 and NaN values are left out
w0=w; w0(3:5:N)=0;
check(maxdiff(spatialWaveWeightedQuantile(X,p,w0),spatialWaveWeightedQuantile(X(:,w0>0),p,w0(w0>0)))<1e-10, ...
    'curves of weight 0 are left out');
Xn=X; Xn(2,7)=NaN;
qn=spatialWaveWeightedQuantile(Xn,p,w);
keep=[1:6 8:N];
check(maxdiff(qn(2,:),spatialWaveWeightedQuantile(X(2,keep),p,w(keep)))<1e-10 && maxdiff(qn([1 3:T],:),q1([1 3:T],:))<1e-10, ...
    'NaN values are left out, row by row');

% 5. equal values (counts repeat): the order of the curves still does not matter
C=round(X/10);
qc=spatialWaveWeightedQuantile(C,p,w);
check(maxdiff(spatialWaveWeightedQuantile(C(:,order),p,w(order)),qc)<1e-10,'equal values: the order of the curves does not matter');
check(maxdiff(spatialWaveWeightedQuantile(C,p,ones(1,N)),quantile(C,p,2))<1e-10,'equal values and equal weights give quantile()');
check(all(all(qc>=min(C,[],2)-1e-12 & qc<=max(C,[],2)+1e-12)) && all(all(diff(qc,1,2)>=-1e-12)), ...
    'quantiles stay within the curves and increase with the probability');

% 6. the ensemble as a mixture: all curves of member k with weight W(k)/M(k)
A=X(:,1:25);          % member 1, 25 curves
B=X(:,26:40)+60;      % member 2, 15 curves, higher, overlapping member 1
perCurve=@(W) [repmat(W(1)/25,1,25) repmat(W(2)/15,1,15)];
mix=@(W) spatialWaveWeightedQuantile([A B],p,perCurve(W));
check(maxdiff(mix([1 0]),quantile(A,p,2))<1e-10,'weight 1 on a member gives that member');
check(maxdiff(mix([0 1]),quantile(B,p,2))<1e-10,'weight 1 on the other member gives the other member');
qm=mix([0.5 0.5]);
check(maxdiff(qm,quantile([A B],p,2))>1e-3,'members with different numbers of curves: equal weights are not the pooled curves');
% equal weights: half the probability on each member, whatever its number of curves
below=zeros(T,1);
for i=1:T
    below(i)=0.5*mean(A(i,:)<=qm(i,p==0.5))+0.5*mean(B(i,:)<=qm(i,p==0.5));
end
check(all(abs(below-0.5)<=0.5/15+1e-12),'the weighted median splits the mixture in halves');
% same number of curves and equal weights: simply all curves together
A15=X(:,1:15);
check(maxdiff(spatialWaveWeightedQuantile([A15 B],p,repmat(0.5/15,1,30)),quantile([A15 B],p,2))<1e-10, ...
    'equal weights and equal numbers of curves give the pooled curves');
% a member above the other: every quantile rises with its weight
Bhigh=X(:,26:40)+200;
mixHigh=@(W) spatialWaveWeightedQuantile([A Bhigh],p,perCurve(W));
m1=mixHigh([0.9 0.1]); m2=mixHigh([0.5 0.5]); m3=mixHigh([0.1 0.9]);
check(all(all(m1<=m2+1e-10 & m2<=m3+1e-10)) && maxdiff(m1,m3)>1,'more weight on the higher member raises every quantile');
% a member with a very small weight is kept (the random draw of round(M*w) curves dropped it)
ends=spatialWaveWeightedQuantile([A Bhigh],[0 1],perCurve([1-1e-6 1e-6]));
check(maxdiff(ends,[min(A,[],2) max(Bhigh,[],2)])<1e-10,'a member with a very small weight still enters the ensemble');

% 7. weights of the member models
RMSES=[3 100 1 -1500.2; 4 80 1 -1498.7; 2 150 2 -1490.1; 4 200 3 -1481.0]; % [npatches onset_thr typedecline AICc]
check(maxdiff(spatialWaveEnsembleWeights(-1,RMSES,1:4),[1;1;1;1]/4)<1e-12,'code -1: equal weights');
a=RMSES(:,4); rl=exp(-(a-min(a))/2);
wA=spatialWaveEnsembleWeights(1,RMSES,1:4);
check(maxdiff(wA,rl/sum(rl))<1e-12 && abs(sum(wA)-1)<1e-12,'code 1: Akaike weights');
shifted=RMSES; shifted(:,4)=shifted(:,4)+5000;
check(maxdiff(spatialWaveEnsembleWeights(1,shifted,1:4),wA)<1e-9,'code 1: a constant added to AICc changes nothing');
w31=spatialWaveEnsembleWeights(1,RMSES,[3 1]);
rl31=exp(-(a([3 1])-min(a([3 1])))/2);
check(maxdiff(w31,rl31/sum(rl31))<1e-12,'weights follow the order of the selected ranks');
[w00,warningId]=warningOf(@() spatialWaveEnsembleWeights(0,RMSES,1:4));
check(maxdiff(w00,wA)<1e-12,'code 0: Akaike weights (AICc below zero no longer stops the run)');
check(strcmp(warningId,'SpatialWave:WeightCodeZero'),'code 0 warns');
[~,warningId]=warningOf(@() spatialWaveEnsembleWeights(1,RMSES,1:4));
check(isempty(warningId),'code 1 does not warn');
score=[2;4;4;8];
check(maxdiff(spatialWaveEnsembleWeights(2,RMSES,1:4,score,[]),[4;2;2;1]/9)<1e-12,'code 2: proportional to 1/WISC');
check(maxdiff(spatialWaveEnsembleWeights(3,RMSES,1:4,[],score),[4;2;2;1]/9)<1e-12,'code 3: proportional to 1/WISF');
check(fails(@() spatialWaveEnsembleWeights(2,RMSES,1:4,[2;0;4;8],[]),'SpatialWave:EnsembleWeights'),'code 2: a score of 0 is an error (scoring off)');
check(fails(@() spatialWaveEnsembleWeights(2,RMSES,1:4,[2;NaN;4;8],[]),'SpatialWave:EnsembleWeights'),'code 2: a NaN score is an error');
check(fails(@() spatialWaveEnsembleWeights(2,RMSES,1:4,[2;4;8],[]),'SpatialWave:EnsembleWeights'),'code 2: one score per member');
check(fails(@() spatialWaveEnsembleWeights(4,RMSES,1:4),'SpatialWave:EnsembleWeights'),'an unknown code is an error');
check(fails(@() spatialWaveEnsembleWeights(1,RMSES,[1 9]),'SpatialWave:InvalidModelSelection'),'a rank outside the ranking is an error');

% 8. scores: equal weights give the scores of the unweighted call
tc=6;                              % calibration points; T-tc forecast points
curves1=abs(X);                    % model trajectories
curves2=abs(X+3*cos(2.3*tt+nn));   % predictive curves
data1=[(0:tc-1)' 100+10*(1:tc)'];
datalatest=[(0:T-1)' 100+10*(1:T)'];
h=T-tc;
[a1,a2,a3]=computeQuantiles(data1,curves2,h);
[b1,b2,b3]=computeQuantiles(data1,curves2,h,ones(1,N)/N);
check(maxdiff(a1,b1)<1e-9 && maxdiff(a2,b2)<1e-9 && isequal(a3,b3),'computeQuantiles: equal weights give the unweighted quantiles');
[WISCa,WISFa,qca,qfa]=quiet(@computeWIS,data1,datalatest,curves2,h);
[WISCb,WISFb,qcb,qfb]=quiet(@computeWIS,data1,datalatest,curves2,h,ones(1,N)/N);
check(maxdiff(WISCa,WISCb)<1e-9 && maxdiff(WISFa,WISFb)<1e-9 && maxdiff(qca,qcb)<1e-9 && maxdiff(qfa,qfb)<1e-9, ...
    'computeWIS: equal weights give the unweighted scores');
sa=cell(1,10); sb=cell(1,10);
[sa{1:10}]=quiet(@computeforecastperformance,data1,datalatest,curves1,curves2,h);
[sb{1:10}]=quiet(@computeforecastperformance,data1,datalatest,curves1,curves2,h,ones(1,N)/N,ones(1,N)/N);
same=true;
for i=1:10
    same=same && maxdiff(sa{i},sb{i})<1e-9;
end
check(same,'computeforecastperformance: equal weights give the unweighted scores');

% 9. scores: all the weight on a member gives the scores of that member alone
wm=[ones(1,25)/25 zeros(1,15)];
[WISCm,WISFm]=quiet(@computeWIS,data1,datalatest,curves2,h,wm);
[WISC1,WISF1]=quiet(@computeWIS,data1,datalatest,curves2(:,1:25),h);
check(maxdiff(WISCm,WISC1)<1e-9 && maxdiff(WISFm,WISF1)<1e-9,'computeWIS: all the weight on a member gives the WIS of that member');
sm=cell(1,10); s1=cell(1,10);
[sm{1:10}]=quiet(@computeforecastperformance,data1,datalatest,curves1,curves2,h,wm,wm);
[s1{1:10}]=quiet(@computeforecastperformance,data1,datalatest,curves1(:,1:25),curves2(:,1:25),h);
same=true;
for i=1:10
    same=same && maxdiff(sm{i},s1{i})<1e-9;
end
check(same,'computeforecastperformance: all the weight on a member gives the scores of that member');

% 10. scores: a weighted ensemble scores differently from the pooled curves
wu=perCurve([0.8 0.2]);
[WISCu,WISFu]=quiet(@computeWIS,data1,datalatest,[A B],h,wu);
[WISCp,WISFp]=quiet(@computeWIS,data1,datalatest,[A B],h);
check(maxdiff(WISCu,WISCp)>1e-6 && maxdiff(WISFu,WISFp)>1e-6,'computeWIS: the weights reach the score');

fprintf('spatialWaveEnsembleCheck: %d checks passed.\n',numel(report));

    function check(condition,name)
        if ~condition
            error('SpatialWave:EnsembleCheck','Check failed: %s.',name);
        end
        report{end+1}=name; %#ok<AGROW>
        fprintf('  ok  %s\n',name);
    end
end

% -----------------------------------------------------------------------------
function d=maxdiff(a,b)
% largest absolute difference; Inf when the sizes differ or a NaN is not matched
if ~isequal(size(a),size(b)) || ~isequal(isnan(a),isnan(b))
    d=Inf;
    return
end
difference=abs(a(:)-b(:));
difference=difference(~isnan(difference));
if isempty(difference)
    d=0;
else
    d=max(difference);
end
end

function yes=fails(call,identifier)
% true when the call stops with the given error identifier
yes=false;
try
    call();
catch problem
    yes=strcmp(problem.identifier,identifier);
end
end

function varargout=quiet(f,varargin)
% calls f without its console output (the scoring functions print as they go)
% (one output of evalc only: Octave returns no more)
varargout=cell(1,nargout);
evalc('[varargout{1:nargout}]=f(varargin{:});');
end

function [value,identifier]=warningOf(call)
% value of the call and the identifier of the warning it gave ('' when none).
% The warning is left on: Octave does not record a warning that is turned off.
lastwarn('');
value=[];
evalc('value=call();');
[~,identifier]=lastwarn;
end
