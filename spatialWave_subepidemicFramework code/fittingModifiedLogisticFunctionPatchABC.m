% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [RMSES,PS,npatches,onset_thr,typedecline1,P]=fittingModifiedLogisticFunctionPatchABC( ...
    datafilename1,data1,DT,epidemic_period,M,flagX,typedecline2,numstartpoints)

global flag1 method1 timevect ydata yfit
global I0 npatches onset_thr typedecline1

flag1=flagX;
close all

global invasions
global timeinvasions
global Cinvasions

global npatches_fixed
global onset_fixed

global dist1
global factor1

global smoothfactor1
global calibrationperiod1

global LBe UBe

% <=============================================================================================>
% <========= Set bounds for parameters associated with error structure (alpha and d) ===========>
% <=============================================================================================>
switch method1
    case {0,1,2}
        LBe=[0 0];
        UBe=[0 0];
    case {3,4}
        LBe=[10^-8 0];
        UBe=[10^5 0];
    case 5
        LBe=[10^-8 0.2];
        UBe=[10^5 10^2];
end

% <==============================================================================>
% <============ Load data and proceed to parameter estimation ===================>
% <==============================================================================>
data1=data1(epidemic_period,:);
if data1(1,2)==0
    data1=data1(2:end,:);
end
data = data1(:,2);

% <==============================================================================>
% <============ Set time vector (timevect) and initial condition (I0) ===========>
% <==============================================================================>
timevect = (data1(:,1));
I0 = data(1); % initial condition

% <==============================================================================>
% <===================== Set initial parameter guesses ==========================>
% <==============================================================================>
r=0.2;
if flag1==3 % logistic model (p=1)
    p=1;
else
    p=0.9;
end
a=1;
K=sum(data1(:,2));
q=0.3;

if flag1==5
    r=1-I0/K;            % r=1-C0/K
    a=r/log(K/I0);       % r/log(K/C0)
end

P0=[r p a K q 1 1];

% <==============================================================================>
% <================= Set range of C_thr values (onset_thrs) =====================>
% <==============================================================================>
cumcurve1 = cumsum(smooth(data1(:,2),smoothfactor1));
onset_thrs = linspace(cumcurve1(1),cumcurve1(end),length(data1(:,2)));
onset_thrs = [0 onset_thrs(1:end-1)];

% <==============================================================================>
% <===== Set range of the possible number of subepidemics (1:npatches_fixed)=====>
% <==============================================================================>
if onset_fixed==1
    npatchess = 1:1:npatches_fixed;
elseif onset_fixed==0
    npatchess = [1 npatches_fixed];
end

if (onset_fixed==1 || all(npatchess==1))
    onset_thrs = 0;
end
onset_thrs2 = onset_thrs;

RMSES = sparse(1000,4);
PS    = sparse(1000,7);
count1 = 1;

% <================================================================================================>
% <==== Evaluate AICc across models with different number of subepidemics and C_thr values ========>
% <================================================================================================>
ydata = smooth(data,smoothfactor1);

for npatches2 = npatchess
    npatches = npatches2;

    if (onset_fixed==1 || npatches==1)
        onset_thrs = 0;
    else
        onset_thrs = onset_thrs2;
    end

    % <================================================================================================>
    % <=========================== Set initial parameter guesses and bounds ===========================>
    % <================================================================================================>
    r  = P0(1); p = P0(2); a = P0(3); K = P0(4); q = P0(5);
    alpha = P0(6); d = P0(7);

    z = [r p a K q alpha d];  % seed in (r,p,a,K,q,alpha,d) order

    rlb = mean(abs(data(1:min(2,end),1)))/200;
    rub = max(abs(data(1:min(2,end),1)))*5;

    switch flag1
        case 0 % GGM
            LB=[rlb  0 1  1 0 LBe];             % r p a K q alpha d
            UB=[rub  1 1  1 0 UBe];
        case 1 % GLM
            LB=[rlb  0 1 20 0 LBe];
            UB=[rub  1 1 1e8 5 UBe];
        case 2 % GRM
            LB=[rlb  0 0 20 0 LBe];
            UB=[rub  1 10 1e8 5 UBe];
        case 3 % Logistic
            LB=[rlb  1 1 20 0 LBe];
            UB=[rub  1 1 1e8 5 UBe];
        case 4 % Richards
            LB=[rlb  1 0 20 0 LBe];
            UB=[rub  1 10 1e8 5 UBe];
        case 5 % Gompertz
            LB=[1e-4 1 0 20 0 LBe];
            UB=[1e3  1 10 1e8 5 UBe];
    end

    % --- safety: flip any LB>UB just in case --------------------------------
    flipMask = LB > UB;
    if any(flipMask)
        tmp = LB(flipMask); LB(flipMask) = UB(flipMask); UB(flipMask) = tmp;
    end
    % ------------------------------------------------------------------------

    % === CHANGED: generate & reuse smarter MultiStart seeds (LHS + user seed) per (npatches, bounds)
    z0   = min(max(z,LB),UB);                         % clamp seed into bounds      % === CHANGED
    dDim = numel(LB);                                                             % === CHANGED
    span = (UB - LB);                                                            % === CHANGED
    nStarts = max(0, floor(numstartpoints));                                     % === CHANGED

    if nStarts > 0                                                                % === CHANGED
        if exist('lhsdesign','file') == 2                                         % === CHANGED
            X = lhsdesign(nStarts, dDim, 'criterion','maximin','iterations',50);  % === CHANGED
        else                                                                       % === CHANGED
            X = rand(nStarts, dDim);                                              % === CHANGED
        end                                                                        % === CHANGED
        starts_base = LB + X .* span;                                             % === CHANGED
        starts_base = [starts_base; z0];     % include user seed                  % === CHANGED
    else                                                                          % === CHANGED
        starts_base = z0;                                                         % === CHANGED
    end                                                                           % === CHANGED
    starts_base = unique(round(starts_base,6),'rows');                            % === CHANGED
    inB = all(starts_base >= (LB - 1e-12) & starts_base <= (UB + 1e-12), 2);      % === CHANGED
    finiteReal = all(isfinite(starts_base),2) & isreal(starts_base);              % === CHANGED
    starts_base = starts_base(inB & finiteReal, :);                                % === CHANGED
    if isempty(starts_base), starts_base = z0; end                                 % === CHANGED
    sp = CustomStartPointSet(starts_base);                                        % === CHANGED
    useParallel = ~isempty(gcp('nocreate'));                                      % === CHANGED
    % ================================================================================

    for onset_thr = onset_thrs
        for j = 1:length(typedecline2)
            typedecline1 = typedecline2(j);

            % ******** MLE estimation method  *********
            ydata = smooth(data,smoothfactor1);

            options = optimoptions('fmincon', ...
                'Algorithm','sqp', ...
                'StepTolerance',1.0e-9, ...
                'MaxFunEvals',20000, ...
                'MaxIter',20000);

            % Keep your original objective handle
            f = @plotModifiedLogisticGrowthPatchMethods1;

            % === CHANGED: MultiStart configured to reuse LHS seeds and parallel if available
            problem = createOptimProblem('fmincon', ...
                        'objective',f, 'x0',z0, 'lb',LB, 'ub',UB, 'options',options);  % === CHANGED
            ms = MultiStart('Display','off', 'UseParallel',useParallel, ...
                            'StartPointsToRun','bounds-ineqs');                         % === CHANGED
            [P,fval,flagg,outpt,allmins] = run(ms,problem,sp);                          % === CHANGED
            % ==========================================================================

            % ---- Solve dynamics with best params to validate sub-epidemic count -------
            r_hat = P(1,1);
            p_hat = P(1,2);
            a_hat = P(1,3);
            K_hat = P(1,4);
            q_hat = P(1,5);
            alpha_hat = P(1,end-1);
            d_hat     = P(1,end);

            IC = zeros(npatches,1);
            IC(1,1) = I0;
            IC(2:end,1) = 1;

            invasions   = zeros(npatches,1);
            timeinvasions = zeros(npatches,1);
            Cinvasions  = zeros(npatches,1);

            invasions(1) = 1;
            timeinvasions(1) = 0;
            Cinvasions(1) = 0;

            [~,x] = ode15s(@modifiedLogisticGrowthPatch,timevect,IC,[], ...
                           r_hat,p_hat,a_hat,K_hat,npatches,onset_thr,q_hat,flag1,typedecline1);

            if sum(invasions) < npatches
                npatches = sum(invasions);
            end

            AICc = getAICc(method1,dist1,npatches,flag1(1),1,fval,length(ydata),onset_fixed);

            RMSES(count1,:) = [npatches onset_thr typedecline1 AICc];
            PS(count1,:)    = P;
            count1 = count1 + 1;
        end % typedecline1
    end % onset_thr
end % npatches2

% <=============================================================================================>
% <======================== Sort the results by AICc (lowest to highest) =======================>
% <=============================================================================================>
RMSES = RMSES(1:count1-1,:);
PS    = PS(1:count1-1,:);

RMSES = full(RMSES);
PS    = full(PS);

[RMSES,index1] = sortrows(RMSES,[4 1]);
PS = PS(index1,:);

[~, index1] = min(RMSES(:,4));

npatches     = RMSES(index1,1);
onset_thr    = RMSES(index1,2);
typedecline1 = RMSES(index1,3);
AICmin       = RMSES(index1,4);
relativelik_i = exp((AICmin - RMSES(:,4))/2);

% <=============================================================================================>
% <============================ Get the best fit resultls ======================================>
% <=============================================================================================>
P = PS(index1,:);

r_hat=P(1,1);
p_hat=P(1,2);
a_hat=P(1,3);
K_hat=P(1,4);
q_hat=P(1,5);

alpha_hat=P(1,end-1);
d_hat=P(1,end);

if method1==3
    dist1=3;  factor1=alpha_hat; % VAR=mean+alpha*mean;
elseif method1==4
    dist1=4;  factor1=alpha_hat; % VAR=mean+alpha*mean^2;
elseif method1==5
    dist1=5;  factor1=alpha_hat; % VAR=mean+alpha*mean^2;
end

IC=zeros(npatches,1);
IC(1,1)=I0;
IC(2:end,1)=1;

invasions=zeros(npatches,1);
timeinvasions=zeros(npatches,1);
Cinvasions=zeros(npatches,1);

invasions(1)=1;
timeinvasions(1)=0;
Cinvasions(1)=0;

[~,x]=ode15s(@modifiedLogisticGrowthPatch,timevect,IC,[], ...
             r_hat,p_hat,a_hat,K_hat,npatches,onset_thr,q_hat,flag1,typedecline1);

% <=============================================================================================>
% <================================= Plot best model fit =======================================>
% <=============================================================================================>
figure(100)
for j=1:npatches
    incidence1=[x(1,j); diff(x(:,j))];
    plot(timevect,incidence1); hold on
end
y=sum(x,2);
totinc=[y(1,1); diff(y(:,1))];

if onset_thr>0
    totinc(1)=totinc(1)-(npatches-1);
end
bestfit=totinc;

plot(timevect,totinc,'r'); hold on
plot(timevect,data,'bo');
xlabel('Time (days)'); ylabel('Incidence');
title('best fit'); [npatches onset_thr];

% overdispersion estimate if needed
if (method1==0 && dist1==2)
    binsize1=7;
    [ratios,~]=getMeanVarianceRatio(data,binsize1,2);
    index1=find(ratios(:,1)>0);
    factor1=mean(ratios(index1,1));
end

'ABC estimates:'
'[npatches onset_thr typedecline1]'
[npatches onset_thr typedecline1]

% <=============================================================================================>
% <===================================  Save the results  ======================================>
% <=============================================================================================>
save(strcat('./output/ABC-original-npatchesfixed-',num2str(npatches_fixed), ...
    '-onsetfixed-',num2str(onset_fixed),'-typedecline-',num2str(sum(typedecline2)), ...
    '-smoothing-',num2str(smoothfactor1),'-',datafilename1(1:end-4), ...
    '-flag1-',num2str(flag1(1)),'-method-',num2str(method1), ...
    '-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'.mat'),'-mat')

end
