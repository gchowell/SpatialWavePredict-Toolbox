function [LB,UB]=getSpatialWaveBounds(flag1,data,LBe,UBe)
%GETSPATIALWAVEBOUNDS Search bounds of [r p a K q alpha d], one rule for search and bootstrap.
%   [LB,UB] = getSpatialWaveBounds(flag1,data,LBe,UBe)
%
%   flag1     growth model, 0 to 5
%   data      calibration incidence (the unsmoothed series), a vector
%   LBe, UBe  bounds of [alpha d], set by the search from method1
%
%   These are the bounds the search (fittingModifiedLogisticFunctionPatchABC.m)
%   used. The selected-model refit and the bootstrap
%   (fittingModifiedLogisticFunctionPatchMultiple.m) had their own copy,
%   which differed for the Gompertz model (r from 0, not from 1e-4) and read
%   data(1:2) without checking the length. Both now call this function.
%
%   Gompertz (flag1 = 5): a is fixed at 1. The equation r*C*log(K/C) does
%   not contain a, and get_nparams.m counts two parameters (r and K).
%
%   See also fittingModifiedLogisticFunctionPatchABC, get_nparams

data=data(:);
first=data(1:min(2,end));
rlb=mean(abs(first))/200;
rub=max(abs(first))*5;

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
        LB=[1e-4 1 1 20 0 LBe];
        UB=[1e3  1 1 1e8 5 UBe];
    otherwise
        error('SpatialWave:InvalidModelFlag','flag1 must be an integer from 0 to 5.');
end

% safety: swap any reversed pair
flipMask=LB>UB;
if any(flipMask)
    tmp=LB(flipMask); LB(flipMask)=UB(flipMask); UB(flipMask)=tmp;
end
end
