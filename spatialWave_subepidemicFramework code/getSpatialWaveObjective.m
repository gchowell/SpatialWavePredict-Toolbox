function [objfunction,yfit] = getSpatialWaveObjective(ydata,yfit,method1,alpha,d)
%GETSPATIALWAVEOBJECTIVE Shared legacy scoring of one incidence curve.
% This is the original objective arithmetic, including zero replacement,
% the negative-binomial sums (now vectorized) and omitted likelihood constants.
% It intentionally does not change the estimator or its data conventions.
% The second output is the effective curve used for scoring (zeros -> .001).

eps=0.001;

%%MLE expression
%This is the negative log likelihood, name is legacy from least squares code
%Note that a term that is not a function of the params has been excluded so to get the actual
%negative log-likliehood value you would add: sum(log(factorial(sum(casedata,2))))

if sum(yfit)==0
    objfunction=10^10;%inf;
else
    %    z
    yfit(yfit==0)=eps; %set zeros to eps to allow calculation below.  Shouldn't affect solution, just keep algorithm going.
    
    
    switch method1
        case 0
            
            %Least squares
            objfunction=sum((ydata-yfit).^2);
            
        case 1
            %MLE Poisson (negative log-likelihood)
            objfunction=-sum(ydata.*log(yfit)-yfit);
            
        case 2
            %Pearson chi squred
            objfunction=sum(((ydata-yfit).^2)./yfit);
            
        case {3,4,5}
            % MLE negative binomial (negative log-likelihood without its
            % data-only constant):
            %   3  sigma^2 = mean + alpha*mean
            %   4  sigma^2 = mean + alpha*mean^2
            %   5  sigma^2 = mean + alpha*mean^d
            % Revision SW-1: the same sums as before, evaluated with gammaln
            % and no loop over the counts (spatialWaveNegativeBinomialNLL.m,
            % from the SubEpiPredict toolbox). The loop made one pass per
            % case per observation, which was very slow for large counts.
            objfunction=spatialWaveNegativeBinomialNLL(ydata,yfit,alpha,d,method1);

        case 6 % Sum of Absolute Deviations (SAD)

            objfunction=sum(abs(ydata-yfit));

    end

    %     if ~isreal(objfunction)
    %         dbstop
    %     end
end


