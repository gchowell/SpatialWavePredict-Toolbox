function [objfunction,yfit] = getSpatialWaveObjective(ydata,yfit,method1,alpha,d)
%GETSPATIALWAVEOBJECTIVE Shared legacy scoring of one incidence curve.
% This is the original objective arithmetic, including zero replacement,
% count-wise negative-binomial sums and omitted likelihood constants.
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
            
        case 3
            % MLE Negative binomial (negative log-likelihood) where sigma^2=mean+alpha*mean;
            
            sum1=0;
            
            for i=1:length(ydata)
                for j=0:(ydata(i)-1)
                    
                    sum1=sum1+log(j+(1/alpha)*yfit(i));
                    
                end
                sum1=sum1+ydata(i)*log(alpha)-(ydata(i)+(1/alpha)*yfit(i))*log(1+alpha);
                
            end
            
            objfunction=-sum1;
            
        case 4
            % MLE Negative binomial (negative log-likelihood) where sigma^2=mean+alpha*mean^2;
            
            sum1=0;
            
            for i=1:length(ydata)
                for j=0:(ydata(i)-1)
                    
                    sum1=sum1+log(j+(1/alpha));
                    
                end
                
                sum1=sum1+ydata(i)*log(alpha*yfit(i))-(ydata(i)+(1/alpha))*log(1+alpha*yfit(i));
                
            end
            
            objfunction=-sum1;
           

        case 5
            % MLE Negative binomial (negative log-likelihood) where sigma^2=mean+alpha*mean^d;

            sum1=0;

            for i=1:length(ydata)
                for j=0:(ydata(i)-1)

                    sum1=sum1+log(j+(1/alpha)*yfit(i).^(2-d));

                end

                sum1=sum1+ydata(i)*log(alpha*(yfit(i).^(d-2)).*yfit(i))-(ydata(i)+(1/alpha)*yfit(i).^(2-d))*log(1+alpha*(yfit(i).^(d-2)).*yfit(i));

            end

            objfunction=-sum1;

        case 6 % Sum of Absolute Deviations (SAD)

            objfunction=sum(abs(ydata-yfit));

    end

    %     if ~isreal(objfunction)
    %         dbstop
    %     end
end


