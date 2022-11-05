function numparams=get_nparams(method1,npatches,flag1,fixI0)

%flag 0=GGM 1=GLM 2=GRM 3=Logistic 4=RIC 5=Gompertz

switch flag1 % model indicator
    
    case 0
        numparams=2;
        
    case 1
        numparams=3;
        
    case 2
        numparams=4;
        
    case 3
        numparams=2;
        
    case 4
        numparams=3;
        
    case 5
        numparams=2;
        
end

if npatches>1
    
    numparams=numparams+2; %q, Cthr
    
end

if fixI0==0 %fix initial datum or estimated
    
    numparams=numparams+1;
    
end


if method1==0
    
    numparams=numparams+1;
    
elseif method1==3 | method1==4 %Neg. Binomial requires one more parameter (alpha)
    
    numparams=numparams+1;
    
elseif method1==5
    
    numparams=numparams+2;  %Neg. Binomial requires 2 more parameters (alpha,d)
    
end
