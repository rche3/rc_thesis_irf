function [liny,confidencey] = SVAR(data,x,hor,rpos,transformation, clevel, opt, bias_corrected)
    %BVAR Summary of this function goes here
    %   Detailed explanation goes here
    
    % preparations
    res_autocorr_nlags = 2
    
    
    % estimate VAR
    if bias_corrected == 0
        [Bc,By,Sigma,Sxx,Res,Beta] = VAR(Y,nlags); % no bias correction
    else
        [Bc,By,Sigma] = VAR_CorrectBias(Y,nlags); % with bias correction
    end

    G = chol(Sigma, 'lower'); % Warning: correspond to matrix C in our paper
    ShockVector = G(:,recursiveShock);

    outputArg1 = inputArg1;
    outputArg2 = inputArg2;
end

