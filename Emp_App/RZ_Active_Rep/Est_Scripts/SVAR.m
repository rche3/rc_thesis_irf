function [liny,confidencey] = SVAR(data,x,hor,rpos,transformation, clevel, opt, nlags, bias_corrected)
    % SVAR Summary of this function goes here
    %   Detailed explanation goes here

    % VAR preparations
    res_autocorr_nlags = 2;
    normalizeV = 3; % we are normalising with one unit of shock
    recursiveShock = 3;
    Y = [data, x(:, 2)]; % construct the three-variable Y vector with gov, gdp, and shock
     
    % begin loop
    [dr,dsize]=size(data);
    for j=1:dsize
        responseV = j;
        % estimate VAR
        if bias_corrected == 0
            [Bc,By,Sigma,Sxx,Res,Beta] = VAR(Y,nlags); % no bias correction
        else
            [Bc,By,Sigma] = VAR_CorrectBias(Y,nlags); % with bias correction
        end

        G = chol(Sigma, 'lower'); % Warning: correspond to matrix C in our paper
        ShockVector = G(:,recursiveShock);
    
        % estimate IRF
        IRF = IRF_SVAR(By,ShockVector,hor-1); % IRF to one unit of shock
        IRF = IRF(responseV,:) / IRF(normalizeV,1); % normalize by response of normalization variable

        % get SEs - PLACEHOLDER, NEED FIX
        se = abs(0.01 * IRF(1));

        % store variables
        liny(j,:)=IRF; % only take the beta results
        sey(j,:)=se;
        confidencey(1,:,j)=liny(j,:)-(sey(j,:)*clevel);
        confidencey(2,:,j)=liny(j,:)+(sey(j,:)*clevel);
    end
end

