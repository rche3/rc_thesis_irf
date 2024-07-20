function [liny,confidencey, lambda_opt]=linlp_penalised(data,x,hor,rpos,transformation, clevel, opt, nlag) 

% Code for penalised LP - code from Li et al. (2024+)

% preparations and settings for penalised LP
lambdaRange = [0.001:0.005:0.021, 0.05:0.1:1.05, 2:1:19, 20:20:100, 200:200:2000]; % cross validation grid, scaled up by T
irfLimitOrder = 2; % shrink towards polynomial of that order
CV_folds = 5; % number of folds used for cross validation

H_min = 0; % min horizon in IRF
H_max = hor-1; % max horizon in IRF
r = irfLimitOrder + 1; % order of finite difference operator in the penalty term

nT = size(data,1);
lambdaRange = lambdaRange * nT;
lambdaRange = [1e-4, lambdaRange, 1e10];

% settings to allow integration of LPW2024+ code
recursiveShock = 1; % per "observed shock" of LPW
normalizeV = 1; % we are just normalising with respect to the shock itself (let this be ordered first column in the entire matrix of observed series)
nlags = nlag;

% prep, init varibles
[dr,dsize]=size(data);
IRF = nan(dsize, hor);

for j=1:dsize
    yy = data(:, j);
    xx = x(:, rpos); % shock variable
    w = x(:, rpos+1:end); % controls

    % leave-out-out cross validation
    rss_cv = locproj_cv(yy, xx, w, H_min, H_max, r, lambdaRange, CV_folds);
    [~,lambda_opt_loc] = min(rss_cv);
    lambda_opt = lambdaRange(lambda_opt_loc); % optimally tuned lambda

    % re-estimate IRf via penalized LP using full sample
    IRF_resp = locproj(yy, xx, w, H_min, H_max, r, lambda_opt);
    % Y = [xx yy w];

    % IRF_normalize = IRF_LP(yy, xx, w, recursiveShock, normalizeV, nlags, 0); 
    IRF_normalize = 1; % the penalised LP is usually normalised at the end, but since we are normalising with respect to the (i.e.  regressing the shock on itself essentially)
    
    IRF(j,:) = (IRF_resp / IRF_normalize)';
    
    % loop to get standard errors - NOTE, THIS NEEDS TO BE FIXED
    for i=1:hor
        [r,nnn]=size(x);
    
        if transformation==1
            yy=data(i:end,j);
        elseif transformation==2
            yy=(data(i+1:end,j)-data(1:end-i,j))./data(1:end-i,2);
        end

        results=nwest(yy, x(1:end-i+1,:),i);
        reglin(:,i)=results.beta;
        if opt==0
        se(:,i)=results.se';
        else
       [EstCov, hacse, coeff]=hac_alt(x(1:end-i+1,:), yy, 'intercept', false, 'smallT', false, 'display', 'off'); 
        se(:,i)=hacse';
        end
    end

    liny(j,:)=IRF(j, :); % 
    sey(j,:)=se(rpos,:); % needs edit
    confidencey(1,:,j)=liny(j,:)-(sey(j,:)*clevel);
    confidencey(2,:,j)=liny(j,:)+(sey(j,:)*clevel);
end

end

