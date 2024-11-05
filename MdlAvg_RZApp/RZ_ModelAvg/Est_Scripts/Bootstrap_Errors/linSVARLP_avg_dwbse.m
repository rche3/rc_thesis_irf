function irf_bs_beta = linSVARLP_avg_dwbse(y,beta,p,B,hor,lambda, ...
    rind, sind, ...
    rpos, transformation, clevel, opt)
% LINSVARLP_AVG_DWBSE computes dependent wild boostrap errors for the linear SVAR/LP model average estimator
    % Outputs: 
        % irf_bs_beta which is a (B x ri x hi) matrix of bootstrapped irfs
    % Inputs:
        % 
    
    % Generate residuals and make them dependent wild to bootstrap with
    % (setup)
    [T,k] = size(y);
    lagy = lagmatrix(y, [1:p]);
    % (drop observations for all)
    lagy = lagy(p+1:end,:);
    y = y(p+1:end,:);
    TT = size(lagy,1);
    X = [ones(TT,1) lagy];

    % (generate fitted values and residuals)
    yhat = X*beta;
    res = y - yhat;

    % (convert to dependent wild)
    run('bootstrap_settings.m')

    % Prepare placeholders to fill in
    r = length(rind);
    irf_bs_beta = zeros(B, r, hor);

    for b=1:B
%         disp(['Bootstrap rep: ', num2str(b), ' out of ', num2str(B)])
        y_bs = zeros(T,k);
        y_bs(1:p,:) = y(1:p,:);

        eps = generate_dwb_resid(res, bw, dwb_settings);

        for j=p+1:T
            lagy_bs = [];
            for l=1:p
                lagy_bs = [lagy_bs, y_bs(j-l,:)];
            end
            x_bs = [1, lagy_bs];
            eps_temp = eps(j-l,:); % the first item in the eps vector corresopnds to the residual on the "(p+1)th" error term
            y_bs(j,:) = x_bs * beta + eps_temp;
        end

        data_temp = y_bs(:,rind);
        shock_temp = y_bs(:, sind);
        x_temp = lagmatrix(y_bs, [1:p]);
        
        % (drop p observations to allow for the x_temp with p lags)
        shock_temp = shock_temp(p+1:end,:);
        data_temp = data_temp(p+1:end,:);
        x_temp = x_temp(p+1:end,:);
        constant = ones(TT,1);
        xx_temp = [constant, shock_temp, x_temp]; 
        bootstrap_temp = 0;
        [lp_temp, ~, ~, ~] = linlp( ...
            data_temp, xx_temp, hor, rpos, transformation, clevel, opt, bootstrap_temp, p, 0, 'normal');
        normalise = 1;
        [svar_temp, ~, ~] = linSVAR(y_bs, hor, p, rind, sind, clevel, bootstrap_temp, B, normalise, 'normal');
        
        % store
        irf_bs_beta(b,:,:) = svar_temp.* lambda + lp_temp.*(1-lambda);
    end
end

