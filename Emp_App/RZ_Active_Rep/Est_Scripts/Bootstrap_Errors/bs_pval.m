function p = bs_pval(bs_dist,beta, h0)
% Compute two sided bootstrap p values for a certain point estimate given
    % Inputs
        % B x 1 bootstrap distribution of point estimates
        % scalar point estimate
        % h0, null hypothesis for beta = h0
    % Outputs
        % single p value

    % compute tau_b
    bs_se = std(bs_dist);
    tau_b = (beta - beta + h0)/bs_se;

    % compute p as the proportion of 
end