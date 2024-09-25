function p = bs_pval(bs_dist,beta, h0)
% Compute two sided bootstrap p values for a certain point estimate given
    % Inputs
        % B x 1 bootstrap distribution of point estimates
        % scalar point estimate
        % h0, null hypothesis for beta = h0
    % Outputs
        % single p value

    % setup
    [B,~] = size(bs_dist);

    % compute |tau_b|
    bs_se = std(bs_dist);
    tau_b = (bs_dist - (beta - h0))./bs_se;
    abstau_b = abs(tau_b);

    % compute |tau|
    tau = (beta-h0)/bs_se;
    abstau = abs(tau);

    % compute p as the proportion of 
    p = (1/B) * sum(abstau_b > abstau);
end