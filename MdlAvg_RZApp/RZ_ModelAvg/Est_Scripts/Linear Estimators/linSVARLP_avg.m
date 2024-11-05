function [liny, confidencey] = linSVARLP_avg(hor, p, method, clevel, nstraps, lambda, ... % general inputs
    data,x,rpos,transformation, opt, ... % LP inputs
    y,rind, sind, normalise) % VAR inputs
% linSVARLP_avg computes the averaged estimator for SVAR and LP irfs 
    
    var_bootstrap = 0;
    lp_bootstrap = 0;
    [irf_var, ~, beta] = linSVAR(y,hor,p,rind,sind,clevel,var_bootstrap,nstraps,normalise, 'normal');
    [irf_lp, ~] = linlp(data,x,hor,rpos,transformation, clevel, opt, lp_bootstrap,p,nstraps, method);
    
    var_weight = lambda;
    liny = irf_var.* var_weight + irf_lp.* (1-var_weight);
    
    [dr, dsize] = size(data);

    confidencey = zeros(2,dsize,hor);
    r = length(rind);

    % get bootstrap distribution from external function
    irf_bs_dist = linSVARLP_avg_dwbse(y,beta,p,nstraps,hor, lambda, ...
        rind, sind, ...
        rpos, transformation, clevel, opt);
    
    for ri = 1:r
        for hi = 1:hor
            bs_vec = squeeze(irf_bs_dist(:,ri,hi)); % vector of bootstrapped irfs for specific reponse/horizon
            switch method
                case "normal"
                    bs_se = var(bs_vec)^(1/2);
                    irf_est = liny(ri,hi);
                    confidencey(1,ri,hi) = irf_est - clevel * bs_se;
                    confidencey(2,ri,hi) = irf_est + clevel * bs_se;
                case "percentile"
                    alpha = 1 - normcdf(clevel);
                    upperq = quantile(bs_vec, (1-alpha));
                    lowerq = quantile(bs_vec, alpha);
                    confidencey(1,ri,hi) = lowerq;
                    confidencey(2,ri,hi) = upperq;
                otherwise 
                    error('Invalid method specified. Use "normal" or "percentile"')
            end
        end
    end
end

