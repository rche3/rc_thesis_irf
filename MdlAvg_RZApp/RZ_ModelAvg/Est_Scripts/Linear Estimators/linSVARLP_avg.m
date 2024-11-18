function [liny, confidencey, mult_lpu] = linSVARLP_avg(hor, p, method, clevel, nstraps, lambda, bootstrap,... % general inputs
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
    mult_lpu = zeros(3,hor);
    bs_irf_dist = zeros(nstraps, dsize, hor);
    
    if bootstrap == 1
        % get bootstrap distribution from external function
        bs_irf_dist = linSVARLP_avg_dwbse(y,beta,p,nstraps,hor, lambda, ...
            rind, sind, ...
            rpos, transformation, clevel, opt);
        
        for ri = 1:r
            for hi = 1:hor
                bs_vec = squeeze(bs_irf_dist(:,ri,hi)); % vector of bootstrapped irfs for specific reponse/horizon
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
    else
        % pass
    end
    
    % compute multipliers
    if size(liny,1) > 1

        mult_lpu(3,:) = cumsum(liny(2,:))./cumsum(liny(1,:));
        % note bs_irf_dist is nstraps x dsize x hor;
        
        mult_bs_dist = cumsum(squeeze(bs_irf_dist(:,2,:)),2)./cumsum(squeeze(bs_irf_dist(:,1,:)),2);
        
        switch method
            case "normal"
                mult_se = var(mult_bs_dist,1,1).^(1/2);
                mult_lpu(1,:) = mult_lpu(3,:) - clevel * mult_se;
                mult_lpu(2,:) = mult_lpu(3,:) + clevel * mult_se;
            case "percentile"
                alpha = 1-normcdf(clevel);
                mult_lpu(1,:) = quantile(mult_bs_dist, 1-alpha);
                mult_lpu(2,:) = quantile(mult_bs_dist, alpha);
            otherwise 
                error('Invalid method specified. Use "normal" or "percentile"')
        end
    else
    end


end

