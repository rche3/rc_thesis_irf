function [liny,confidencey, bs_beta_means, bs_beta_dist]=linlp( ...
    data,x,hor,rpos,transformation, clevel, opt, bootstrap, nlag, nstraps, method) 

[dr,dsize]=size(data);
confidencey = zeros(2,dsize,hor);
[~, nvar]=size(x);
bs_beta_means = zeros(dsize,hor);
bs_ci = zeros(2,dsize,hor);
bs_beta_dist = zeros(nstraps, dsize, hor);
    % pass

for j=1:dsize
    y_pos_control = j+1;
    [r,nnn]=size(x);
    for i=1:hor
        if transformation==1
            yy=data(i:end,j);
        elseif transformation==2
            yy=(data(i+1:end,j)-data(1:end-i,j))./data(1:end-i,2);
        end

        results=nwest_rc(yy, x(1:end-i+1,:),i);
        reglin(:,i)=results.beta;
        irf_est = reglin(rpos, i);

        % compute SEs
        if bootstrap == 1
            ctrlstart = rpos + 1; % the controls right after the shock position
            bs_beta = linlp_dwbse(results.resid, yy, x(1: end-i +1, :), results.beta, i, nlag, y_pos_control, ctrlstart, nstraps);
            irf_bs_dist = bs_beta(rpos, :); % returns 1 x B bootstrapped estimates
            bs_beta_dist(:,j,i) = irf_bs_dist';
            
            switch method
                case "normal"
                    bs_se = var(irf_bs_dist)^(1/2);
                    bs_ci(1,j,i) = irf_est - clevel * bs_se;
                    bs_ci(2,j,i) = irf_est + clevel * bs_se;
                case "percentile"
                    upperq = quantile(irf_bs_dist, 0.975);
                    lowerq = quantile(irf_bs_dist, 0.025);
                    bs_ci(1,j,i) = 2*irf_est - upperq;
                    bs_ci(2,j,i) = 2*irf_est + lowerq;                
                otherwise 
                    error('Invalid method specified. Use "normal" or "percentile"')
            end
        else
            % standard analytic SE
            if opt==0
                se(:,i)=results.se';
            else
               [EstCov, hacse, coeff]=hac_alt(x(1:end-i+1,:), yy, 'intercept', false, 'smallT', false, 'display', 'off'); 
                se(:,i)=hacse';
            end
        end

    end

    liny(j,:)=reglin(rpos,:); % only take the beta results
    
    if bootstrap == 1
        confidencey = bs_ci;
    else
        sey(j,:)=se(rpos,:);
        confidencey(1,j,:)=liny(j,:)-(sey(j,:)*clevel);
        confidencey(2,j,:)=liny(j,:)+(sey(j,:)*clevel);
    end
end

