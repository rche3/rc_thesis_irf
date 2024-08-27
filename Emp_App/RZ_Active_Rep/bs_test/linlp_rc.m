function [liny,confidencey, bs_beta_means, bs_beta_dist]=linlp_rc(data,x,hor,rpos,transformation, clevel, opt, bootstrap, nlag, nstraps) 

[dr,dsize]=size(data);
if bootstrap == 1
    [~, nvar]=size(x);
    bs_beta_means = nan(dsize, hor);
    bs_ci = nan(2,hor,dsize);
    bs_beta_dist = nan(nstraps, dsize, hor);
else
    % pass
end

for j=1:dsize
    if j == 1
        y_pos_control = 3; % rgov is ordered 3rd in the control vector
    else
        y_pos_control = 2; % rgdp is ordered 2nd in the control vector
    end
    [r,nnn]=size(x);
    for i=1:hor
        disp([j, i]);
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
%             [bs_beta_se, bs_beta_mean] = compute_bootstrap_se1(results.resid, yy, x(1: end-i +1, :), results.beta);
            [bs_beta_se, bs_beta_mean, bs_beta] = compute_dwb_se(results.resid, yy, x(1: end-i +1, :), results.beta, i, nlag, y_pos_control, nstraps);
%             [bs_beta_se, bs_beta_mean] = copute_bootstrap_se(results.resid, yy, x(1: end-i +1, :), results.beta, i, nlag, y_pos_control, wild);
%             [bs_beta_se, bs_beta_mean] = compute_bootstrap_se_plain(results.resid, yy, x(1: end-i +1, :), results.beta);
%             [bs_beta_se, bs_beta_mean] = compute_block_bootstrap_se(results.resid, yy, x(1: end-i +1, :), results.beta, i, nlag, y_pos_control);
%             [bs_beta_se, bs_beta_mean] = compute_block_bootstrap_se1(results.resid, yy, x(1: end-i +1, :), results.beta);
            se(:, i) = bs_beta_se;
            bs_beta_means(j, i) = bs_beta_mean(rpos);
            irf_bs_dist = bs_beta(rpos, :); % returns 1 x B bootstrapped estimates
            bs_beta_dist(:,j,i) = irf_bs_dist';
            upperq = quantile(irf_bs_dist, 0.975);
            lowerq = quantile(irf_bs_dist, 0.025);
            bs_ci(1,i,j) = 2*irf_est - upperq;
            bs_ci(2,i,j) = 2*irf_est + lowerq;
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
    sey(j,:)=se(rpos,:);

    if bootstrap == 1
        confidencey = bs_ci;
    else
        confidencey(1,:,j)=liny(j,:)-(sey(j,:)*clevel);
        confidencey(2,:,j)=liny(j,:)+(sey(j,:)*clevel);
    end
end

