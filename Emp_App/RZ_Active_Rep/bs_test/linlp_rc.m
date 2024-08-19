function [liny,confidencey, bs_beta_means]=linlp_rc(data,x,hor,rpos,transformation, clevel, opt, bootstrap, nlag, wild) 

[dr,dsize]=size(data);
bs_beta_means = nan(dsize, hor);
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

        % compute SEs
        if bootstrap == 1
%             [bs_beta_se, bs_beta_mean] = compute_bootstrap_se1(results.resid, yy, x(1: end-i +1, :), results.beta);
            [bs_beta_se, bs_beta_mean] = compute_dwb_se(results.resid, yy, x(1: end-i +1, :), results.beta, i, nlag, y_pos_control, wild);
%             [bs_beta_se, bs_beta_mean] = copute_bootstrap_se(results.resid, yy, x(1: end-i +1, :), results.beta, i, nlag, y_pos_control, wild);
%             [bs_beta_se, bs_beta_mean] = compute_bootstrap_se_plain(results.resid, yy, x(1: end-i +1, :), results.beta);
%             [bs_beta_se, bs_beta_mean] = compute_block_bootstrap_se(results.resid, yy, x(1: end-i +1, :), results.beta, i, nlag, y_pos_control);
%             [bs_beta_se, bs_beta_mean] = compute_block_bootstrap_se1(results.resid, yy, x(1: end-i +1, :), results.beta);

            se(:, i) = bs_beta_se;
            bs_beta_means(j, i) = bs_beta_mean(rpos);
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
    confidencey(1,:,j)=liny(j,:)-(sey(j,:)*clevel);
    confidencey(2,:,j)=liny(j,:)+(sey(j,:)*clevel);
end

