function [stateay, stateby, confidenceya, confidenceyb, bs_beta_dist, pval]=statelp_rc(data,x,hor,rpost,transformation, clevel, opt, bootstrap, emp, nlag, nstraps) 

% settings
h0 = 0; % hypothesis for the pvalues

% constants
[dr,dsize]=size(data);
[~, nvar]=size(x);

% initialise matrices for later use
bs_cia = nan(2,hor,dsize);
bs_cib = nan(2,hor,dsize);
bs_beta_dista = nan(nstraps, dsize, hor);
bs_beta_distb = nan(nstraps, dsize, hor);
pval = nan(2, dsize, hor);


for j=1:dsize
    [r,nnn]=size(x);
    if j == 1
        y_pos_control = 3; % rgov is ordered 3rd in the control vector
    else
        y_pos_control = 2; % rgdp is ordered 2nd in the control vector
    end
    
    for i=1:hor
%         disp([j,i])
        if transformation==1
            yy=data(i:end,j);
         elseif transformation==2
            yy=(data(i+1:end,j)-data(1:end-i,j))./data(1:end-i,2);
        end
        results=nwest(yy, x(1:end-i+1,:),i);
        regy(:,i)=results.beta;
        irf_esta = regy(rpost, i);
        irf_estb = regy(rpost+1, i);
        
        if bootstrap == 1
            [bs_beta_se, bs_beta_mean, bs_beta_dist] = compute_dwb_se(results.resid, yy, x(1:end-i+1, :), results.beta, i, nlag, y_pos_control, nstraps);
            irf_bs_dista = bs_beta(rpost, :); % returns 1 x B bootstrapped estimates
            irf_bs_distb = bs_beta(rpost+1, :); % returns 1 x B bootstrapped estimates

            if emp == 1
                % state a
                uppera = quantile(irf_bs_dista, 0.975);
                lowera = quantile(irf_bs_distb, 0.025);
                % state b
                upperb = quantile(irf_bs_distb, 0.975);
                lowerb = quantile(irf_bs_distb, 0.025);
                % store 
                bs_cia(1,i,j) = 2*irf_esta - uppera;
                bs_cia(2,i,j) = 2*irf_esta + lowera;
                bs_cib(1,i,j) = 2*irf_estb - upperb;
                bs_cib(2,i,j) = 2*irf_estb + lowerb;
            else
                bs_se_a = var(irf_bs_dista)^(1/2);
                bs_se_b = var(irf_bs_distb)^(1/2);
                % store
                bs_cia(1,i,j) = irf_esta - clevel * bs_se_a;
                bs_cia(2,i,j) = irf_esta + clevel * bs_se_a;
                bs_cib(1,i,j) = irf_estb - clevel * bs_se_b;
                bs_cib(2,i,j) = irf_estb + clevel * bs_se_b;
            end

            % computes tau_b = (beta_b - \hat{beta})/SE_b(beta_b)
            % bs_pval takes in the bootstrapped estimate distribution, 
            pval(1,j,i) = bs_pval(bs_beta_dista, irfesta, h0); 
            pval(2,j,i) = bs_pval(bs_beta_distb, irfestb, h0);

            % store everything
            bs_beta_dista(:,j,i) = irf_bs_dista';
            bs_beta_distb(:,j,i) = irf_bs_distb';

        else
            if opt==0
                se(:,i)=results.se';
            else
                [EstCov, hacse, coeff]=hac_alt(x(1:end-i+1,:), yy, 'intercept', false, 'smallT', false, 'display', 'off'); 
                se(:,i)=hacse';
            end
            
            irfse_a = results.se(rpost);
            irfse_b = results.se(rpost+1);
            % compute p values
            tau_a = (irf_esta - h0) / irfse_a;
            tau_b = (irf_estb - h0) / irfse_b;
            pval(1,j,i) = 2*(1-normcdf(abs(tau_a)));
            pval(2,j,i) = 2*(1-normcdf(abs(tau_b)));
        end
    end

    % retrieve point estimates
    stateay(j,:)=regy(rpost,:);
    stateby(j,:)=regy(rpost+1,:);
    
    if bootstrap == 1
        confidenceya = bs_cia;
        confidenceyb = bs_cib;        
    else
        seay(j,:)=se(rpost,:);
        seby(j,:)=se(rpost+1,:);
        confidenceya(1,:,j)=stateay(j,:)-(seay(j,:)*clevel);
        confidenceya(2,:,j)=stateay(j,:)+(seay(j,:)*clevel);
        confidenceyb(1,:,j)=stateby(j,:)-(seby(j,:)*clevel);
        confidenceyb(2,:,j)=stateby(j,:)+(seby(j,:)*clevel);
    end
    
end

bs_beta_dist(1,:,:,:) = bs_beta_dista;
bs_beta_dist(2,:,:,:) = bs_beta_distb;
