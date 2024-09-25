function [stateay, stateby, confidenceya, confidenceyb, bs_beta_dist, pval]=statelp( ...
    data,shock,z,In,hor,transformation, clevel, opt, ...
    bootstrap,method,nlag,nstraps) 
% computes the state dependent local projection, built off Ramey and
% Zubairy (2018) code, but assume no time trends

% settings
h0 = 0; % hypothesis for the pvalues

% dimensions
[dr,dsize]=size(data);
[~, nvar]=size(z);
[T_ind,~] = size(In);

% data setup
constant = ones(T_ind,1);
x=[constant, (1-In).*constant, In.*shock, (1-In).*shock, repmat(In,1,size(z,2)).*z, repmat((1-In),1,size(z,2)).*z];
rposta=3; % position of shock in state A
rpostb=4; % position of shock in state B

% initialise matrices for later use
bs_cia = nan(2,dsize,hor);
bs_cib = nan(2,dsize,hor);
bs_beta_dist = nan(2,nstraps,dsize,hor);
pval = nan(2, dsize, hor);

for j=1:dsize
    y_pos_control = j+1; % since the shock variable is positioned first
    for i=1:hor
        disp([j,i])
        if transformation==1
            yy=data(i:end,j);
         elseif transformation==2
            yy=(data(i+1:end,j)-data(1:end-i,j))./data(1:end-i,2);
        end
        results=nwest(yy, x(1:end-i+1,:),i);
        regy(:,i)=results.beta;
        irf_esta = regy(rposta, i);
        irf_estb = regy(rpostb, i);
        
        if bootstrap == 1
            ctrlstart = rpostb + 1; % the controls start 2 indices after the shock (there are 2 shock terms due to A / B)
            effective_In = In(1:end-i+1,:);
            bs_beta_dist = statelp_dwbse(results.resid, yy, x(1:end-i+1, :), effective_In, results.beta,i,nlag,y_pos_control,ctrlstart, nstraps);
            irf_bs_dista = bs_beta_dist(rposta, :); % returns 1 x B bootstrapped estimates
            irf_bs_distb = bs_beta_dist(rpostb, :); % returns 1 x B bootstrapped estimates

            switch method
                case "normal"
                    bs_se_a = std(irf_bs_dista);
                    bs_se_b = std(irf_bs_distb);
                    % store
                    bs_cia(1,j,i) = irf_esta - clevel * bs_se_a;
                    bs_cia(2,j,i) = irf_esta + clevel * bs_se_a;
                    bs_cib(1,j,i) = irf_estb - clevel * bs_se_b;
                    bs_cib(2,j,i) = irf_estb + clevel * bs_se_b;
                case "percentile"
                    alpha = 1 - normcdf(clevel);
                    uppera = quantile(irf_bs_dista, 1-alpha);
                    lowera = quantile(irf_bs_distb, alpha);
                    upperb = quantile(irf_bs_distb, 1-alpha);
                    lowerb = quantile(irf_bs_distb, alpha);
                    % store 
                    bs_cia(1,j,i) = 2*irf_esta - uppera;
                    bs_cia(2,j,i) = 2*irf_esta + lowera;
                    bs_cib(1,j,i) = 2*irf_estb - upperb;
                    bs_cib(2,j,i) = 2*irf_estb + lowerb;
                otherwise
                    error('Invalid method specified. Use "normal" or "percentile"')
            end        
            % bs_pval takes in the bootstrapped estimate distribution ad
            % creates pval
            pval(1,j,i) = bs_pval(irf_bs_dista, irf_esta, h0); 
            pval(2,j,i) = bs_pval(irf_bs_distb, irf_estb, h0);

            % store everything
            bs_beta_dist(1,:,j,i) = irf_bs_dista';
            bs_beta_dist(2,:,j,i) = irf_bs_distb';
        else
            if opt==0
                se(:,i)=results.se';
            else
                [EstCov, hacse, coeff]=hac_alt(x(1:end-i+1,:), yy, 'intercept', false, 'smallT', false, 'display', 'off'); 
                se(:,i)=hacse';
            end
            
            irfse_a = results.se(rposta);
            irfse_b = results.se(rpostb);
            % compute p values
            tau_a = (irf_esta - h0) / irfse_a;
            tau_b = (irf_estb - h0) / irfse_b;
            pval(1,j,i) = 2*(1-normcdf(abs(tau_a)));
            pval(2,j,i) = 2*(1-normcdf(abs(tau_b)));
        end
    end

    % retrieve point estimates
    stateay(j,:)=regy(rposta,:);
    stateby(j,:)=regy(rpostb,:);
    
    % store CIs
    if bootstrap == 1
        confidenceya = bs_cia;
        confidenceyb = bs_cib;        
    else
        seay(j,:)=se(rposta,:);
        seby(j,:)=se(rpostb,:);
        confidenceya(1,j,:)=stateay(j,:)-(seay(j,:)*clevel);
        confidenceya(2,j,:)=stateay(j,:)+(seay(j,:)*clevel);
        confidenceyb(1,j,:)=stateby(j,:)-(seby(j,:)*clevel);
        confidenceyb(2,j,:)=stateby(j,:)+(seby(j,:)*clevel);
    end
    
end