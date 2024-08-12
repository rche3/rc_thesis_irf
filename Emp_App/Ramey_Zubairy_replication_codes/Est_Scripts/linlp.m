function [liny,confidencey]=linlp(data,x,hor,rpos,transformation, clevel, opt, nlag, bootstrap) 

% code based on Ramey and Zubairy original - edited to allow for
% bootstrapped SEs
[dr,dsize]=size(data);
for j=1:dsize
    if j == 1
        y_pos_control = 3; % rgov is ordered 3rd in the control vector
    else
        y_pos_control = 2; % rgdp is ordered 2nd in the control vector
    end
    [r,nnn]=size(x);
    for i=1:hor
        if transformation==1
            yy=data(i:end,j);
        elseif transformation==2
            yy=(data(i+1:end,j)-data(1:end-i,j))./data(1:end-i,2);
        end

        results=nwest(yy, x(1:end-i+1,:),i);
        reglin(:,i)=results.beta;

        % compute SEs
        if bootstrap == 1
            % bootstrapped SEs
            [bs_beta_se, ~] = lp_bserrors(results.resid, yy, x(1: end-i +1, :), results.beta, i, nlag, y_pos_control);
            se(:, i) = bs_beta_se;
        else
            % standard analytic SE,
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

