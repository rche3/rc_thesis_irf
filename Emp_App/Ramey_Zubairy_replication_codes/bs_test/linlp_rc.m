function [liny,confidencey]=linlp_rc(data,x,hor,rpos,transformation, clevel, opt, bootstrap, nlag) 

[dr,dsize]=size(data);
for j=1:dsize
    [r,nnn]=size(x);
    for i=1:hor
        if transformation==1
            yy=data(i:end,j);
        elseif transformation==2
            yy=(data(i+1:end,j)-data(1:end-i,j))./data(1:end-i,2);
        end

        results=nwest(yy, x(1:end-i+1,:),i);
        reglin(:,i)=results.beta;
        
        % disp(i)
        % disp(yy(30:50))
        % xx = x(1:end-i+1, :);
        % disp(round(xx(30:50, rpos), 4))

        % compute SEs
        if bootstrap == 1
            % bootstrapped SEs
            se(:, i) = compute_bootstrap_se(results.resid, yy, x(1: end-i +1, :), results.beta);
        else
            % standard analytic SE
            if opt==0
                se(:,i)=results.se';
            else
               [EstCov, hacse, coeff]=hac_alt(x(1:end-i+1,:), yy, 'intercept', false, 'smallT', false, 'display', 'off'); 
                se(:,i)=hacse';
            end
        end
        % disp([se, results.se])
        
    end
    liny(j,:)=reglin(rpos,:); % only take the beta results
    sey(j,:)=se(rpos,:);
    disp(se(rpos, :));
    confidencey(1,:,j)=liny(j,:)-(sey(j,:)*clevel);
    confidencey(2,:,j)=liny(j,:)+(sey(j,:)*clevel);
end

