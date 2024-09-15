function [liny,confidencey]=linlp_lagaug(data,x,hor,rpos,transformation, clevel, opt, bootstrap) 

% for lag-aug just need to change the regressor to include extra lag.
[dr,dsize]=size(data);
for j=1:dsize % looping over the response parameters
    [r,nnn]=size(x);
    
    for i=1:hor % looping over the horizons
        if transformation==1
            yy=data(i:end,j);
        elseif transformation==2
            yy=(data(i+1:end,j)-data(1:end-i,j))./data(1:end-i,2);
        end

        results=nwest(yy, x(1:end-i+1,:),i);
        reglin(:,i)=results.beta;
        if bootstrap == 1
            % bootstrapped SEs
            se = lplagaug_bserrors(results.resid, yy, x(1: end-i +1, :), results.beta);
        else
            % standard analytic SE
            if opt==0
                se(:,i)=results.se';
            else
               [~, hacse, ~]=hac_alt(x(1:end-i+1,:), yy, 'intercept', false, 'smallT', false, 'display', 'off'); 
                se(:,i)=hacse';
            end
        end
        
    end

    
    liny(j,:)=reglin(rpos,:);% random placeholder for testing
    sey(j,:)=se(rpos,:);
    confidencey(1,:,j)=liny(j,:)-(sey(j,:)*clevel);
    confidencey(2,:,j)=liny(j,:)+(sey(j,:)*clevel);
end

