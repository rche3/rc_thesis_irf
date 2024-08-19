function [liny,confidencey]=linlp_biascorrect(data,x,hor,rpos,transformation, clevel, opt) 

% bias-corrected LP, code from Li et. al (2024+)
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
        if opt==0
        se(:,i)=results.se';
        else
       [EstCov, hacse, coeff]=hac_alt(x(1:end-i+1,:), yy, 'intercept', false, 'smallT', false, 'display', 'off'); 
        se(:,i)=hacse';
        end
        
    end
    liny(j,:)=reglin(rpos,:);% random placeholder for testing
    
    % do bias correction
    IRF = liny(j, :);
    w = x(:, rpos+1:end); % define control matrix as everything ordered after the shock in the "x" matrix of regressor
    
    liny(j, :) = bias_correction(IRF, w);

    sey(j,:)=se(rpos,:);
    confidencey(1,:,j)=liny(j,:)-(sey(j,:)*clevel);
    confidencey(2,:,j)=liny(j,:)+(sey(j,:)*clevel);
end

end
