function [stateay, stateby, confidenceya, confidenceyb]=statelp_rz(data,x,hor,rpost,transformation, clevel, opt) 

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
        regy(:,i)=results.beta;
        
        if opt==0
            se(:,i)=results.se';
        else
            [EstCov, hacse, coeff]=hac_alt(x(1:end-i+1,:), yy, 'intercept', false, 'smallT', false, 'display', 'off'); 
            se(:,i)=hacse';
        end
    end
    stateay(j,:)=regy(rpost,:);
    stateby(j,:)=regy(rpost+1,:);
    
    seay(j,:)=se(rpost,:);
    seby(j,:)=se(rpost+1,:);
    
    confidenceya(1,:,j)=stateay(j,:)-(seay(j,:)*clevel);
    confidenceya(2,:,j)=stateay(j,:)+(seay(j,:)*clevel);
    confidenceyb(1,:,j)=stateby(j,:)-(seby(j,:)*clevel);
    confidenceyb(2,:,j)=stateby(j,:)+(seby(j,:)*clevel);
end
