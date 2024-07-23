function [liny, confidencey] = VARLP_Avg(data,x,hor,rpos,transformation, clevel, opt, nlags, bias_corrected, bootstrap)
%VARLP_AVG Summary of this function goes here
%   Detailed explanation goes here

    [irf_var, ~] = SVAR(data,x,hor,rpos,transformation, clevel, opt, nlags, bias_corrected);
    [irf_lp, ~] = linlp(data,x,hor,rpos,transformation, clevel, opt, bootstrap);
    
    % S.E. TBD
    var_weight = 0.5
    liny = irf_var * var_weight + irf_lp * (1-var_weight);
    se = abs(0.01 * liny(1, 1));
    
    [dr, dsize] = size(data)
    for j = 1:dsize
        sey(j, :) = se;
        confidencey(1,:,j)=liny(j,:)-(sey(j,:)*clevel);
        confidencey(2,:,j)=liny(j,:)+(sey(j,:)*clevel);
    end
end

