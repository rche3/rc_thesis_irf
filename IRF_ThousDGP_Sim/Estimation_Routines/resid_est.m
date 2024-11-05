function [IRF, n_lags_est] = resid_est(data_sim,settings)
%RESID_EST Summary of this function goes here
%   Detailed explanation goes here

run('Estimation_Setup'); % common setup for all estimation methods

%%% compute the LP residuals (LP errors) at max horizon H
res = IRF_LP_res(Y,recursiveShock,nlags,IRF_hor - 1); % residual matrix for LP-form of VAR(p)
% disp(Y)
% disp(res)
% disp(size(res))
% (we subtract 1 because our regression above only includes non-IV dependent variables) 
theta_1_tH = res(:, normalizeV-1); % this is theta_1 per BL 2022
theta_i_tH = res(:, responseV-1); % this is theta_i, residuals for the response variable

%%% compute the estimator using the proxy z_t
% restrict the z_t to p+1:T
z_t = Y(nlags+1:end,1)'; % the proxy z_t is stored in the first column, but we must remove nlags to make it consistent size
% disp(size(z_t))

T_effective = length(z_t) - (IRF_hor - 1); % because the 'real' H is actually IRF_hor - 1 since no zero-indexing in matlab
% T_effective is mathematically T-H-p because that is our available sample
% for constructing resid estimator

% Replace fliplr(hankel()) with explicit matrix construction
z_t_matrix = zeros(T_effective, IRF_hor);
for i = 1:T_effective
    for j = 1:IRF_hor
        z_t_matrix(i,j) = z_t(i + IRF_hor - j);
    end
end

% Replace bsxfun multiplication and sum with explicit loop
numerator = zeros(1, IRF_hor);
for i = 1:T_effective
    htemp_irfs = theta_i_tH(i) * z_t_matrix(i,:);
    numerator = numerator + htemp_irfs;
end

% for t = 1:IRF_hor
%     sum_product = 0;
%     for i = 1:T_effective
%         sum_product = sum_product + theta_i_tH(i) * z_t_matrix(i,t);
%     end
%     numerator(t) = sum_product;
% end

% z_t_matrix = fliplr(hankel(z_t(1:T_effective), z_t(T_effective:end)));
% 
% % multiply each row of z_t_matrix by the corresponding scalar resid
% numerator = sum(bsxfun(@times, theta_i_tH, z_t_matrix), 1);
% % ^bsxfun is binary singleton expansion function - applies
% % element-by-element operation to two arrays, in this case, for each "t",
% % multiplies the line of the z_t_matrix by the theta_i_tH scalar

denominator = theta_1_tH' * z_t(IRF_hor:end)';

IRF = (1/denominator) * numerator;
IRF = IRF';

end


