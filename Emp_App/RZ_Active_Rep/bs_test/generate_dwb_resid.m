function innov_vec = generate_dwb_resid(samp_resid,p,wild_settings)
% generates a T x 1 vector of dependent wild bootstrap residuals using the
% Bartlett Kernel.

% inputs:
    % sampled residuals (sampled with replacement)
    % max bandwidth
    % wild settings (parameters for the discrete probability distribution to generate wild errors
   
% override bandwidth manually
[T, ~] = size(samp_resid);

% generate diagonalised residual matrix
U = diag(samp_resid);

% generate Barlett Kernel K weighing matrix
K = zeros(T, T); % this will be Cholesky decomposed into LL'

for i=1:T
    for j=1:T 
        if abs(i-j) >= p+1
            K(i,j) = 0;
        else
            K(i,j) = 1 - abs(i-j)/(p+1);
        end
    end
end
L = chol(K, 'lower');


% generate wild epsilon vector
values = wild_settings(1, :);
probabilities = wild_settings(2,:);
e_vec = disc_rv(values, probabilities, T);
% e_vec = randn(T, 1);

% finish
innov_vec = U * L * e_vec;

end

