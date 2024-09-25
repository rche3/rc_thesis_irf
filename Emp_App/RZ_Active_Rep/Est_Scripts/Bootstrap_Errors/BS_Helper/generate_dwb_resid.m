function innov_vec = generate_dwb_resid(resid_matrix,bw,wild_settings)
% generates a T x K vector of dependent wild bootstrap residuals using the
% Bartlett Kernel.

% inputs:
    % sampled residuals (sampled with replacement), T x K
    % max bandwidth
    % wild settings (parameters for the discrete probability distribution to generate wild errors
   
% override bandwidth manually
[T, nvar] = size(resid_matrix);

% generate wild epsilon vector
values = wild_settings(1, :);
probabilities = wild_settings(2,:);
e_vec = disc_rv(values, probabilities, T);

% storage for end
innov_vec = zeros(T, nvar);

for n=1:nvar
    % generate diagonalised residual matrix
    samp_resid = resid_matrix(:,n);
    U = diag(samp_resid);
    
    % generate Barlett Kernel K weighing matrix
    K = zeros(T, T); % this will be Cholesky decomposed into LL'
    
    for i=1:T
        for j=1:T 
            if abs(i-j) >= bw+1
                K(i,j) = 0;
            else
                K(i,j) = 1 - abs(i-j)/(bw+1);
            end
        end
    end
    L = chol(K, 'lower');
    
    % finish
    innov_vec(:,n) = U * L * e_vec;
end

