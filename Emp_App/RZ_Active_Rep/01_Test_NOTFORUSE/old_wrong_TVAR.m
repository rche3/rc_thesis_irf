    %%%% FUNCTION START %%%%
    
    % setup
    [Ta, nvara] = size(ya);
    [Tb, nvarb] = size(yb);
    K = nvara;
    state = {'a', 'b'};
    num_states = length(state);  % assuming 'state' is already defined
    beta = cell(1, num_states);
    resid = cell(1, num_states);
    sigma_u = cell(1, num_states);
    irfa = zeros(2, hor); % store irfs to two response variables
    irfb = zeros(2, hor);
    irf_comb = zeros(2,2,hor);

    for i=1:2
        % Step 1: Complete the reduced form regression
        yi = ['y' state{i}];
        y = eval(yi);
        ylag = lagmatrix(y,[1:p]);
        ylag = ylag(p+1:end,:);
        TT = size(ylag, 1);
        X = [ones(TT,1) ylag];
        Y = y(p+1:end,:);
        beta_temp = inv(X'*X) * X' * Y;

        % Step 2: Obtain residuals and perform Cholesky decomposition
        resid_temp = Y - X * beta_temp;
        sigma_u_temp = (resid_temp'* resid_temp)/(TT-K*(p+1)+1);
        P = chol(sigma_u_temp, 'Lower'); % P is equal to B^{-1}

        % Step 3: Compute the IRFs recursively
        theta = zeros(hor+1,K,K); % final H x (KxK) matrix of impulse responses
        phi = zeros(hor+1,K,K);
        phi(1,:,:) = eye(K); % VMA coefficients of RF-VAR
        % store the RF-VAR coeffs for easy access
        A = zeros(p,K,K);
        for j=1:p
            Aj_ind = (j-1)*K+2:j*K+1;
            A(j,:,:) = beta_temp(Aj_ind,:)'; % transpose since the beta matrix holds transposed A_1, .. A_p
        end

        % recursive orthogonalised IRFs
        for h=2:hor+1
            phi_temp = zeros(K,K);
            for s = 1:min(h,p)
                A_slice = squeeze(A(s,:,:));
                phi_slice = squeeze(phi(h-s+1,:,:));
                phi_temp = phi_temp + A_slice*phi_slice;
            end            
            phi(h,:,:) = phi_temp;
        end

        for h=1:hor
            phi_slice_h = squeeze(phi(h,:,:));
            theta(h,:,:) = phi_slice_h * P;
        end

        % Step 4: store IRFs
        theta_i = zeros(2,hor); 
        for j=1:2
            for h=1:hor
                theta_i(j,i) = theta(h+1,j,1); % third dimension always 1 since we ordered military first
            end
        end
        % extract the irfs for our two response variables

        % normalise IRFs
        norm = theta(1,1,1);
        theta_i = theta_i / norm;
        irf_comb(i,:,:) = theta_i;

        % (optional) store auxiliary variables
        beta{i} = beta_temp;
        resid{i} = resid_temp;
        sigma_u{i} = sigma_u_temp;
    end

    % retrieve IRFs
    irfa = squeeze(irf_comb(1,:,:));
    irfb = squeeze(irf_comb(2,:,:));


    % placeholders CIs
    cia(1,:,:) = irfa * 1.1;
    cia(2,:,:) = irfa * 0.9;
    cib(1,:,:) = irfb * 1.1;
    cib(2,:,:) = irfb * 0.9;
end

