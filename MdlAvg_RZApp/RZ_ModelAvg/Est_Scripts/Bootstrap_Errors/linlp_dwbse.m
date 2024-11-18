function beta_bs = linlp_dwbse( ...
    resid, y, x, beta, h, nlag, y_pos_control, ctrl_start, B)
% computes dependent wild bootstrap for LINEAR LP
% Inputs:
    % resid is original residual T x k computed from the estimated coefficients
    % y is lp response variable
    % x is lp regressor matrix
    % beta is estimated beta
    % h is desired lp horizon 
    % nlag is number of lagged obs for lp controls
    % ctrl_start is the position in x where controls begin
    % B is number of bootstrap replications
% Outputs:
    % beta_bs is size(x,1) x B empirical bootstrap distribution

% check the size is the same for resid, y, x
if size(resid, 1) == size(y, 1) && size(x, 1) == size(y,1)
    % pass
else
    disp('Error, dimensions of resid, regressors, and dep. var not equal')
end

% bootstrap settings
% setup
[TT,k] = size(resid); % by design, this should be T - h + 1 in length, h \in [1,20]
beta_bs = nan(size(beta, 1), B);
T = TT;
[xrow, xcol] = size(x);

% DWB bootstrap settings
e1 = 1; e2 = -1;
p1 = 0.5; p2 = 1-p1;
values = [e1, e2];
probabilities = [p1, p2];
wild_settings = [values; probabilities];

% create residual bootstrap samples

for b = 1:B
    % disp(['Bootstrap rep: ', num2str(b), ' out of ', num2str(B)])

    % create the bootstrap dependent variable y and x
    y_bs = zeros(T,1); % this is the storage for our bootstrap LP response var
    x_bs = x; % iteratively updated bootstrap LP regressor matrix

    % create the DWB residuals
    dwb_innov = generate_dwb_resid(resid, nlag, wild_settings);
    
    for i = 1:T
        % compute and store y_bs variables
        y_bs(i) = x_bs(i,:) * beta + dwb_innov(i);
        
        % isolate the control vector
        static = x_bs(:, 1:ctrl_start-1);
        ctrls = x_bs(:, ctrl_start:end);
        
        % compute the "unsplit" control vector - T x kp
        z_temp = ctrls;

        % insert the current y value as the value of the:
            % 1st lag of response variable - h+1 periods ahead of now
            % 2nd lag of response variable - h+2 periods ahead of now
            % ...
            % pth lag of response variable - h+p periods ahead of now
                 
        for pi = 1:nlag
            y_insert_idx = y_pos_control + (pi-1)*k;
            future_ctrl_idx = i+h+(pi-1);
            if future_ctrl_idx <= T
                z_temp(future_ctrl_idx,y_insert_idx);
            else
                % pass
            end
            % note that it is not "h+1" since RZ uses h=1 as contemp
            % note (pi-1) since we add +1 to the (h+1) for how much (pi-1)
        end

        % recombine static and recurisvely computed parts of regressor matrix
        x_bs = [static z_temp];
    end

    % estimate bootstrapped beta
    results=nwest_rc(y_bs, x_bs, 0); % note we don't need to lag x as input arg x is already lagged
    beta_bs(:, b) = results.beta;
end

%function end
end

