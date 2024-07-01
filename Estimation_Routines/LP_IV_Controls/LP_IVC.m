function [Bc,Br,Bx,By,Sigma,Sxx,w] = LP_IVC(Y,recurShock,respV,nlags,horizon)

% Function for computing the Stock, Watson 2018 LP-IV projection with
% extral controls

nv = size(Y,2);

%%% Below is a work-in-progress
%{ 
%%% first, getting the control-projection residual datapoints

% generate the control vector

y_lags = Y(: end)
W_t = [1, y_lags];

% define projection matrix
P_w = W_t * inv(W_t' * W_t) * W_t';

% create the residuals for z_t and y_t


%%% end work-in-progress section, commented out for now
%}


% Function for h-step ahead LP

% data for LP routine
[y,x,w] = LP_gen_data_controls_yz(Y,recurShock,respV,nlags,horizon);
X = [ones(size(x,1),1), x, w];

 % least-squares LP
[Beta,Sigma,Sxx,~] = LS(y,X);

% store LP coefficients
Bc = Beta(1); % constant
if recurShock == 1 % check if there are contemperaneous controls
    Br = []; % contemperaneous control
else
    Br = Beta(3:(recurShock+1));
end
Bx = Beta(2); % impulse variable
By = Beta((recurShock + 2):end); % lagged controls
By = reshape(By,[nv,nlags]);

end

