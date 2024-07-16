function [Bc,Br,Bx,By,Sigma,Sxx,w] = LP(y, x, w,recurShock,respV,nlags,horizon)
% Function for h-step ahead LP

% manually generate the response variable series
Y = [x y];
nv = size(Y,2);

% since we are just normalising with respect to the shock, TBD how to
% impact

% data for LP routine
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
Bx = Beta(2); % impulse variable, effect of 
By = Beta((recurShock + 2):end); % lagged controls
By = reshape(By,[nv,nlags]);

end

