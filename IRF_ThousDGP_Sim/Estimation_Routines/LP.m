function [Bc,Br,Bx,By,Sigma,Sxx,w] = LP(Y,recurShock,respV,nlags,horizon)
% Function for h-step ahead LP

nv = size(Y,2);

% data for LP routine
% disp('This is the Y data for LP_gen_data properly via LP_est')
% disp(size(Y))
% disp(Y(1:10,:))
[y,x,w] = LP_gen_data(Y,recurShock,respV,nlags,horizon);
% if horizon==0
%     disp(size(y))
%     disp(y(1:10,:))
% else
%     %pass
% end
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

