function res = LP_res(Y,recurShock,respV,nlags,horizon)
% Function for h-step ahead LP

nv = size(Y,2);

% data for LP routine
[y,x,w] = LP_gen_data(Y,recurShock,respV,nlags,horizon);

X = [ones(size(x,1),1), x, w];

 % least-squares LP to generate residuals
[~,~,~,res] = LS(y,X);


end

