function [y, x] = LP_gen_data_res(Y,recurShock,nlags,h)
% LP_GEN_DATA_RES generates data (Y and X) to estimate the LP form of
% VAR(p)

nv = size(Y,2);
nT = size(Y,1);

y  = Y(:, recurShock+1:end); % everything except the recursiveShock (ordered 1st)
if recurShock > 1 % check if there are contemperaneous controls
    w  = [ lagmatrix(Y(:,1:(recurShock - 1)), h) , lagmatrix( Y , (1:nlags) + h ) ]; % control variables (contemporaneous vars, lagged vars, no constant)
else
    x  = lagmatrix(y, (1:nlags) + h );
end
y = y((nlags + h + 1):end, :);
x = x((nlags + h + 1):end, :); % Warning: here w include both contemperaneous and lagged controls
end
