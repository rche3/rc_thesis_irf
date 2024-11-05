function res = IRF_LP_res(Y,recurShock,nlags, nhorizons)
% Computes a matrix of residuals (T-H-p) x K from estimating the LP form of
% VAR(p) model.

nT = size(Y,1);

% error checking
if (nhorizons + nlags) >= nT
    error("Number of horizons too large! No obs in sample!")
end

% compute the residuals matrix for max horizon, return (K x (T-H-p))
res = LP_res(Y,recurShock,nlags,nhorizons); 

end

