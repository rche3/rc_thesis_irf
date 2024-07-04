function res = IRF_LP_res(Y,recurShock,respV,nlags, nhorizons)
% Computes a matrix of residuals K x T from which the relevant ones can be
% selected for IRF

nT = size(Y,1);

% error checking
if (nhorizons + nlags) >= nT
    error("Number of horizons too large! No obs in sample!")
end

% compute the residuals for the max horizon H
res = LP_res(Y,recurShock,respV,nlags,nhorizons); 

end

