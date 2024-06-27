function res = IRF_LP_res(Y,recurShock,respV,nlags,nhorizons)
% Auxiliary function for estimating IRFs using least-squares LP

nT = size(Y,1);

% error checking
if (nhorizons + nlags) >= nT
    error("Number of horizons too large! No obs in sample!")
end

% compute the residuals for the max horizon H
res = LP_res(Y,recurShock,respV,nlags,nhorizons); % this produces the residuals per horizon

end

