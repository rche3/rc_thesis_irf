function res = IRF_LP_res(Y,recurShock,respV,nlags,nhorizons)
% Auxiliary function for estimating IRFs using least-squares LP

nT = size(Y,1);

% error checking
if (nhorizons + nlags) >= nT
    error("Number of horizons too large! No obs in sample!")
end

res = zeros(1,nhorizons + 1);

% go thru horizon 0 to horizon max
for h = 0:nhorizons
    res = LP_res(Y,recurShock,respV,nlags,h); % this produces the residuals per horizon. 
    % res should be a (T-p-H) x 1 vector
end

end

