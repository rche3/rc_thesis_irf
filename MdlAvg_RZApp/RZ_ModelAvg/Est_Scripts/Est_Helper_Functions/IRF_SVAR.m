function [irf] = IRF_SVAR(By,ShockVector,nhorizons)
% Auxiliary function for estimating IRFs using VAR coefficients
% returns a K x h matrix of IRFs (K is number of dimensions

nv = size(By,1);
nlags = size(By,3);
irf = zeros(nv,nhorizons);
irf(:,1) = ShockVector; % IRF at horizon 0

% iterate thru horizon 1 to horizon max
for istep = 1:(nhorizons)
    nlagsToUse = min(istep - 1,nlags);
    for ilag = 1:nlagsToUse
        irf(:,istep) = irf(:,istep) + By(:,:,ilag) * irf(:,istep - ilag); % h-step ahead IRF
    end
end

end

