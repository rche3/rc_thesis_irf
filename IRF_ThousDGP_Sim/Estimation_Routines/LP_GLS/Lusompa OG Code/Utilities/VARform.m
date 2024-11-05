
function [y,x] = VARform(data,lags,deterministic_terms)

% Builds vectors of dependent variables and (independent) lags for a VAR.
% Usage:   [y,x] = VARform(data,lags,deterministic_terms)
% Inputs:  data:   Matrix of data (assumed: # periods > # variables)
%          lags:   # lags to be included in the VAR
%          determ: (opt.) Deterministic regressors: 0=none, 1=const,
%                  2=const&trend; 3=const&trend^2 (default is 1)
% Outputs: y:      dependent variable
%          x:      independent variable (in the form [determ. terms,
%                  data(lag 1), data(lag 2), ... data(lag L)], where
%                  determ. terms = [constant trend trend^2] (iff included)
%
% by B. Kolb, Jan. 2015

%% Housekeeping
if size(data,2) > size(data,1)
    data = data';
end

[T N] = size(data);

%% Create matrices
% Dependent variables vector
y = data(lags+1:T,:);

% Independent variables vector/matrix
x = zeros(T-lags,N*lags);
for ii=1:lags
    x(:,(ii-1)*N+1:ii*N) = data(lags-ii+1:T-ii,:);
end

% add deterministic terms to regressors
switch deterministic_terms
    case 0
        return
    case 1 % constant
        const = ones(T-lags,1);
        x = [const x];
    case 2 % time trend and constant
        const = ones(T-lags,1);
        trend = (1:size(x,1))';
        x = [const trend x];
    case 3 % squared time trend, linear time trend, and constant
        const = ones(T-lags,1);
        trend = (1:size(x,1))';
        x = [const trend trend.^2 x];
    otherwise % not valid
        disp('Please select between 0 and 3 deterministic terms')
end
