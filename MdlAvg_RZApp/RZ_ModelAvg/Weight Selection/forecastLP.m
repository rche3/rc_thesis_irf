function lp_forecast = forecastLP(train_y, nlag, horizon, fcast_pos)        
    % Prepare dependent variable (y_{t+h})
    y_h = train_y(horizon+nlag:end, fcast_pos);
    % y_h is (T-h-nlag) x 1
    Teff = size(y_h,1);

    % Prepare regressors (y_t)
    X = ones(Teff, 1);
    for p = 1:nlag
        X = [X train_y(nlag-p+1:end-horizon-(p-1), :)];
    end
    
    % OLS estimation
    beta = inv(X'*X)*X'*y_h;
    
    % Generate forecast
    last_X = [1];
    for p=1:nlag
        last_X = [last_X train_y(end-(p-1),:)];
    end
    lp_forecast = last_X * beta;
end

