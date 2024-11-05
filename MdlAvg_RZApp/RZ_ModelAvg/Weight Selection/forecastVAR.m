function var_forecast = forecastVAR(train_y, nlag, horizon, fcast_pos)
    % Estimate VAR parameters
    Y = train_y(nlag+1:end, :);
    [Teff, n_vars] = size(Y);

    X = ones(Teff,1);
    for i=1:nlag
        X = [X train_y(nlag+1-i:end-i,:)];
    end
    beta = inv(X'*X)*X'*Y;
    
    % Generate forecasts iteratively
    last_X = [1];
    for i=1:nlag
        last_X = [last_X, train_y(end-(i-1),:)];
    end
    
    total_forecast = zeros(horizon, n_vars);
    % Iterate forward using the VAR(p) structure

    temp_data = train_y;
    TT = size(train_y,1);
    for h = 1:nlag
        XX = [1];
        for p=1:nlag
            XX = [XX, temp_data(end-(p-1),:)];
        end
        temp_data(TT+h,:) = XX * beta;
        total_forecast(h,:) = temp_data(TT+h,:);
    end

    for h = nlag+1:horizon
        XX = [1];
        for p=1:nlag
            XX = [XX, total_forecast(h-(p-1),:)];
        end
        total_forecast(h,:) = XX *  beta;
    end
    final_forecast = total_forecast(horizon,:);
    var_forecast = final_forecast(:,fcast_pos);
end

