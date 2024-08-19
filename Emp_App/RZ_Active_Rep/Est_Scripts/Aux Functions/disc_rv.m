function X = disc_rv(values,probabilities, num_samples)
    % Generates a random vector of i.i.d. variables from a specified discrete probability distirbution
    %   Detailed explanation goes here
    if numel(values) ~= numel(probabilities)
        error('The number of values must match the number of probabilities.');
    end
    
    if sum(probabilities) ~= 1
        error('Probabilities must sum to 1');
    end
    
    % Create cumulative probabilities
    cum_prob = cumsum(probabilities);
    
    % Generate random numbers
    r = rand(num_samples, 1); % returns an nx1 matrix of random scalars in [0, 1]
    
    % Initialize output
    selected_value = zeros(num_samples, 1);
    
    % Select values based on the random numbers
    for i = 1:num_samples
        val_index = find(r(i) <= cum_prob, 1, 'first'); % this finds the first occurence of a "true" in the array r(i) <= cum_prob
        selected_value(i) = values(val_index); % selects the corresponding  
    end

    X = selected_value;
end
