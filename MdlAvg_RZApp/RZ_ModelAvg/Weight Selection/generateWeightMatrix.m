function W = generateWeightMatrix()
    lambda_values = 0:0.1:1;
    
    % Use ndgrid to create all combinations
    [l1, l2, l3, l4] = ndgrid(lambda_values);
    
    % Reshape the matrices into column vectors
    l1 = l1(:);
    l2 = l2(:);
    l3 = l3(:);
    l4 = l4(:);
    
    % Combine into a matrix where each row is a combination
    combinations = [l1 l2 l3 l4];
    
    % Transpose to get 4x(11^4) matrix where each column is a combination
    combinations = combinations';
    
    final_W = zeros(size(combinations,1)*5, size(combinations, 2));
    for i = 1:4
        row_indices = ((i-1)*5 + 1):(i*5);
        i_blockweight = combinations(i,:);
        final_W(row_indices,:) = repmat(i_blockweight,5,1);
    end

    W = final_W;
end