function W = generateWeightMatrix_v2(num_blocks, block_size, lambda_step)
    % Input validation with defaults
    if nargin < 3
        lambda_step = 0.1;  % Default step size
    end
    if nargin < 2
        block_size = 5;     % Default block size
    end
    
    % Generate lambda values from 0 to 1
    lambda_values = 0:lambda_step:1;
    
    % Create cell array for ndgrid inputs
    grid_inputs = repmat({lambda_values}, 1, num_blocks);
    
    % Use ndgrid with variable number of outputs
    [grid_outputs{1:num_blocks}] = ndgrid(grid_inputs{:});
    
    % Initialize combinations matrix
    combinations = zeros(num_blocks, numel(grid_outputs{1}));
    
    % Reshape each output into a row of combinations
    for i = 1:num_blocks
        combinations(i, :) = grid_outputs{i}(:)';
    end
    
    % Generate final weight matrix
    final_W = zeros(num_blocks * block_size, size(combinations, 2));
    
    % Fill the blocks
    for i = 1:num_blocks
        row_indices = ((i-1) * block_size + 1):(i * block_size);
        i_blockweight = combinations(i, :);
        final_W(row_indices, :) = repmat(i_blockweight, block_size, 1);
    end
    
    W = final_W;
end