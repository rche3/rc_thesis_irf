% Start the timer
tic;

estimand_vars = {'Recursive'};
% estimand_vars = {'IV'};
% Define the values for the DGP types
dgp_types = {'G','MP'};

% Iterate over all combinations of estimand variables and DGP types
for i = 1:length(estimand_vars)
    for j = 1:length(dgp_types)
        % Set the current estimand variable and DGP type
        estimand = estimand_vars{i};
        dgp = dgp_types{j};
        
        % Display the current combination
        disp(['Running run_dfm.m with estimand: ' estimand ' and DGP: ' dgp]);
        
        % Assign the current estimand and dgp to the base workspace
        assignin('base', 'estimand_type', estimand);
        assignin('base', 'dgp_type', dgp);
        
        % Run your script (run_dfm.m)
        try
            run('run_dfm.m');
        catch ME
            disp(['Error running run_dfm.m: ' ME.message]);
            % Optionally, you can display more error details:
            % disp(getReport(ME));
        end
        
        % clearvars
        % Add any necessary post-processing or saving of results here
        
        % Pause for a short duration to avoid overwhelming the system (optional)
        pause(1); % Pause for 1 second
    end
end

% Stop the timer and display the total time with variable details
total_time = toc;
fprintf('Total time to run run_multiple_dfm with:\n');
fprintf('Estimand variables: %s\n', strjoin(estimand_vars, ', '));
fprintf('DGP types: %s\n', strjoin(dgp_types, ', '));
fprintf('Total time: %.2f seconds\n', total_time);
