%----------------------------------------------------------------
% File/Folder Names
%----------------------------------------------------------------

exper_filename = exper_files{ne}; % Name of current experiment
exper_plotname = exper_names{ne};
file_name = fullfile(mode_folders{n_mode}, lags_folders{nf}, exper_filename); % Name of .mat results file
folder_name = fullfile(mode_folders{n_mode}, lags_folders{nf}, exper_folders{ne}); % Name of figure folder
output_folder = fullfile(output_dir, folder_name); % Name of output folder  

%----------------------------------------------------------------
% Load Simulation Results
%----------------------------------------------------------------

% see if this exper is start of group
if (ne == 1) || (exper_group_end(ne-1) == 1)
    res = load(fullfile(rootfolder, strcat(file_name, '.mat'))); % Directly load
    mkdir(output_folder); % Create output folder
else
    res_part = load(fullfile(rootfolder, strcat(file_name, '.mat'))); % Load
    res = combine_struct(res, res_part, [], []); % Merge
end

horzs = res.settings.est.IRF_select; % Impulse response horizons
methods_names_plot = methods_names{ne};