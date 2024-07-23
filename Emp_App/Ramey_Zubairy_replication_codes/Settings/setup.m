line_colors = {
    [0, 0, 0],       % Black
    [1, 0, 0],       % Red
    [0, 0, 1],       % Blue
    [0, 0.5, 0],     % Dark Green
    [0.5, 0, 0.5],   % Purple
    [1, 0.5, 0],     % Orange
    [0, 0.75, 0.75], % Teal
    [0.75, 0, 0.75]  % Magenta
};

% Predefined colors for patches (lighter versions)
patch_colors = {
    [0.7, 0.7, 0.7], % Light Gray
    [1, 0.8, 0.8],   % Light Red
    [0.8, 0.8, 1],   % Light Blue
    [0.8, 1, 0.8],   % Light Green
    [1, 0.8, 1],     % Light Purple
    [1, 0.9, 0.8],   % Light Orange
    [0.8, 1, 1],     % Light Teal
    [1, 0.8, 1]      % Light Magenta
};

plotStyles = {
    {'b-', 'LineWidth', 1.2},     % Solid blue line
    {'r--', 'LineWidth', 1.2},    % Dashed red line
    {'g-.', 'LineWidth', 1.2},    % Dash-dot green line
    {'m:', 'LineWidth', 1.2},     % Dotted magenta line
    {'k-', 'LineWidth', 1.5},     % Thick solid black line
    {'c--', 'LineWidth', 1.0},    % Thin dashed cyan line
    {'y-.', 'LineWidth', 1.3},    % Medium dash-dot yellow line
    {'[0.5 0.5 0.5]', 'LineStyle', '-', 'LineWidth', 1.2},  % Solid gray line
    {'[0.8 0.4 0]', 'LineStyle', '--', 'LineWidth', 1.2},   % Dashed orange line
    {'[0 0.5 0]', 'LineStyle', ':', 'LineWidth', 1.4}       % Thick dotted dark green line
};

set(0, 'DefaultTextInterpreter', 'none');
set(0, 'DefaultLegendInterpreter', 'none');
set(0, 'DefaultAxesTickLabelInterpreter', 'none');
set(0, 'DefaultColorbarTickLabelInterpreter', 'none');