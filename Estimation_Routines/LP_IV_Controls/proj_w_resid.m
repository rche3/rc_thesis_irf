function resid = proj_w_resid(input_vec_position, Y, inp_rel_cont, IV_location, nlags, h, include_proxy)
% input_vec_position defines the column of Y (the data matrix) that is going to be projected on the control
% include_proxy, 0 = only set the lagged y_{t-1}... y_{t-p} as control, 1 = include "future" z_t
% inp_rel_cont denotes the time subscript on the input vector relative to the control - inp_rel_cont = i means we are projecting y_{t+i} on W_t,


% To get the residuals, we need the input vector, control vector (for proj matrix)
input_vector = Y(:, input_vec_position);
T = size(input_vector, 1);
y = Y(:, IV_location+1:end);
z = Y(:, IV_location);

% create lagged y's
% lagged_y_col = zeros(T-h-nlags, 1+K*nlags+h)

% for i:nlags
lagged_y = lagmatrix(y, 1:nlags);

if include_proxy == 1
    if h==0
        W = [ones(T-h-nlags,1), lagged_y(nlags+1:end-h, :)];
    else 
        lagged_z = lagmatrix(z, 0:h-1);
        W = [ones(T-h-nlags,1), lagged_y(nlags+1:end-h, :), lagged_z(nlags+h+1:end, :)];
    end
else
    W = [ones(T-h-nlags,1), lagged_y(nlags+1:end-h, :)];
end
% ^ column of lagged y's stacked as [y_{t-1}; y_{t-2}; ... y_{t-p}]
input_vector = input_vector(nlags+1+inp_rel_cont:end-(h-inp_rel_cont));

resid = input_vector - W * inv(W' * W) * W' * input_vector;

% figure;
% hold on
% plot(input_vector)
% plot(resid)
% legend('input vec', 'resid')

end

