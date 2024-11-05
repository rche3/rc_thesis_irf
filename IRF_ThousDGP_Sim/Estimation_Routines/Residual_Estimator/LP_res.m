function res = LP_res(Y,recurShock,nlags,h)
% Function for h-step ahead LP

nv = size(Y,2);

% data for LP routine
[y, x] = LP_gen_data_res(Y,recurShock,nlags,h);

X = [ones(size(x,1),1), x];

 % least-squares LP to generate residuals
[~,~,~,res] = LS(y,X);

% figure;
% hold on
% plot(res)
% plot(y)
% legend('residuals', 'y')

end

