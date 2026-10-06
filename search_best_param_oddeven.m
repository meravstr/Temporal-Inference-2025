function [best_param, err_mean_curve, err_std_curve] = search_best_param_oddeven(y, method, gamma, param_range)

% SEARCH_BEST_PARAM_ODDEVEN chooses the free parameter to use 
% with the requested deconvolution method for fluorescence-only recordings. 
% It follows Jewell & Witten 2018.
% The algorithm splits each trace into training and testing sets by assigning odd- and even-indexed time points to one set or the other. 
% It deconvolves one set, reconvolves it to a predicted calcium trace, and compares it to the other set's fluorescence. 
% It repeats this process for each tested parameter. 
% The chosen parameter value is the largest one that still yields an error below the minimal error plus the standard deviation at the minimal error.
% 
% INPUTS: y - fluorescence traces T x N  
% gamma - the calcium decay per time bin; 
% param_range - a vector of candidate parameter values, 
% 
% OUTPUTS best_param - the chosen value; 
% err_curve and err_std_curve - the mean and standard deviation cross-validation error at each param_range value 

n_t = size(y,1);
n_p = size(y,2);
n_pairs = floor(n_t/2);
sample_pairs = floor(n_pairs/3); % timepoints from the end to compare inferred calcium and fluorescence, for methods without beta0 which corrects for the shifts from the very beginning
y_odd  = y(1:2:n_pairs*2,:);
y_odd_part_nodc = y_odd(end-sample_pairs:end,:)-repmat(mean(y_odd(end-sample_pairs:end,:),1),sample_pairs+1,1);
y_even = y(2:2:n_pairs*2,:);
y_even_part_nodc = y_even(end-sample_pairs-1:end-1,:)-repmat(mean(y_even(end-sample_pairs-1:end-1,:),1),sample_pairs+1,1);
gamma2 = gamma^2;   % odd/even split halves the effective sampling rate

n_params = numel(param_range);
err_mean_curve = nan(1, n_params);
err_std_curve = nan(1, n_params);

for i = 1:n_params
    p = param_range(i);

    [r_odd,  start_odd,r0_odd,beta0_odd]  = run_method(y_odd,  gamma2, method, p);
    [r_even, start_even,r0_even,beta0_even] = run_method(y_even, gamma2, method, p);
    c_odd  = rebuild_calcium(r_odd, r0_odd,  start_odd, beta0_odd, gamma2);
    c_even = rebuild_calcium(r_even, r0_even, start_even, beta0_even, gamma2);

    % predict y_even from the (averaged, neighboring) c_odd, and vice verse
    predict_even = (c_odd(1:end-1,:) + c_odd(2:end, :)) / 2;
    predict_even_part_nodc = predict_even(end-sample_pairs:end,:) - repmat(mean(predict_even(end-sample_pairs:end,:),1),sample_pairs+1,1);
    predict_odd = (c_even(1:end-1, :) + c_even(2:end, :)) / 2;
    predict_odd_part_nodc = predict_odd(end-sample_pairs:end,:) - repmat(mean(predict_odd(end-sample_pairs:end,:),1),sample_pairs+1,1);
    
    switch method
        case {1,2}
            err_even = mean((predict_even - y_even(1:end-1,:)).^2, 1);
            err_odd = mean((predict_odd - y_odd(2:end,:)).^2, 1);
        case {3,4,5}
            err_even = mean((predict_even_part_nodc - y_even_part_nodc).^2, 1);
            err_odd = mean((predict_odd_part_nodc - y_odd_part_nodc).^2, 1);
    end

    err_per_trace = (err_even + err_odd) / 2;
    err_mean_curve(i)     = mean(err_per_trace);
    err_std_curve(i) = std(err_per_trace);
end

[min_err, min_idx] = min(err_mean_curve);
threshold = min_err + err_std_curve(min_idx);

good_idx = find(err_mean_curve <= threshold);
best_idx = good_idx(end);   % largest (most noise removed) parameter still within threshold

best_param = param_range(best_idx);

figure;
if method == 5 || method == 3
    plot(param_range, err_mean_curve, 'o-');
    if method == 5 
        fprintf('The find param algorithm does not work for Lucy-Richardson, use the known param 10');
    end
else
    semilogx(param_range, err_mean_curve, 'o-');
end
hold on;
yline(threshold, 'k--');
plot(param_range(best_idx), err_mean_curve(best_idx), 'ro', 'MarkerFaceColor', 'r');
xlabel('parameter range');
ylabel('odd/even cross-validation error');
legend('error curve', 'min + std threshold', 'chosen parameter');
title('Parameter search');

end
