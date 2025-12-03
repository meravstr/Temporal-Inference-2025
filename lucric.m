function [r] = lucric(y,gamma,p_num,iter)
% This function infers the rate r from the fluorescence recordings y, 
% assuming the calcium decays to gamma of its value each time bin with no spikes accrued. 
% The function uses Lucy-Richardson iterative algorithm to calculate r 
% NOTATIONS: INPUTS y is t x n - time by traces matrix ; 
% OUTPUTS r is t x n matrix

t = size(y,1);
n = size(y,2);

r = zeros(size(y));

p = 0:1:p_num;
conv_kernel = gamma.^p;
kernel_for_lucy=[zeros(size(conv_kernel)) conv_kernel]';

if iter == []
    iter = 10;
end

for i = 1:n
    cur_y = y(:,i)-min(y(:,i));
    r(:,i) = deconvlucy(cur_y,kernel_for_lucy,iter);
end

