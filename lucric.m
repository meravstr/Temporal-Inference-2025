function [r] = lucric(y,gamma,p_num,iter)
% This function implements the Lucy-Richardson algorithm for temporal
% inference of aggregated fluorescence recordings
% It infers the rate r from the fluorescence recordings y, 
% assuming the noise follows a Poisson distribution 
% and that calcium decays to gamma of its value if no spikes occurred within a time bin.

% Inputs: 
% y - fluorescence. t x n, time by trials/pixels matrix ; 
% gamma - the decay of the calcium during a time bin. A number between 0 and 1 (typically close to 1); 
% p_num - the number of time points the calcium decay kernel spans (a
% natural number). It has to be smaller than t/2.
% iter - the number of iterations the algorithm perfroms (a natural number;
% iter = 10 if no iter is given)
% Outputs:
% r - the inferred spiking rate from t=1 to t=time. t x n matrix. 

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


