function [r] = fft_weiner(y,gamma,k)
% This function infers the rate r from the fluorescence recordings y, 
% assuming the calcium decays to gamma of its value each time bin with no spikes accrued. 
% The function uses Lucy-Richardson iterative algorithm to calculate r 
% NOTATIONS: INPUTS y is t x n - time by traces matrix ; 
% OUTPUTS r is t x n matrix

% pedding y
t_org = size(y,1);
y = [flipud(y);y;flipud(y)];
t = size(y,1);

t_filter = 1:min(200,t_org);
n = size(y,2);

tau = -1/log(gamma);   % decay time‐constant - in units 1t 
h = exp(-t_filter/tau)';
h = fft(h,t);

f = fft(y);
w = conj(h) ./ (abs(h).^2 + k);
w = repmat(w,1,n);
r_w = w.*f;
r = real(ifft(r_w));
r = r(t_org+1:t_org*2,:);


