function [r] = fft_wiener(y,gamma,k)

% This function implements the Wiener filter deconvolution method. 
% It infers the spiking rate r from the fluorescence recordings y, 
% assuming calcium decays by gamma each time bin when no spikes occur, 
% and that the ratio between noise and signal is fixed. 
% It works in Fourier space.

% Inputs: 
% y - fluorescence. t x n, time by trials/pixels matrix ; 
% gamma - the decay of the calcium during a time bin. A number between 0 and 1 (typically close to 1); 
% k - relates to inverse signal to noise ratio. 
% Outputs:
% r - the inferred spiking rate from t=2 to t=time. t-1 x n matrix. 

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


