function [r, r_idx_start, r0, beta0] = run_method(y, gamma, method, param)
% RUN_METHOD infers the spiking rate from fluorescence y using one of
% the five methods in this repository, chosen by METHOD:
%   1 = Continuously-Varying (Convar)      param = lambda
%   2 = Dynamically-Binning                param = lambda
%   3 = First-Differences                  param = smt (smoothing window)
%   4 = Wiener-Filter                      param = k (related to inverse SNR)
%   5 = Lucy-Richardson                    param = iter (# iterations)
%
% This is the same code as Part 3 of run_inference.m, wrapped as a
% function (it is also used by search_best_param_oddeven.m).
%
% INPUTS y - fluorescence traces T x P (time x pixels/trials
% gamma - the calcium decay per time bin; 
% method - a number 1-5 (above); 
% param - a parameter for the method.
% OUTPUTS r - the inferred rate; 
% r_idx_start - the row of y that the first row of r corresponds to: 
% Convar, Dynamically-Binning and First-Differences do not produce a rate for t=1 (r_idx_start = 2),
% while Wiener-Filter and Lucy-Richardson do (r_idx_start = 1).
% r0, beta0 (Convar and Dynamically-Binning only, empty otherwise) - the
% initial calcium and baseline estimates, needed to rebuild calcium.

n_t = size(y,1); % for method 5
r0 = [];
beta0 = [];

switch method
    case 1  % Continuously-Varying
        [r, r0, beta0] = convar(y, gamma, param);
        r_idx_start = 2;

    case 2  % Dynamically-Binning
        [r, r0, beta0, ~] = dynbin_wstop(y, gamma, param, 0.0001);
        r_idx_start = 2;

    case 3  % First-Differences
        r = firdif(y, gamma, param);
        r_idx_start = 2;
        r0 = y(1,:);

    case 4  % Wiener-Filter
        r = fft_wiener(y, gamma, param);
        r_idx_start = 1;

    case 5  % Lucy-Richardson
        % p_num sets the length of the calcium-decay filter kernel. 
        % It is restricted by the data size, up to half of its length. 
        % It is set here long enough for gamma^p_num to decay below 0.1%, 
        % or third of the trial timepoints, whichever is shorter.
        p_num = min([ceil(log(0.001) / log(gamma)) floor(n_t/3)]);
        r = lucric(y, gamma, p_num, param);
        r_idx_start = 1;

    otherwise
        error('run_method:unknownMethod', 'method must be an integer 1-5.');
end

end
