%% calcium_decay_finder.m
% Finds gamma - the fraction of calcium signal left after one time bin -
% for your indicator at your recording rate.
% gamma = exp(-dt/tau), with dt the time bin and tau the time its takes to decay into 1/e.
% Hence, a gamma known at one rate converts to any other rate as
% gamma = gamma_ref ^ (rate_ref / rate_new)

recording_rate = 30;           % <<< UPDATE HERE: your recording rate (Hz)
indicator = 'GCaMP6s';         % <<< UPDATE HERE: 'GCaMP6f', 'GCaMP6s', or 'custom'

switch indicator
    case 'GCaMP6f'
        gamma_ref = 0.97;  rate_ref = 40;
    case 'GCaMP6s'
        gamma_ref = 0.95;  rate_ref = 10;   
    case 'custom'
        gamma_ref = [];  rate_ref = [];   % <<< UPDATE HERE: a gamma you know, and the rate (Hz) it was measured at
end

tau   = -1 / (rate_ref * log(gamma_ref));          % decay time constant (s)
gamma = gamma_ref ^ (rate_ref / recording_rate);

fprintf('%s: tau = %.3f s, gamma = %.4f at %g Hz.\n', indicator, tau, gamma, recording_rate);
