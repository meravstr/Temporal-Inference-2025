function c = rebuild_calcium(r, r0, r_idx_start, shift, gamma)
% REBUILD_CALCIUM recovers an inferred calcium trace
% from an inferred spiking rate, via the model's own recursion
% c_t = gamma*c_(t-1) + r_t.
%
% Used by search_best_param_oddeven.m to turn an inferred rate back into
% something directly comparable to fluorescence.
%
% INPUTS r - the spiking rate 
% r0 - the spiking rate set artificially for r at t=1 if the r_idx_start is
% 2.
% shift - set to beta0 (calcium fluorescence base line builder) if known,
% otherwise 0.
% gamma - the calcium decay in a time bin.

n_t = size(r,1)+r_idx_start-1;
n_p = size(r,2);

Dinv = zeros(n_t);
insert_vec = 1;
for i_t = 1:length(Dinv)
    Dinv(i_t,1:i_t) = insert_vec;
    insert_vec = [gamma^i_t, insert_vec];
end

if isempty(shift)
    shift = zeros(1,n_p);
end

switch r_idx_start
    case 1  % Wiener-Filter / Lucy-Richardson
        c = Dinv*r;
    case 2  % Continuously-Varying / Dynamically-Binning / First-Differences
        c = Dinv*[r0;r]+repmat(shift,n_t,1);
    
    otherwise
        error('rebuild_calcium:unknownR_idx_start', 'r_idx_start must be an integer 1-2.');
end