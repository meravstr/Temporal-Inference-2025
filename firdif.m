function [r] = firdif(y,gamma,smt)
% This function estimates the rate r from the fluorescence recordings y.
% It assumes the calcium decays to gamma of its value each time bin if no
% spikes accrued. It calculates r_t = c_t - \gamma c_(t-1).
% The results are then smoothed by smt nearest ponts.
% NOTATIONS: INPUTS y is t x n - time by traces matrix ; 
% gamma a number between 0 and 1 (typically close to 1); smt a number
% OUTPUTS r is t-1 x n matrix ;

t = size(y,1);
n = size(y,2);

D = [zeros(1,t); [-gamma*eye(t-1) zeros(t-1,1)]] + eye(t);
r = D*y;

% smoothing the results (without r(1) which is c(1) and not a spiking rate)
r_long = [flipud(r(3:3+floor(smt/2)-1,:)); r(2:end,:); flipud(r(end-floor(smt/2):end-1,:))];
r_smoothed = zeros(size(r_long));
for i = 1:n
    r_smoothed(:,i) = smooth(r_long(:,i),smt,'moving');
end
r = r_smoothed(floor(smt/2)+1:end-floor(smt/2),:);
r = r - min([min(r);zeros(1,n)]);

end

