function [r] = firdif(y,gamma,smt)
% This function implements the first diffrences method
% It infers the spiking rate from the fluorescence by calculating r_t = c_t - \gamma c_(t-1).
% The result is then smoothed by smt nearest points.

% Inputs: 
% y - fluorescence. t x n, time by trials/pixels matrix ; 
% gamma - the decay of the calcium during a time bin. A number between 0 and 1 (typically close to 1); 
% lambda - the penalty weight, a real number. 
% Outputs:
% r - the inferred spiking rate from t=2 to t=time. t-1 x n matrix. 

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

