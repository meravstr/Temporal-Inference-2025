function [r,r0,beta0] = convar(y,gamma,lambda)
% This function implements the continuously varying (convar) method. 
% It infers the spiking rate r from the fluorescence recordings y, 
% assuming calcium decays by gamma each time bin when no spikes occur, 
% and that the penalty for a (squared) change in the spiking rate is lambda. 
% Inputs: 
% y - fluorescence. t x n, time by trials/pixels matrix ; 
% gamma - the decay of the calcium during a time bin. A number between 0 and 1 (typically close to 1); 
% lambda - the penalty weight, a real number. 
% Outputs:
% r - the inferred spiking rate from t=2 to t=time. t-1 x n matrix. 
% r0 - r at t=1, carries no direct biological meaning . A real number.
% beta0 - the shift between the fluorescence and the spiking rate.

t = size(y,1);
n = size(y,2);

Dinv = zeros(t);
insert_vec = 1;
for i_t = 1:t
    Dinv(i_t,1:i_t) = insert_vec;
    insert_vec = [gamma^i_t, insert_vec];
end
P = eye(t)-1/t*ones(t);

ytilde = P*y;
A = P*Dinv;
L = [zeros(t,1) [zeros(1,t-1); [zeros(1,t-1); [-eye(t-2), zeros(t-2,1)] + [zeros(t-2,1), eye(t-2)]]]];
Z = L'*L;
multiplies_r = (A'*A+lambda*Z);
multiplies_r = pinv(multiplies_r);
r_anlytic = multiplies_r*A'*ytilde;

d = -min([r_anlytic(2:end,:); zeros(1,n)],[],1);
rm_ratio = A'*A*ones(t,1)./(A'*A*[1; zeros(t-1,1)]);
rm = -d*rm_ratio(1);
r = r_anlytic+repmat(d,t,1)+[rm;zeros(t-1,n)];

beta0 = 1/t*ones(1,t)*(y-Dinv*r);

r0 = r(1,:);
r = r(2:end,:);

end

