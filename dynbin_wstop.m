function [r_final,r0,beta_0,iter] = dynbin_wstop(y,gamma,lambda,err_or_iter)

% This function implements the dynamically binned (dynbin) algorithm. 
% It infers the spiking rate r from the fluorescence recordings y, 
% using an iterative algorithm.
% It assumes the calcium decays by gamma each time bin when no spikes occur, 
% and that the penalty for an absolute change in the spiking rate is lambda.
% It retrieves constant spiking rates within dynamically determent time
% spans.
% Inputs: 
% y - fluorescence. t x n, time by trials/pixels matrix ; 
% gamma - the decay of the calcium during a time bin. A number between 0 and 1 (typically close to 1); 
% lambda - the penalty weight, a real number. 
% err_or_iter – a stopping criterion for the iterations. 
% If it is a number greater than 1, it specifies the number of iterations to run. 
% If it is a fraction (between 0 and 1), iterations are terminated once the average change in the inferred spiking rate falls below the given fraction from the mean absolute range of r. 
% If left unspecified, the default is 1000 iterations. 
% Accepts either a fraction (typically small like 0.01) or a natural number.
% Outputs:
% r - the inferred spiking rate from t=2 to t=time. t-1 x n matrix. 
% r0 - r at t=1, carries no direct biological meaning . A real number.
% beta0 - the shift between the fluorescence and the spiking rate.
% iter - the number of iterations run by the algorithm

t = size(y,1);
n = size(y,2);

Dinv = zeros(t);
insert_vec = 1;
for i_t = 1:t
    Dinv(i_t,1:i_t) = insert_vec;
    insert_vec = [gamma^i_t, insert_vec];
end

P = eye(t)-1/t*ones(t);
tildey = P*y;
A = P*Dinv;
% largest step size that ensures convergence
s = 0.5*((1-gamma)/(1-gamma^t))^2;

% initializing
r = rand(size(y));
if isempty(err_or_iter)
    err_or_iter = 1000;
end

if err_or_iter>1
    for i = 1:ceil(err_or_iter)
        Ar = A*r;
        tmAr = (tildey-Ar);
        At_tmAr = A'*tmAr;
        x = r + s*At_tmAr;
        for j = 1:size(y,2)
            r(2:end,j) = fTVdenoise(s*lambda,x(2:end,j));
        end
        r(r<0) = 0;
        r(1,:) = x(1,:);
    end
    
else
    i = 1;
    relative_change = zeros(size(y,2),1);
    test_err = 1; % percentage of error compared to rate magnitude
    indx = 1:size(y,2);
    to_update = ones(size(y,2),1);
    while test_err > (err_or_iter/2)
        r_old = r;
        Ar = A*r;
        tmAr = (tildey-Ar);
        At_tmAr = A'*tmAr;
        x = r + s*At_tmAr;
        for j = 1:size(y,2)
            if indx(to_update)
            r(2:end,j) = fTVdenoise(s*lambda,x(2:end,j));
            else 
            r(2:end,j) = r_old(2:end,j);
            end
        end
        r(r<0) = 0;
        for j = indx(to_update)
           relative_change(j) = mean(abs(r_old(2:end,j)-r(2:end,j)))/(1+abs(max((r(2:end,j)-min(r(2:end,j)))))); %mean(abs(r_old(2:end,j)-r(2:end,j)))/mean((r(2:end,j)-min(r(2:end,j))));
           if relative_change(j) < (err_or_iter/2)
               to_update(j) = 0;
           end
        end
        r(1,:) = x(1,:);
        i = i+1;
        test_err = max(relative_change);
    end
end
    
r_final = r(2:end,:);
r0 = r(1,:);
beta_0 = mean(y-Dinv*r);
iter = i;

end

