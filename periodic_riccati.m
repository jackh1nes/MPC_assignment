function [P, K] = periodic_riccati(A, B, Q, R)
% Periodic Riccati solver
%   Finds an approximation to the solution of the algebraic riccati
%   equation in the case where state space matrices are time varying. Here,
%   only B varies in time
arguments (Input)
    A % nxn matrix - doesn't vary in time
    B % nxmxT matrix set - n state dimension, m inputs & periodic about T samples
    Q % nxn symmetrical matrix - state penalty
    R % mxm symmetrical matrix - input penalty
end

arguments (Output)
    P % Solution to the periodic riccati equation nxnxT - use for terminal cost function
    K % Optimal feedback gain, mxnxT - used as x^{dot} = (A+BK)x 
end

% Taken for now from here
% https://elib.dlr.de/12280/1/varga_ifac2005p2.pdf


n = size(A,1); m = size(B,2); T = size(B,3);

% Instantiate P
P = zeros(n,n,T);
for i = 1:T
    P(:,:,i) = eye(n);
end

k = T; 

eps = zeros(T,1);

for iteration = 1:10000
    kp1 = k+1; 
    if kp1 > T
        kp1 = 1;
    end
    P_new = Q + A'*P(:,:,kp1)*A - (A'*P(:,:,kp1)*B(:,:,k))*(R+B(:,:,k)'*P(:,:,kp1)*B(:,:,k))^(-1) * (A'*P(:,:,kp1)*B(:,:,k))';
    
    eps(k) = norm(P(:,:,k) - P_new,'fro');
    P(:,:,k) = P_new;
    k = k - 1;
    if k == 0
        k = T;
    end
    if max(eps) < 1e-3
        break
    end
end

K = zeros(m,n,T);
for k = 1:T
    kp1 = k+1; 
    if kp1 > T 
        kp1 = 1; 
    end
    K(:,:,k) = -(B(:,:,k)' * P(:,:,kp1) * B(:,:,k) + R)^(-1) * B(:,:,k)' * P(:,:,kp1) * A;
end
end