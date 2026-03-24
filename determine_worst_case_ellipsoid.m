function alpha = determine_worst_case_ellipsoid(a,b,P)
% determine_worst_case_ellipsoid find the upper bound alpha such that x^TPx
% <= alpha is the largest control invariant set 
%   This is for the time varying case
arguments (Input)
    a % NxnxT - coefficients of x in inequalities ax<=b
    b % Nx1x - the other part of those inequalities
    P % nxnxT - solution of the discrete time riccati equation, periodic with period T
end

arguments (Output)
    alpha % upper norm bound of largest control invariant set
end

N = size(a,1); T = size(P,3);
alpha = inf*ones(T,1);

for k = 1:T
    Pinv = inv(P(:,:,k));
    for i = 1:N
        aik = a(i,:,k)';
        bik = b(i,:,k);
        
        alpha_ik = (bik^2) / (aik' * Pinv * aik);
        alpha(k) = min(alpha(k), alpha_ik);
    end
end