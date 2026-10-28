function X = sim_var(T, A, Sigma, burn)
%SIM_VAR  Simulate a stable VAR(p) process.
%
%   X = sim_var(T, A, Sigma)
%   X = sim_var(T, A, Sigma, burn)
%
%   Simulates x(t) = sum_{r=1..p} A_r * x(t-r) + e(t), with e(t) ~ N(0, Sigma).
%
%   Inputs
%     T     : number of samples to return (after burn-in)
%     A     : N x N x p VAR coefficient tensor; A(:,:,r) is the lag-r matrix
%     Sigma : N x N innovation covariance (symmetric positive-definite)
%     burn  : burn-in length (default 500). The simulator runs T+burn
%             steps and discards the first `burn`.
%
%   Output
%     X     : N x T multivariate time series

    if nargin < 4 || isempty(burn), burn = 500; end
    [N, ~, p] = size(A);
    Ttot = T + burn;

    L = chol(Sigma, 'lower');
    E = L * randn(N, Ttot);

    X = zeros(N, Ttot);
    for t = (p+1):Ttot
        z = E(:, t);
        for r = 1:p
            z = z + A(:,:,r) * X(:, t-r);
        end
        X(:, t) = z;
    end
    X = X(:, burn+1:end);
end
