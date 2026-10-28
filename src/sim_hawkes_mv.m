function events = sim_hawkes_mv(mu, alpha, beta, T_max)
%SIM_HAWKES_MV  Multivariate Hawkes process with exponential kernels.
%
%   events = sim_hawkes_mv(mu, alpha, beta, T_max)
%
%   Conditional intensity of component i:
%       lambda_i(t) = mu_i + sum_j sum_{s in events_j, s < t}
%                            alpha(i, j) * beta * exp(-beta * (t - s))
%
%   Simulation is by Ogata's thinning algorithm using the (decreasing)
%   intensity as its own upper bound between events (Ogata 1981).
%
%   Inputs
%     mu     : N x 1 baseline rates (events/sec), >= 0
%     alpha  : N x N excitation matrix; alpha(i, j) is the gain of an
%              event in j on the intensity of i. Stability requires the
%              spectral radius of alpha to be < 1. Entries must be >= 0
%              for the thinning bound to be valid.
%     beta   : scalar decay rate of the exponential kernel
%     T_max  : simulation horizon (sec)
%
%   Output
%     events : N x 1 cell; events{i} is a column vector of event times.

    assert(all(alpha(:) >= 0), 'alpha must be non-negative for this thinning scheme');
    assert(beta > 0, 'beta must be positive');

    N = numel(mu);
    events = cell(N, 1);
    for i = 1:N, events{i} = zeros(0, 1); end

    % Pre-allocate growable buffers
    cap = max(round(2 * sum(mu) * T_max), 1024);
    buf  = cell(N, 1);
    nev  = zeros(N, 1);
    for i = 1:N, buf{i} = zeros(cap, 1); end

    S = zeros(N, 1);   % current excitation state
    t = 0;
    while t < T_max
        lambda_bar = sum(mu + S);              % upper bound on next-window intensity
        if lambda_bar <= 0, break; end

        dt    = -log(rand) / lambda_bar;
        t_new = t + dt;
        if t_new >= T_max, break; end

        S          = S * exp(-beta * dt);      % decay
        lambda_new = mu + S;
        lam_tot    = sum(lambda_new);

        if rand <= lam_tot / lambda_bar
            % accept event; choose component
            c = cumsum(lambda_new) / lam_tot;
            i = find(rand < c, 1, 'first');
            nev(i) = nev(i) + 1;
            if nev(i) > numel(buf{i})
                buf{i} = [buf{i}; zeros(numel(buf{i}), 1)];   % double capacity
            end
            buf{i}(nev(i)) = t_new;
            S = S + alpha(:, i) * beta;        % self-update of state
        end
        t = t_new;
    end

    for i = 1:N
        events{i} = buf{i}(1:nev(i));
    end
end
