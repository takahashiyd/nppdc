%% hodge_jack_validate.m -- Monte Carlo calibration of the Hodge jackknife inference
% Per-node and global tests of IPDC_HODGE_JACK on three VAR(2) DGPs
% (M = 100 replicates, 16 jackknife blocks):
%
%   symm  : exchangeable channels with equal coupling both ways and
%           Sigma = I, so every off-diagonal |iPDC| is identical and s = 0
%   indep : independent channels
%   chain : 1 -> 2 -> 3 with couplings 0.90 and 0.75 (as in hodge_demo.m;
%           with equal couplings the middle node would have a zero score)

clear; close all
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src')); cd(root); rng(11)

fs = 200; T = 4096; NW = 4; K = 7; seg = 256; nfft = seg;   % T/seg = 16 blocks
alpha = 0.05; M = 100; N = 3;
fband = @(f) f >= 2 & f <= fs/2 - 2;

dgp.symm.A1 = [0.40 0.20 0.20; 0.20 0.40 0.20; 0.20 0.20 0.40];
dgp.symm.A2 = -0.15 * eye(3);
dgp.indep.A1 = diag([0.55 0.30 0.40]);   dgp.indep.A2 = diag([-0.25 -0.10 -0.20]);
dgp.chain.A1 = [0.55 0 0; 0.90 0.30 0; 0 0.75 0.40];   % true order 1 < 2 < 3 at every in-band f
dgp.chain.A2 = [-0.25 0 0; 0 -0.10 0; 0 0 -0.20];

names = fieldnames(dgp);
for d = 1:numel(names)
    nm = names{d}; A = cat(3, dgp.(nm).A1, dgp.(nm).A2);
    if strcmp(nm, 'symm'), Sig = eye(N); else, Sig = diag([1 0.8 1.2]); end
    comp = [reshape(A, N, 2*N); eye(N) zeros(N)];
    assert(max(abs(eig(comp))) < 1, '%s DGP unstable', nm);
    It = ipdc_from_var(A, Sig, fs, nfft); P = eye(N) - ones(N)/N;
    strue = zeros(nfft, N);
    for k = 1:nfft, R = squeeze(It(k,:,:)).^2; strue(k,:) = ((P*sum(R-R.',2))/N).'; end

    rej_node = 0; cov_node = 0; rej_glob = 0; order_ok = 0; nf = 0;
    for r = 1:M
        X = sim_var(T, A, Sig);
        o = ipdc_hodge_jack(X, fs, NW, K, seg, 0.5, nfft, alpha);
        fk = fband(o.f);
        rej_node = rej_node + mean(o.score_reject(fk, :), 'all');
        inside   = strue(fk,:) >= o.score_ci(fk,:,1) & strue(fk,:) <= o.score_ci(fk,:,2);
        cov_node = cov_node + mean(inside, 'all');
        rej_glob = rej_glob + mean(o.reject_glob(fk));
        sb = mean(o.scores(fk, :), 1);
        order_ok = order_ok + double(sb(1) < sb(2) && sb(2) < sb(3));
    end
    fprintf('%-6s  per-node reject = %.3f   per-node CI coverage = %.3f   global reject = %.3f   1<2<3 order = %.2f\n', ...
            nm, rej_node/M, cov_node/M, rej_glob/M, order_ok/M);
end
