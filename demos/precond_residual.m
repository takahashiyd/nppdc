%% precond_residual.m -- preconditioning and the Wilson factorization residual
%
% Reproduces the worked example of Sec. 4.3: coherency preconditioning
% reduces the factorization residual
%
%       rho = || S - psi psi' ||_F / || S ||_F
%
% by up to seven orders of magnitude, increasingly with the channel-power
% spread and with N, while the scale-invariant iPDC is unchanged.
%
% Two sweeps on stable random VAR(2) systems:
%   A. fixed N=16, vary channel-power spread (cond(S))
%   B. fixed power spread, vary N
%
% Compares WILSON_FACTOR (baseline) vs WILSON_FACTOR_NORMALIZED (preconditioned).

clear; close all; root = fileparts(fileparts(mfilename('fullpath')));   % project root (code/)
addpath(fullfile(root, 'src'), fullfile(root, 'demos')); cd(root);   % data/ and figures/ are relative to root
rng(0)

fs = 200; T = 8192; segLen = 512; overl = 0.5; NW = 4; K = 7; nfft = segLen;

%% ---- helper: a stable, somewhat sparse random VAR(2) ---------------------
make_var = @(N) local_stable_var(N);

%% ---- A. fixed N, vary channel-power spread -------------------------------
N = 16;
[A, ~] = make_var(N);
X0 = sim_var(T, A, eye(N));               % homogeneous innovations
Itrue = ipdc_from_var(A, eye(N), fs, nfft);   % iPDC is gain-invariant

spreads = [1e2 1e4 1e6 1e8 1e10 1e12 1e14 1e16];   % channel-POWER ratio across channels

fprintf('\n=== A. N=%d, vary channel-power spread ===\n', N);
fprintf('%-10s %-9s | %-11s %-11s %-9s | %-9s %-9s\n', ...
    'powerspr','cond(S)','rho_base','rho_prec','base/prec','RMSE_b','RMSE_p');
for sp = spreads
    g = logspace(-log10(sp)/4, log10(sp)/4, N).';   % gains; power ratio = sp
    [rb, rn, cnd, eb, en] = compare(X0 .* g, fs, NW, K, segLen, overl, nfft, Itrue);
    fprintf('%-10.0e %-9.1e | %-11.2e %-11.2e %-9.1e | %-9.4f %-9.4f\n', ...
        sp, cnd, rb, rn, rb/rn, eb, en);
end

%% ---- B. fixed power spread, vary N ---------------------------------------
sp = 1e8;
Ns = [4 8 16 24 32];
fprintf('\n=== B. channel-power spread = %.0e, vary N ===\n', sp);
fprintf('%-4s %-9s | %-11s %-11s %-9s | %-9s %-9s\n', ...
    'N','cond(S)','rho_base','rho_prec','base/prec','RMSE_b','RMSE_p');
for N = Ns
    [A, ~] = make_var(N);
    X0 = sim_var(T, A, eye(N));
    Itrue = ipdc_from_var(A, eye(N), fs, nfft);
    g = logspace(-log10(sp)/4, log10(sp)/4, N).';
    [rb, rn, cnd, eb, en] = compare(X0 .* g, fs, NW, K, segLen, overl, nfft, Itrue);
    fprintf('%-4d %-9.1e | %-11.2e %-11.2e %-9.1e | %-9.4f %-9.4f\n', ...
        N, cnd, rb, rn, rb/rn, eb, en);
end

fprintf(['\nrho_base/rho_prec spans ~10^1 (mild) to ~10^7 (severe) across these\n', ...
         'settings; the iPDC RMSE is essentially identical (scale-invariance).\n']);

%% ---- local functions -----------------------------------------------------
function [rb, rn, cnd, eb, en] = compare(X, fs, NW, K, segLen, overl, nfft, Itrue)
    [S, ~, ~] = multitaper_csd(X, fs, NW, K, segLen, overl, nfft);
    nf = size(S, 1);
    Snrm = 0; cvals = zeros(nf, 1);
    for k = 1:nf
        Sk = squeeze(S(k,:,:));
        Snrm = Snrm + norm(Sk, 'fro');
        cvals(k) = cond(Sk);
    end
    Snrm = Snrm / nf;  cnd = median(cvals);
    [Hb, Sgb, ~, ib] = wilson_factor(S);
    [Hn, Sgn, ~, in] = wilson_factor_normalized(S);
    rb = ib.factor_err / Snrm;            % relative factorization residual
    rn = in.factor_err / Snrm;
    Ib = ipdc(Hb, Sgb);  Inr = ipdc(Hn, Sgn);
    eb = sqrt(mean((Ib(:)  - Itrue(:)).^2));
    en = sqrt(mean((Inr(:) - Itrue(:)).^2));
end

function [A, me] = local_stable_var(N)
    ok = false;
    while ~ok
        A = zeros(N, N, 2);
        A(:,:,1) = (randn(N)/N)     .* (rand(N) < 0.35);
        A(:,:,2) = (randn(N)/(2*N)) .* (rand(N) < 0.20);
        comp = [reshape(A, N, 2*N); eye(N) zeros(N)];
        me = max(abs(eig(comp)));
        if me < 1.5
            A = A * (0.92/me);
            me = max(abs(eig([reshape(A, N, 2*N); eye(N) zeros(N)])));
            if me < 0.97, ok = true; end
        end
    end
end
