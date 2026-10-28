%% calibrate_dof.m -- empirical DOF of the iPDC null distribution
% Under H0, |iota_hat|^2 ~ Beta(1, n-1) with mean 1/n, where n is the
% effective DOF of the estimated (partial) coherence. This script measures
% n_emp = 1 / mean(|iota_hat|^2) on all null edges of a chain DGP
% 1 -> 2 -> ... -> N, over a grid of N, K (tapers) and M (independent
% segments), and fits
%       n = a*K*M + b*M + c*N + d
% to test the form used in the tutorial, n = L - N + 1 with
% L = (2K+1) M - 2, i.e. a = 2, b = 1, c = -1, d = -1.
clear; close all
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src')); cd(root); rng(5)

fs = 200; segLen = 512; nfft = segLen; overlap = 0.5; Mrep = 60;
Ns = [3 6 10]; NWs = [2 3 4 5]; Ms = [4 8 16];
f = (0:nfft-1).' * fs / nfft;
fk = f > 3 & f < fs/2 - 3;            % skip DC/Nyquist neighbourhoods

rows = [];
fprintf('%3s %3s %3s %4s | %8s %8s %8s | %8s\n', 'N','K','M','KM','n_emp','n_paper','ratio','nNull');
for N = Ns
    A = zeros(N, N, 2);
    for i = 1:N
        A(i,i,1) = 0.5 + 0.1*mod(i,3);  A(i,i,2) = -0.2 - 0.05*mod(i,2);
        if i > 1, A(i,i-1,1) = 0.4; end
    end
    Sig = eye(N);
    It = ipdc_from_var(A, Sig, fs, nfft);
    nullm = false(N); for i=1:N, for j=1:N, nullm(i,j) = i~=j && max(It(:,i,j)) < 1e-8; end, end
    for NW = NWs
        K = 2*NW - 1;
        for M = Ms
            T = segLen * M;
            acc = 0; cnt = 0;
            for r = 1:Mrep
                X = sim_var(T, A, Sig);
                [S, ~, mt] = multitaper_csd(X, fs, NW, K, segLen, overlap, nfft);
                [H, Sg] = wilson_factor_normalized(S);
                I2 = ipdc(H, Sg).^2;
                v = I2(fk, :, :); v = v(:, nullm);
                acc = acc + sum(v(:)); cnt = cnt + numel(v);
            end
            n_emp = cnt / acc;
            n_pap = mt.effDOF - N + 1;
            rows(end+1, :) = [N K M n_emp n_pap]; %#ok<SAGROW>
            fprintf('%3d %3d %3d %4d | %8.1f %8d %8.3f | %8d\n', N, K, M, K*M, n_emp, n_pap, n_emp/n_pap, nnz(nullm));
        end
    end
end

N = rows(:,1); K = rows(:,2); M = rows(:,3); n = rows(:,4);
X = [K.*M, M, N, ones(size(N))];
beta = X \ n;
fprintf('\nfree fit:  n = %.3f*K*M + %.3f*M + %.3f*N + %.3f\n', beta);
res = n - X*beta; fprintf('  RMSE %.2f, max |res| %.2f\n', sqrt(mean(res.^2)), max(abs(res)));
res_p = n - rows(:,5);
fprintf('paper form n = (2K+1)M - N - 1:  RMSE %.2f, max |res| %.2f, mean ratio %.3f\n', ...
        sqrt(mean(res_p.^2)), max(abs(res_p)), mean(n ./ rows(:,5)));
% constrained: n = a*K*M + c*N + d (no separate M term)
X2 = [K.*M, N, ones(size(N))]; b2 = X2 \ n; r2 = n - X2*b2;
fprintf('fit without M term:  n = %.3f*K*M + %.3f*N + %.3f   RMSE %.2f\n', b2, sqrt(mean(r2.^2)));
