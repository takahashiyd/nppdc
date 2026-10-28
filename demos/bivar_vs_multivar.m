%% bivar_vs_multivar.m -- demo for paper §3
% Three DGPs designed to expose what bivariate analysis (ordinary
% coherence + bivariate iPDC on a pair) misses or fabricates relative
% to multivariate iPDC computed on the full system.
%
%   DGP A : chain  1 -> 2 -> 3                (spurious 1-3 from cascade)
%   DGP B : common driver  3 -> {1, 2}        (spurious 1-2 from confounding)
%   DGP C : cancellation: direct 1 -> 2 plus opposite indirect 1 -> 3 -> 2
%           (real direct edge that bivariate coherence misses)

clear; close all; root = fileparts(fileparts(mfilename('fullpath')));   % project root (code/)
addpath(fullfile(root, 'src'), fullfile(root, 'demos')); cd(root);   % data/ and figures/ are relative to root
rng(0)

fs     = 200;
T      = 8192;
nfft   = 512;
NW     = 4;
K      = 7;
segLen = 512;
overl  = 0.5;

%% ----- DGPs ------------------------------------------------------------
% (A) Chain 1 -> 2 -> 3
A_A = zeros(3,3,2);
A_A(:,:,1) = [ 0.55  0    0;
               0.50  0.30 0;
               0    0.50 0.30];
A_A(:,:,2) = [-0.25 0    0;
               0   -0.10 0;
               0    0   -0.20];
S_A = diag([1 1 1]);

% (B) Common driver 3 -> 1, 3 -> 2, no 1<->2.
%     Strong driving so the spurious 1-2 coupling is well above noise.
A_B = zeros(3,3,1);
A_B(:,:,1) = [ 0.20 0    0.85;
               0    0.20 0.85;
               0    0    0.40];
S_B = diag([1 1 1]);

% (C) Cancellation: direct 1->2 (lag 1) plus indirect 1->3->2 that exactly
%     cancels the direct effect at the bivariate level for DC frequencies.
%     Construction: x_3(n) = -x_1(n-1) + w_3(n), x_2(n) = a*x_1(n-1) +
%     a*x_3(n-1) + w_2(n). Substituting yields x_2(n) = a*x_1(n-1) -
%     a*x_1(n-2) + correlated noise, so the bivariate transfer
%     H_{21}^biv(f) = a*e^{-iw}(1 - e^{-iw}) vanishes at f = 0 and peaks
%     at Nyquist. Multivariate iPDC shows both direct edges 1->2 and 1->3
%     constant across f.
A_C = zeros(3,3,1);
a   = 0.5;
A_C(:,:,1) = [ 0   0   0;
               a   0   a;
              -1   0   0];
S_C = diag([1 1 1]);

cases = struct( ...
  'name',  {'Chain 1\rightarrow 2\rightarrow 3', ...
            'Common driver 3\rightarrow\{1,2\}', ...
            'Cancellation: 1\rightarrow 2 direct + 1\rightarrow 3\rightarrow 2'}, ...
  'A',     {A_A, A_B, A_C}, ...
  'Sigma', {S_A, S_B, S_C}, ...
  'pair',  {[1 3], [1 2], [1 2]} ...   % the pair where the issue appears
);

%% ----- analyze each DGP ------------------------------------------------
res = repmat(struct(), 1, numel(cases));
for c = 1:numel(cases)
    cs   = cases(c);
    X    = sim_var(T, cs.A, cs.Sigma);

    % --- multivariate (3-channel) iPDC ---
    [S3, fW, mt3]  = multitaper_csd(X, fs, NW, K, segLen, overl, nfft);
    [H3, Sg3]      = wilson_factor_normalized(S3);
    I3             = ipdc(H3, Sg3);

    % --- bivariate iPDC on the pair only ---
    Xp             = X(cs.pair, :);
    [S2, fW2, mt2] = multitaper_csd(Xp, fs, NW, K, segLen, overl, nfft);
    [H2, Sg2]      = wilson_factor_normalized(S2);
    I2             = ipdc(H2, Sg2);

    % --- bivariate ordinary coherence on the pair ---
    s11 = real(S2(:,1,1));
    s22 = real(S2(:,2,2));
    s12 =      S2(:,1,2);
    coh = abs(s12).^2 ./ max(s11 .* s22, eps);

    % --- parametric ground truth on full system ---
    [Itrue, ftrue] = ipdc_from_var(cs.A, cs.Sigma, fs, nfft);

    res(c).f         = fW;
    res(c).coh       = coh;
    res(c).I_biv     = I2;        % nfft x 2 x 2 (over the pair)
    res(c).I_mult    = I3;        % nfft x 3 x 3 (full)
    res(c).I_true    = Itrue;
    res(c).pair      = cs.pair;
    res(c).name      = cs.name;

    % Console summary
    p = cs.pair;
    fprintf('\nDGP %s  (pair = [%d %d])\n', cs.name, p(1), p(2));
    fprintf('  max bivariate ordinary coherence on pair  : %.3f\n', max(coh));
    fprintf('  max bivariate iPDC %d->%d                  : %.3f\n', p(1), p(2), max(squeeze(I2(:,2,1))));
    fprintf('  max bivariate iPDC %d->%d                  : %.3f\n', p(2), p(1), max(squeeze(I2(:,1,2))));
    fprintf('  max multivariate iPDC %d->%d (estimate)    : %.3f\n', p(1), p(2), max(squeeze(I3(:,p(2),p(1)))));
    fprintf('  max multivariate iPDC %d->%d (estimate)    : %.3f\n', p(2), p(1), max(squeeze(I3(:,p(1),p(2)))));
    fprintf('  max multivariate iPDC %d->%d (truth)       : %.3f\n', p(1), p(2), max(squeeze(Itrue(:,p(2),p(1)))));
    fprintf('  max multivariate iPDC %d->%d (truth)       : %.3f\n', p(2), p(1), max(squeeze(Itrue(:,p(1),p(2)))));
end

%% ----- figure ----------------------------------------------------------
fig = figure('Position', [50 50 1200 850], 'Color', 'w');
nyq = fs/2;

for c = 1:numel(res)
    r  = res(c);
    p  = r.pair;
    f  = r.f;
    msk = f <= nyq;

    % (col 1) bivariate ordinary coherence
    subplot(3, 3, (c-1)*3 + 1)
    plot(f(msk), r.coh(msk), 'k', 'LineWidth', 1.6); hold on
    ylim([0 1]); xlim([0 nyq]); grid on
    title(sprintf('%s\nBivar. ord. coherence |C_{%d%d}(f)|^2', r.name, p(1), p(2)))
    xlabel('Hz')

    % (col 2) bivariate iPDC (on the pair only)
    subplot(3, 3, (c-1)*3 + 2)
    plot(f(msk), squeeze(r.I_biv(msk, 2, 1)), 'r', 'LineWidth', 1.6); hold on
    plot(f(msk), squeeze(r.I_biv(msk, 1, 2)), 'b', 'LineWidth', 1.6);
    ylim([0 1]); xlim([0 nyq]); grid on
    legend({sprintf('%d\\rightarrow%d', p(1), p(2)), ...
            sprintf('%d\\rightarrow%d', p(2), p(1))}, 'Location', 'Best')
    title('Bivariate iPDC (pair only)')
    xlabel('Hz')

    % (col 3) multivariate iPDC + true
    subplot(3, 3, (c-1)*3 + 3)
    plot(f(msk), squeeze(r.I_mult(msk, p(2), p(1))), 'r', 'LineWidth', 1.6); hold on
    plot(f(msk), squeeze(r.I_mult(msk, p(1), p(2))), 'b', 'LineWidth', 1.6);
    % truth in thick black on top; ylim starts below 0 so a zero line stays visible
    plot(f(msk), squeeze(r.I_true(msk, p(2), p(1))), 'k--', 'LineWidth', 2.2);
    plot(f(msk), squeeze(r.I_true(msk, p(1), p(2))), 'k:',  'LineWidth', 2.2);
    ylim([-0.04 1]); xlim([0 nyq]); grid on
    legend({sprintf('%d\\rightarrow%d (est.)', p(1), p(2)), ...
            sprintf('%d\\rightarrow%d (est.)', p(2), p(1)), ...
            sprintf('%d\\rightarrow%d (truth)', p(1), p(2)), ...
            sprintf('%d\\rightarrow%d (truth)', p(2), p(1))}, ...
            'Location', 'Best', 'FontSize', 7)
    title('Multivariate iPDC (full system)')
    xlabel('Hz')
end

style_paper_figure(fig);
set(fig, 'PaperPositionMode', 'auto')
print(fig, 'figures/bivar_vs_multivar.png', '-dpng', '-r150')
fprintf('\nSAVED figures/bivar_vs_multivar.png\n');
