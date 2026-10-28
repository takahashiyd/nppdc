%% hodge_demo.m -- inferring directed hierarchy from iPDC (Hodge decomposition)
% Left:  simulated chain 1->2->3, ranking scores +/- 95% CI vs frequency
%        against the closed-form truth (jackknife inference, ipdc_hodge_jack).
% Right: Allen V1<->LM LFP (16 channels), band-averaged ranking score per
%        channel +/- 95% CI (needs data/allen/subset_movie.mat).
% Reproduces Figure fig:hodge of the tutorial.

close all; root = fileparts(fileparts(mfilename('fullpath')));   % project root (code/)
addpath(fullfile(root, 'src'), fullfile(root, 'demos')); cd(root);   % data/ and figures/ are relative to root
rng(8)
fs = 200; T = 4096; NW = 4; K = 7;
za = norminv(0.975);

%% ---- A. simulated chain 1 -> 2 -> 3 -----------------------------------
% Unequal couplings so that node 2 has a nonzero net inflow (same chain as
% hodge_jack_validate.m).
A1 = [0.55 0 0; 0.90 0.30 0; 0 0.75 0.40];
A2 = [-0.25 0 0; 0 -0.10 0; 0 0 -0.20];
A = cat(3, A1, A2); Sig = diag([1 0.8 1.2]); N = 3;
segA = 256; nfftA = segA;                 % 16 independent jackknife blocks

X = sim_var(T, A, Sig);
oA = ipdc_hodge_jack(X, fs, NW, K, segA, 0.5, nfftA, 0.05);

% closed-form true ranking scores
It = ipdc_from_var(A, Sig, fs, nfftA); P = eye(N) - ones(N)/N;
strue = zeros(nfftA, N);
for k = 1:nfftA, R = squeeze(It(k,:,:)).^2; strue(k,:) = ((P*sum(R-R.',2))/N).'; end

fkA = oA.f >= 2 & oA.f <= fs/2 - 2;
col = lines(N);

fig = figure('Color','w','Position',[60 60 1180 470]);
tl = tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');

nexttile; hold on
for i = 1:N
    lo = oA.score_ci(:,i,1); hi = oA.score_ci(:,i,2);
    fill([oA.f(fkA); flipud(oA.f(fkA))], [lo(fkA); flipud(hi(fkA))], col(i,:), ...
         'FaceAlpha', 0.15, 'EdgeColor','none');
    plot(oA.f(fkA), oA.scores(fkA,i), '-',  'Color', col(i,:), 'LineWidth', 1.6);
    plot(oA.f(fkA), strue(fkA,i),    '--', 'Color', col(i,:), 'LineWidth', 1.0);
end
yline(0,'k:'); grid on; xlim([0 fs/2])
xlabel('Hz'); ylabel('ranking score  s_i  (high = sink)')
title('Chain 1\rightarrow2\rightarrow3: ranking \pm 95% CI (dashed = truth)')
legend({'','ch1 (source)','','','ch2','','','ch3 (sink)',''}, 'Location','best')

%% ---- B. Allen V1 <-> LM LFP ----------------------------------------------
SA = load('data/allen/subset_movie.mat');
Tmin = min(size(SA.V1_lfp,2), size(SA.LM_lfp,2));
Xe = double([SA.V1_lfp(:,1:Tmin); SA.LM_lfp(:,1:Tmin)]);
fse = double(SA.fs_lfp_v1); Ne = size(Xe,1);
NV1 = size(SA.V1_lfp,1);
[bH, aH] = butter(4, 1/(fse/2), 'high');  Xe = filtfilt(bH, aH, Xe.').';   % as in allen_lfp_analysis.m
[bN, aN] = iirnotch(60/(fse/2), 60/(fse/2)/30); Xe = filtfilt(bN, aN, Xe.').';
Xe = Xe - mean(Xe, 2);
eshort = [arrayfun(@(k) sprintf('V1-%d',k), 1:NV1, 'uni', 0), ...
          arrayfun(@(k) sprintf('LM-%d',k), 1:Ne-NV1, 'uni', 0)];
segE = 2048; ngroup = 40; bandE = [2 120];
grp  = [ones(NV1,1); 2*ones(Ne-NV1,1)];  % V1 circles, LM squares

oE = ipdc_hodge_jack(Xe, fse, NW, K, segE, 0.5, segE, 0.05, ngroup);

% band-averaged ranking scores, with jackknife CI from the divergence replicates
fkE = oE.f >= bandE(1) & oE.f <= bandE(2); Pe = eye(Ne)-ones(Ne)/Ne; njk = oE.ngrp;
sbar = mean(oE.scores(fkE,:), 1);
sJ = zeros(njk, Ne);
for g = 1:njk, dg = squeeze(mean(oE.divJ(g,fkE,:),2)); sJ(g,:) = (Pe*dg/Ne).'; end
seE = sqrt(((njk-1)/njk) * sum((sJ - mean(sJ,1)).^2, 1));

nexttile; hold on
[ss, ord] = sort(sbar);
mk = {'o', 's'}; cl = {[0.2 0.3 0.8], [0.85 0.3 0.2]};
for gi = unique(grp).'
    sel = find(grp(ord) == gi).';
    errorbar(ss(sel), sel, za*seE(ord(sel)), 'horizontal', 'LineStyle','none', ...
             'Marker', mk{gi}, 'Color', cl{gi}, 'MarkerFaceColor', cl{gi}, 'LineWidth', 1.0);
end
set(gca,'YTick',1:Ne,'YTickLabel',eshort(ord),'TickLabelInterpreter','none','FontSize',7);
xline(0,'k:'); grid on; ylim([0 Ne+1])
xlabel('ranking score  s_i  (low = source, high = sink)')
title('Allen V1/LM LFP: hierarchy \pm 95% CI (V1 circles, LM squares)')
fprintf('Allen block means: V1 %+.4f, LM %+.4f\n', mean(sbar(grp==1)), mean(sbar(grp==2)));

if ~exist('figures','dir'), mkdir('figures'); end
style_paper_figure(fig); set(fig,'PaperPositionMode','auto');
print(fig, 'figures/hodge_hierarchy.png', '-dpng', '-r150');
fprintf('SAVED figures/hodge_hierarchy.png\n');

%% summary
fprintf('\nChain: mean band-avg scores [%.3f %.3f %.3f] (expect increasing)\n', ...
        mean(oA.scores(fkA,:),1));
fprintf('Allen (%d blocks -> %d groups): score range [%.4f, %.4f]\n', ...
        oE.nseg, oE.ngrp, min(sbar), max(sbar));
fprintf('  source: %s (%.4f);  sink: %s (%.4f)\n', ...
        eshort{ord(1)}, ss(1), eshort{ord(end)}, ss(end));
