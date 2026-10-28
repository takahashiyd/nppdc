%% extract_allen_subset.m -- extract the V1 + LM LFP used in the tutorial
% from Allen Brain Observatory Neuropixels session 715093703: 8 V1 and 8 LM
% channels during the first natural_movie_three block (600 s). Reads the
% session and probe NWB files (download them into data/allen/ first) and
% saves data/allen/subset_movie.mat.

clear; close all
root = fileparts(fileparts(mfilename('fullpath')));   % project root (code/)
cd(root);
sessionNWB = 'data/allen/session_715093703.nwb';
probeC_lfp = 'data/allen/probe_C_V1_lfp.nwb';   % V1 (VISp)
probeD_lfp = 'data/allen/probe_D_LM_lfp.nwb';   % LM (VISl)

% Window: natural_movie_three block #1 (600 s of stimulus-driven activity)
t_start = 2205.3;
t_stop  = 2805.8;
T_dur   = t_stop - t_start;
fprintf('Window: [%.1f, %.1f] s   duration = %.1f s\n', t_start, t_stop, T_dur);

%% ---------------- Session-level: electrodes table + units --------------
ele_id   = h5read(sessionNWB, '/general/extracellular_ephys/electrodes/id');
ele_loc  = string(h5read(sessionNWB, '/general/extracellular_ephys/electrodes/location'));
ele_pid  = h5read(sessionNWB, '/general/extracellular_ephys/electrodes/probe_id');
ele_y    = h5read(sessionNWB, '/general/extracellular_ephys/electrodes/probe_vertical_position');

unit_id        = h5read(sessionNWB, '/units/id');
unit_pkchan    = h5read(sessionNWB, '/units/peak_channel_id');
unit_quality   = string(h5read(sessionNWB, '/units/quality'));
unit_firing    = h5read(sessionNWB, '/units/firing_rate');
unit_snr       = h5read(sessionNWB, '/units/snr');
unit_isi_viol  = h5read(sessionNWB, '/units/isi_violations');

%% ---------------- Identify V1 + LM electrodes per probe ----------------
probeC_id = unique(ele_pid(ele_loc == "VISp"));    % presumed unique to probe C
probeD_id = unique(ele_pid(ele_loc == "VISl"));    % presumed unique to probe D
assert(numel(probeC_id) == 1 && numel(probeD_id) == 1, 'unexpected probe-area mapping');
fprintf('V1 probe id: %d   LM probe id: %d\n', probeC_id, probeD_id);

v1_mask = (ele_loc == "VISp") & (ele_pid == probeC_id);
lm_mask = (ele_loc == "VISl") & (ele_pid == probeD_id);
v1_ele  = ele_id(v1_mask);   v1_y = ele_y(v1_mask);
lm_ele  = ele_id(lm_mask);   lm_y = ele_y(lm_mask);
fprintf('  V1 channels in session: %d   LM channels in session: %d\n', numel(v1_ele), numel(lm_ele));

%% ---------------- Probe LFP: which channels are LFP-saved --------------
% The `electrodes` field of an LFP ElectricalSeries is a DynamicTableRegion
% that points into the *probe*-NWB electrodes table, not the session table.
% Read the probe electrodes table directly to get IDs and locations.
pC_ele_id    = h5read(probeC_lfp, '/general/extracellular_ephys/electrodes/id');
pC_ele_loc   = string(h5read(probeC_lfp, '/general/extracellular_ephys/electrodes/location'));
pC_ele_y     = h5read(probeC_lfp, '/general/extracellular_ephys/electrodes/probe_vertical_position');
pD_ele_id    = h5read(probeD_lfp, '/general/extracellular_ephys/electrodes/id');
pD_ele_loc   = string(h5read(probeD_lfp, '/general/extracellular_ephys/electrodes/location'));
pD_ele_y     = h5read(probeD_lfp, '/general/extracellular_ephys/electrodes/probe_vertical_position');

v1_lfp_idx = h5read(probeC_lfp, '/acquisition/probe_810755801_lfp/probe_810755801_lfp_data/electrodes');
lm_lfp_idx = h5read(probeD_lfp, '/acquisition/probe_810755803_lfp/probe_810755803_lfp_data/electrodes');
% Probe-NWB DynamicTableRegion: 0-based indices into the probe electrodes table.
pC_loc_per_lfp = pC_ele_loc(v1_lfp_idx + 1);
pC_id_per_lfp  = pC_ele_id(v1_lfp_idx + 1);
pC_y_per_lfp   = pC_ele_y(v1_lfp_idx + 1);
pD_loc_per_lfp = pD_ele_loc(lm_lfp_idx + 1);
pD_id_per_lfp  = pD_ele_id(lm_lfp_idx + 1);
pD_y_per_lfp   = pD_ele_y(lm_lfp_idx + 1);

v1_idx_in_lfp = find(pC_loc_per_lfp == "VISp");
lm_idx_in_lfp = find(pD_loc_per_lfp == "VISl");
fprintf('  V1 channels with LFP saved: %d   LM channels with LFP saved: %d\n', ...
        numel(v1_idx_in_lfp), numel(lm_idx_in_lfp));

% Take an evenly-spaced subset of 8 per area, across cortical depth
nKeep = 8;
v1_pick_local = round(linspace(1, numel(v1_idx_in_lfp), nKeep));
lm_pick_local = round(linspace(1, numel(lm_idx_in_lfp), nKeep));
v1_lfp_rows   = v1_idx_in_lfp(v1_pick_local);
lm_lfp_rows   = lm_idx_in_lfp(lm_pick_local);

v1_lfp_eid_real = pC_id_per_lfp(v1_lfp_rows);
lm_lfp_eid_real = pD_id_per_lfp(lm_lfp_rows);
v1_y_pick       = pC_y_per_lfp(v1_lfp_rows);
lm_y_pick       = pD_y_per_lfp(lm_lfp_rows);

%% ---------------- LFP timestamps, find window indices ------------------
v1_ts = h5read(probeC_lfp, '/acquisition/probe_810755801_lfp/probe_810755801_lfp_data/timestamps');
lm_ts = h5read(probeD_lfp, '/acquisition/probe_810755803_lfp/probe_810755803_lfp_data/timestamps');

v1_i0 = find(v1_ts >= t_start, 1, 'first');
v1_i1 = find(v1_ts <= t_stop,  1, 'last');
lm_i0 = find(lm_ts >= t_start, 1, 'first');
lm_i1 = find(lm_ts <= t_stop,  1, 'last');
nT_v1 = v1_i1 - v1_i0 + 1;
nT_lm = lm_i1 - lm_i0 + 1;
fprintf('  V1 LFP samples in window: %d   LM LFP samples in window: %d\n', nT_v1, nT_lm);

fs_v1 = 1 / median(diff(v1_ts(v1_i0:v1_i0+1000)));
fs_lm = 1 / median(diff(lm_ts(lm_i0:lm_i0+1000)));
fprintf('  LFP sampling rates: V1 = %.2f Hz, LM = %.2f Hz\n', fs_v1, fs_lm);

%% ---------------- Read LFP windows  (chunked, channel-by-channel) ------
fprintf('Reading LFP windows...\n');
V1_lfp = zeros(nKeep, nT_v1, 'single');
for k = 1:nKeep
    row = v1_lfp_rows(k);
    V1_lfp(k, :) = h5read(probeC_lfp, ...
        '/acquisition/probe_810755801_lfp/probe_810755801_lfp_data/data', ...
        [row, v1_i0], [1, nT_v1]);
end
fprintf('  V1 LFP loaded: %d x %d\n', size(V1_lfp));

LM_lfp = zeros(nKeep, nT_lm, 'single');
for k = 1:nKeep
    row = lm_lfp_rows(k);
    LM_lfp(k, :) = h5read(probeD_lfp, ...
        '/acquisition/probe_810755803_lfp/probe_810755803_lfp_data/data', ...
        [row, lm_i0], [1, nT_lm]);
end
fprintf('  LM LFP loaded: %d x %d\n', size(LM_lfp));

%% ---------------- Spikes for good-quality units in V1 and LM -----------
% Map units by peak_channel_id -> electrode area
[~, idx_in_ele] = ismember(unit_pkchan, ele_id);
unit_area = repmat("", numel(unit_pkchan), 1);
ok = idx_in_ele > 0;
unit_area(ok) = ele_loc(idx_in_ele(ok));
unit_probe = zeros(size(unit_pkchan));
unit_probe(ok) = ele_pid(idx_in_ele(ok));

good_quality = unit_quality == "good";
good_rate    = unit_firing >= 1.0 & unit_isi_viol < 0.5;
good         = good_quality & good_rate;

is_v1 = good & (unit_area == "VISp") & (unit_probe == probeC_id);
is_lm = good & (unit_area == "VISl") & (unit_probe == probeD_id);
fprintf('  V1 good units: %d   LM good units: %d\n', sum(is_v1), sum(is_lm));

% pick top-N by SNR for clean spike trains
nUnits = 10;
v1_uidx = find(is_v1);
[~, ord] = sort(unit_snr(v1_uidx), 'descend');
v1_uidx = v1_uidx(ord(1:min(nUnits, end)));
lm_uidx = find(is_lm);
[~, ord] = sort(unit_snr(lm_uidx), 'descend');
lm_uidx = lm_uidx(ord(1:min(nUnits, end)));
fprintf('  Selected V1 units: %d   LM units: %d\n', numel(v1_uidx), numel(lm_uidx));

% spike_times + spike_times_index (cumulative end indices, 0-based exclusive in NWB)
% Convert to 1-based MATLAB cumulative ends.
sti = h5read(sessionNWB, '/units/spike_times_index');   % cumulative end (NWB)
% Spikes for unit i (1-based MATLAB): spike_times(start_i : sti(i)) where
%   start_i = 1 if i==1 else sti(i-1)+1
% We only need spikes for our selected units within [t_start, t_stop].
selected_uidx = [v1_uidx(:); lm_uidx(:)];
spike_cells = cell(numel(selected_uidx), 1);
for k = 1:numel(selected_uidx)
    i = selected_uidx(k);
    if i == 1, s0 = 1; else, s0 = double(sti(i-1)) + 1; end
    s1 = double(sti(i));
    n  = s1 - s0 + 1;
    if n > 0
        st = h5read(sessionNWB, '/units/spike_times', double(s0), double(n));
        st = st(st >= t_start & st <= t_stop);
        spike_cells{k} = st;
    else
        spike_cells{k} = [];
    end
end
v1_spikes = spike_cells(1:numel(v1_uidx));
lm_spikes = spike_cells(numel(v1_uidx)+1:end);
fprintf('  V1 spike counts: '); fprintf('%d ', cellfun(@numel, v1_spikes)); fprintf('\n');
fprintf('  LM spike counts: '); fprintf('%d ', cellfun(@numel, lm_spikes)); fprintf('\n');

%% ---------------- Pack output ------------------------------------------
out = struct();
out.session_id      = 715093703;
out.area_v1         = 'VISp';
out.area_lm         = 'VISl';
out.window          = [t_start, t_stop];
out.V1_lfp          = V1_lfp;          % nKeep x nT
out.LM_lfp          = LM_lfp;
out.V1_lfp_ts       = v1_ts(v1_i0:v1_i1);
out.LM_lfp_ts       = lm_ts(lm_i0:lm_i1);
out.V1_lfp_eid      = v1_lfp_eid_real;
out.LM_lfp_eid      = lm_lfp_eid_real;
out.V1_lfp_y        = v1_y_pick;
out.LM_lfp_y        = lm_y_pick;
out.fs_lfp_v1       = fs_v1;
out.fs_lfp_lm       = fs_lm;
out.V1_spikes       = v1_spikes;
out.LM_spikes       = lm_spikes;
out.V1_unit_id      = unit_id(v1_uidx);
out.LM_unit_id      = unit_id(lm_uidx);
out.V1_unit_snr     = unit_snr(v1_uidx);
out.LM_unit_snr     = unit_snr(lm_uidx);
out.V1_unit_rate    = unit_firing(v1_uidx);
out.LM_unit_rate    = unit_firing(lm_uidx);
out.session_dataset = 'Allen Brain Observatory Visual Coding - Neuropixels';

if ~exist('data/allen', 'dir'), mkdir('data/allen'); end
out.stimulus = 'natural_movie_three, block 1';
save('data/allen/subset_movie.mat', '-struct', 'out', '-v7');
fprintf('\nSaved data/allen/subset_movie.mat (%.1f MB)\n', ...
        (dir('data/allen/subset_movie.mat').bytes)/1e6);
