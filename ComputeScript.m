

clear; clc; close all;

%% =========================================================
%  SECTION 1 — CONFIGURATION
%% =========================================================

phase_labels = {'EO1','EC1','G1','M1','G2','EO2','EC2','M2'};
phase_colors = [
    0.36 0.79 0.65;   % EO1  teal
    0.94 0.62 0.15;   % EC1  amber
    0.22 0.54 0.85;   % G1   blue
    0.85 0.35 0.19;   % M1   coral
    0.22 0.54 0.85;   % G2   blue
    0.36 0.79 0.65;   % EO2  teal
    0.94 0.62 0.15;   % EC2  amber
    0.85 0.35 0.19;   % M2   coral
];

% Frequency bands
bands.delta = [0.5  4];
bands.theta = [4    8];
bands.alpha = [8   13];
bands.beta  = [13  30];
bands.gamma = [30 100];
band_names  = fieldnames(bands);
n_bands     = numel(band_names);

% Welch parameters
window_sec   = 1;
overlap_pct  = 50;
nfft         = 1024;

% --- Permutation entropy parameters ---
% m = embedding dimension (pattern length). 3-7 typical for EEG.
%     Higher m → finer grained but needs longer signal (need m! << n_samples)
%     At 2.5s × 256Hz = 640 samples, m=5 gives 5!=120 patterns — well sampled.
% tau = time delay (samples). 1 = consecutive samples.
%     Can set tau = round(fs/f_dominant) to match a frequency of interest.
PE_m   = 5;   % embedding dimension
PE_tau = 1;   % time delay in samples

% File matching patterns (case-insensitive prefix match)
phase_patterns = {
    {'EO1','E01'}, {'EC1'}, {'G1'}, {'M1'}, ...
    {'G2'},        {'EO2','E02'}, {'EC2'}, {'M2'}
};

n_phases     = numel(phase_labels);
group_labels = {'Meditators', 'Controls'};
group_colors = [0.53 0.29 0.72;   % purple — meditators
                0.40 0.40 0.40];  % gray   — controls

%% =========================================================
%  SECTION 2 — SUBJECT IDs AND FOLDER LOADING
%% =========================================================

med_ids = { ...
    '019CKa','096MS', '040VS','012GK','095KM','090AV','054MP','015RK', ...
    '038DK', '041AG', '053DR','045SP','025RK','035SS','044PN','046ME', ...
    '031BK', '056PR', '052PR','059MS','013AR','074KS','006SR','094SR', ...
    '017KG', '042VA', '060GV','050UR','030SH','089AB'};

ctrl_ids = { ...
    '022SSP','026HM','100UK','093AK','075AD','028HB','043AK','098GS', ...
    '071GK', '077LK', '003S', '078BM','070TB','080RP','101PB','085BM', ...
    '073SK', '086AB', '082MS','064PK','102AS','084AK','062MT','049KK', ...
    '072DK', '097SV', '099SP','083SP','076BH','079SG'};

pair_gender = [ ...
    repmat({'M'}, 1, 16), ...
    repmat({'F'}, 1, 14)];

n_pairs = numel(med_ids);

% --- Select single root folder ---
fprintf('\n=== Select root folder containing ALL subject subfolders ===\n');
root_folder = uigetdir('', 'Select root folder (contains all subject subfolders)');
if isequal(root_folder, 0), error('No folder selected.'); end

all_dirs  = dir(root_folder);
all_dirs  = all_dirs([all_dirs.isdir] & ~startsWith({all_dirs.name}, '.'));
all_names = {all_dirs.name};

fprintf('Found %d subfolders in root.\n\n', numel(all_names));

% --- Verify all subject folders exist ---
fprintf('=== Verifying subject folders ===\n');
med_folders  = cell(1, n_pairs);
ctrl_folders = cell(1, n_pairs);
missing = {};

for s = 1:n_pairs
    try
        med_folders{s}  = find_subject_folder(root_folder, all_names, med_ids{s});
    catch
        missing{end+1} = sprintf('MEDITATOR  pair %d: %s', s, med_ids{s});
    end
    try
        ctrl_folders{s} = find_subject_folder(root_folder, all_names, ctrl_ids{s});
    catch
        missing{end+1} = sprintf('CONTROL    pair %d: %s', s, ctrl_ids{s});
    end
end

if ~isempty(missing)
    fprintf('\n*** Missing subject folders ***\n');
    fprintf('  %s\n', missing{:});
    error('%d subject folder(s) not found. Check root folder.', numel(missing));
end

fprintf('\n=== Confirmed pairs ===\n');
fprintf('  %-4s  %-12s  %-12s  %s\n', 'Pair','Meditator','Control','Sex');
fprintf('  %s\n', repmat('-',1,40));
for s = 1:n_pairs
    fprintf('  %-4d  %-12s  %-12s  %s\n', s, med_ids{s}, ctrl_ids{s}, pair_gender{s});
end
fprintf('\n');

%% =========================================================
%  SECTION 3 — COMPUTE SpEn FOR ALL SUBJECTS
%% =========================================================

tmp = probe_subject(struct('folder',med_folders{1},'name',med_ids{1}), phase_patterns);
n_ch       = tmp.n_ch;
all_labels = tmp.labels;
fs         = tmp.fs;

window_samp  = round(window_sec * fs);
overlap_samp = round(window_samp * overlap_pct / 100);

fprintf('Channels: %d  |  fs: %d Hz  |  Welch window: %d samples\n\n', n_ch, fs, window_samp);
fprintf('=== Computing SpEn ===\n\n');

SE_med_raw       = nan(n_pairs, n_phases, n_ch);
SE_band_med_raw  = nan(n_pairs, n_phases, n_ch, n_bands);
SE_ctrl_raw      = nan(n_pairs, n_phases, n_ch);
SE_band_ctrl_raw = nan(n_pairs, n_phases, n_ch, n_bands);
% Trial-to-trial SpEn variability (std across trials per phase per channel)
SE_var_med_raw   = nan(n_pairs, n_phases, n_ch);
SE_var_ctrl_raw  = nan(n_pairs, n_phases, n_ch);
% Permutation entropy [n_pairs x n_phases x n_ch]
PE_med_raw       = nan(n_pairs, n_phases, n_ch);
PE_ctrl_raw      = nan(n_pairs, n_phases, n_ch);
bad_elecs_med    = cell(1, n_pairs);
bad_elecs_ctrl   = cell(1, n_pairs);

for s = 1:n_pairs
    fprintf('--- Meditator %d/%d: %s ---\n', s, n_pairs, med_ids{s});
    [SE_med_raw(s,:,:), SE_band_med_raw(s,:,:,:), bad_elecs_med{s}, SE_var_med_raw(s,:,:), PE_med_raw(s,:,:)] = ...
        compute_subject_spen( ...
            struct('folder',med_folders{s},'name',med_ids{s}), ...
            phase_patterns, phase_labels, bands, band_names, ...
            window_samp, overlap_samp, nfft, fs, n_ch, PE_m, PE_tau);
end

for s = 1:n_pairs
    fprintf('--- Control %d/%d: %s ---\n', s, n_pairs, ctrl_ids{s});
    [SE_ctrl_raw(s,:,:), SE_band_ctrl_raw(s,:,:,:), bad_elecs_ctrl{s}, SE_var_ctrl_raw(s,:,:), PE_ctrl_raw(s,:,:)] = ...
        compute_subject_spen( ...
            struct('folder',ctrl_folders{s},'name',ctrl_ids{s}), ...
            phase_patterns, phase_labels, bands, band_names, ...
            window_samp, overlap_samp, nfft, fs, n_ch, PE_m, PE_tau);
end

% NaN-mask bad channels per subject (nanmean handles the rest)
for s = 1:n_pairs
    SE_med_raw(s, :, bad_elecs_med{s})         = NaN;
    SE_ctrl_raw(s, :, bad_elecs_ctrl{s})        = NaN;
    SE_band_med_raw(s, :, bad_elecs_med{s},  :) = NaN;
    SE_band_ctrl_raw(s, :, bad_elecs_ctrl{s}, :) = NaN;
    SE_var_med_raw(s, :, bad_elecs_med{s})      = NaN;
    SE_var_ctrl_raw(s, :, bad_elecs_ctrl{s})    = NaN;
    PE_med_raw(s, :, bad_elecs_med{s})          = NaN;
    PE_ctrl_raw(s, :, bad_elecs_ctrl{s})        = NaN;
end

good_ch   = 1:n_ch;
n_good    = n_ch;
good_labs = all_labels;

n_good_per_ch = sum(~isnan(squeeze(mean(SE_med_raw(:,1,:), 2))), 1) + ...
                sum(~isnan(squeeze(mean(SE_ctrl_raw(:,1,:), 2))), 1);
fprintf('\nChannel coverage:\n');
fprintf('  Min: %d / %d   Median: %.0f / %d   <50%% coverage: %d\n\n', ...
    min(n_good_per_ch), 2*n_pairs, median(n_good_per_ch), 2*n_pairs, ...
    sum(n_good_per_ch < n_pairs));

SE_med       = SE_med_raw;
SE_ctrl      = SE_ctrl_raw;
SE_band_med  = SE_band_med_raw;
SE_band_ctrl = SE_band_ctrl_raw;


EO1_idx = 1;

EO1_med  = SE_med(:, EO1_idx, :);
EO1_ctrl = SE_ctrl(:, EO1_idx, :);

SE_med_norm  = SE_med  ./ EO1_med;
SE_ctrl_norm = SE_ctrl ./ EO1_ctrl;

EO1_band_med  = SE_band_med(:, EO1_idx, :, :);
EO1_band_ctrl = SE_band_ctrl(:, EO1_idx, :, :);

SE_band_med_norm  = SE_band_med  ./ EO1_band_med;
SE_band_ctrl_norm = SE_band_ctrl ./ EO1_band_ctrl;

diff_pairs = SE_med_norm - SE_ctrl_norm;

gm_med      = squeeze(nanmean(SE_med_norm,  1));   % [n_phases x n_ch]
gm_ctrl     = squeeze(nanmean(SE_ctrl_norm, 1));
gm_diff     = squeeze(nanmean(diff_pairs,   1));
gm_diff_sem = squeeze(nanstd(diff_pairs, 0, 1)) / sqrt(n_pairs);

traj_med  = squeeze(nanmean(SE_med_norm,  3));   % [n_pairs x n_phases]
traj_ctrl = squeeze(nanmean(SE_ctrl_norm, 3));
traj_diff = squeeze(nanmean(diff_pairs,   3));

diff_band_pairs    = SE_band_med_norm - SE_band_ctrl_norm;
traj_band_med      = squeeze(nanmean(SE_band_med_norm,  3));
traj_band_ctrl     = squeeze(nanmean(SE_band_ctrl_norm, 3));
traj_band_diff     = squeeze(nanmean(diff_band_pairs,   3));

%% =========================================================
%  SECTION 5 — TOPOGRAPHIC MAP SETUP
%% =========================================================

topo_ch_names = { ...
  'Fp1','Fp2','Fz','F3','F4','F7','F8', ...
  'FC1','FC2','FC5','FC6','FT7','FT8','FT9','FT10', ...
  'C3','C4','Cz','T7','T8', ...
  'CP1','CP2','CP5','CP6','TP7','TP8','TP9','TP10', ...
  'P3','P4','Pz','P7','P8','O1','O2','Oz', ...
  'AF3','AF4','AF7','AF8','AFz', ...
  'F1','F2','F5','F6','FC3','FC4', ...
  'C1','C2','C5','C6', ...
  'CP3','CP4','CPz','P1','P2','P5','P6', ...
  'PO3','PO4','PO7','PO8','POz','Iz'};

topo_xy = [ ...
 -0.18  0.86;  0.18  0.86;  0.00  0.54; -0.33  0.58;  0.33  0.58; -0.65  0.54;  0.65  0.54; ...
 -0.17  0.35;  0.17  0.35; -0.55  0.38;  0.55  0.38; -0.76  0.35;  0.76  0.35; -0.92  0.25;  0.92  0.25; ...
 -0.51  0.00;  0.51  0.00;  0.00  0.00; -0.86  0.00;  0.86  0.00; ...
 -0.17 -0.35;  0.17 -0.35; -0.55 -0.38;  0.55 -0.38; -0.76 -0.35;  0.76 -0.35; -0.92 -0.25;  0.92 -0.25; ...
 -0.33 -0.58;  0.33 -0.58;  0.00 -0.54; -0.65 -0.54;  0.65 -0.54; -0.18 -0.86;  0.18 -0.86;  0.00 -0.90; ...
 -0.20  0.72;  0.20  0.72; -0.45  0.73;  0.45  0.73;  0.00  0.72; ...
 -0.17  0.58;  0.17  0.58; -0.50  0.56;  0.50  0.56; -0.36  0.35;  0.36  0.35; ...
 -0.26  0.00;  0.26  0.00; -0.68  0.00;  0.68  0.00; ...
 -0.36 -0.35;  0.36 -0.35;  0.00 -0.35; -0.17 -0.58;  0.17 -0.58; -0.50 -0.56;  0.50 -0.56; ...
 -0.20 -0.72;  0.20 -0.72; -0.45 -0.73;  0.45 -0.73;  0.00 -0.72;  0.00 -1.00];

[~, data_idx, topo_idx] = intersect(good_labs, topo_ch_names, 'stable');
ch_xy    = topo_xy(topo_idx, :);
res      = 300;
[gx, gy] = meshgrid(linspace(-1.1,1.1,res), linspace(-1.1,1.1,res));
head_mask = (gx.^2 + gy.^2) > 1.02^2;
th_head   = linspace(0, 2*pi, 300);

all_norm_vals = [gm_med(:); gm_ctrl(:)];
all_norm_vals = all_norm_vals(isfinite(all_norm_vals));
if numel(all_norm_vals) < 2 || range(all_norm_vals) == 0
    clim_topo = [0.9  1.1];
else
    lo = prctile(all_norm_vals, 2);
    hi = prctile(all_norm_vals, 98);
    if lo >= hi, lo = hi - 0.01; end
    clim_topo = [lo  hi];
end

diff_vals = gm_diff(isfinite(gm_diff(:)));
if isempty(diff_vals) || max(abs(diff_vals)) == 0
    clim_diff = [-0.05  0.05];
else
    clim_diff = max(abs(diff_vals)) * [-1  1];
end

n_c      = 256;
cmap_div = [linspace(0.22,1.00,n_c/2)', linspace(0.54,1.00,n_c/2)', linspace(0.85,1.00,n_c/2)'; ...
            linspace(1.00,0.85,n_c/2)', linspace(1.00,0.35,n_c/2)', linspace(1.00,0.19,n_c/2)'];

%% =========================================================
%  SECTION 6 — SAVE
%% =========================================================

save_path = fullfile(root_folder, 'SpEn_group_computed.mat');
save(save_path, ...
    'SE_med_norm','SE_ctrl_norm', ...
    'SE_band_med_norm','SE_band_ctrl_norm', ...
    'diff_pairs','diff_band_pairs', ...
    'gm_med','gm_ctrl','gm_diff','gm_diff_sem', ...
    'traj_med','traj_ctrl','traj_diff', ...
    'traj_band_med','traj_band_ctrl','traj_band_diff', ...
    'SE_med_raw','SE_ctrl_raw', ...
    'SE_band_med_raw','SE_band_ctrl_raw', ...
    'SE_var_med_raw','SE_var_ctrl_raw', ...
    'PE_med_raw','PE_ctrl_raw','PE_m','PE_tau', ...
    'all_labels','good_labs','good_ch', ...
    'n_good_per_ch','n_pairs','n_ch','phase_labels', ...
    'band_names','bands','med_ids','ctrl_ids', ...
    'bad_elecs_med','bad_elecs_ctrl','pair_gender', ...
    '-v7.3');
fprintf('\nComputation complete. Saved to:\n  %s\n', save_path);
fprintf('Now run SpEn_group_plot.m to generate figures.\n');

%% =========================================================
%  LOCAL FUNCTIONS
%% =========================================================

function fpath = find_subject_folder(root, all_names, subj_id)
    hits = all_names(strncmpi(all_names, subj_id, numel(subj_id)));
    if isempty(hits)
        error('Could not find subfolder for subject: %s', subj_id);
    end
    fpath = fullfile(root, hits{1});
end

% ---- ----

function info = probe_subject(subj, phase_patterns)
    % Reads first available phase file to discover channel count, labels, fs
    all_files = dir(fullfile(subj.folder, '*.mat'));
    all_names = {all_files.name};
    fname = '';
    for k = 1:numel(phase_patterns{1})
        pat  = phase_patterns{1}{k};
        hits = all_names(strncmpi(all_names, pat, numel(pat)));
        if ~isempty(hits), fname = hits{1}; break; end
    end
    if isempty(fname)
        error('probe_subject: could not find EO1 file in %s', subj.folder);
    end
    S = load(fullfile(subj.folder, fname));
    % Find FieldTrip struct
    fns = fieldnames(S);
    d = [];
    for fi = 1:numel(fns)
        if isstruct(S.(fns{fi})) && isfield(S.(fns{fi}), 'trial')
            d = S.(fns{fi}); break;
        end
    end
    if isempty(d)
        error('probe_subject: no FieldTrip struct with .trial found in %s', fname);
    end
    info.n_ch   = size(d.trial{1}, 1);
    info.labels = cellstr(d.label);
    info.fs     = double(d.fsample);
end

% ---- ----

function [SE_broad, SE_band_out, bad_elecs, SE_var, PE_broad] = compute_subject_spen( ...
        subj, phase_patterns, phase_labels, bands, band_names, ...
        window_samp, overlap_samp, nfft, fs, n_ch, PE_m, PE_tau)
    % -----------------------------------------------------------------
    % Computes broadband and band-limited spectral entropy for ONE
    % subject across ALL phases.
    %
    

    n_phases   = numel(phase_labels);
    n_bands    = numel(band_names);
    SE_broad   = nan(n_phases, n_ch);
    SE_band_out = nan(n_phases, n_ch, n_bands);
    SE_var     = nan(n_phases, n_ch);   % std of per-trial SpEn across trials
    PE_broad   = nan(n_phases, n_ch);   % permutation entropy (averaged across trials)
    bad_elecs  = [];   % read once from EO1; same list applied to all phases

    all_files = dir(fullfile(subj.folder, '*.mat'));
    all_names = {all_files.name};

   
    eo1_fname = '';
    for k = 1:numel(phase_patterns{1})
        pat  = phase_patterns{1}{k};
        hits = all_names(strncmpi(all_names, pat, numel(pat)));
        if ~isempty(hits), eo1_fname = hits{1}; break; end
    end
    if ~isempty(eo1_fname)
        S_eo1 = load(fullfile(subj.folder, eo1_fname));
        fns_e = fieldnames(S_eo1);
        for fi = 1:numel(fns_e)
            candidate = S_eo1.(fns_e{fi});
            if isstruct(candidate) && isfield(candidate, 'badElecs') && ~isempty(candidate.badElecs)
                bad_elecs = double(candidate.badElecs(:)');
                break;
            end
        end
    end
    if ~isempty(bad_elecs)
        fprintf('    bad_elecs (from EO1): %s\n', num2str(bad_elecs));
    end

    % ---- Phase loop — ONE loop, correctly closed at the bottom ----
    for p = 1:n_phases

        % Find file for this phase
        fname = '';
        for k = 1:numel(phase_patterns{p})
            pat  = phase_patterns{p}{k};
            hits = all_names(strncmpi(all_names, pat, numel(pat)));
            if ~isempty(hits), fname = hits{1}; break; end
        end

        if isempty(fname)
            warning('Subject %s: no file for phase %s — filled with NaN', ...
                subj.name, phase_labels{p});
            continue;   % leaves SE_broad(p,:) as NaN
        end

        % Load and locate FieldTrip struct
        S = load(fullfile(subj.folder, fname));
        d = [];
        fns = fieldnames(S);
        for fi = 1:numel(fns)
            candidate = S.(fns{fi});
            if isstruct(candidate) && isfield(candidate, 'trial')
                d = candidate; break;
            end
        end

        if isempty(d)
            warning('Subject %s phase %s: no FieldTrip struct with .trial — skipping.', ...
                subj.name, phase_labels{p});
            continue;
        end

        
        trials = d.trial;
        n_t    = numel(trials);
        good   = setdiff(1:n_ch, bad_elecs);

        % Clamp Welch window to actual trial length to avoid pwelch error
        trial_len   = size(trials{1}, 2);
        eff_win     = min(window_samp, trial_len);
        eff_win     = max(eff_win, 2);
        eff_lap     = round(eff_win * overlap_samp / max(window_samp, 1));
        eff_nfft    = max(nfft, eff_win);

        % ---- Channel loop ----
        for ci = good
            % Compute SE for each trial individually, then:
            %   mean PSD → SE_broad  (same as before, unaffected)
            %   std of per-trial SE  → SE_var  (new variability metric)
            pxx_acc   = zeros(eff_nfft/2 + 1, 1);
            se_trials = nan(1, n_t);   % per-trial broadband SpEn

            for t = 1:n_t
                seg            = double(trials{t}(ci, :));
                [pxx_t, freqs] = pwelch(seg, hann(eff_win), eff_lap, eff_nfft, fs);
                pxx_acc        = pxx_acc + pxx_t;
                se_trials(t)   = compute_SE(pxx_t);   % SE of this trial's PSD
            end
            pxx = pxx_acc / n_t;   % mean PSD across trials (unchanged)

            SE_broad(p, ci) = compute_SE(pxx);         % SE of mean PSD
            SE_var(p, ci)   = nanstd(se_trials);       % trial-to-trial variability

            % Permutation entropy — averaged across trials
            pe_trials = nan(1, n_t);
            for t = 1:n_t
                seg_pe       = double(trials{t}(ci, :));
                pe_trials(t) = compute_PE(seg_pe, PE_m, PE_tau);
            end
            PE_broad(p, ci) = nanmean(pe_trials);

            for b = 1:n_bands
                brange = bands.(band_names{b});
                idx    = freqs >= brange(1) & freqs <= brange(2);
                if sum(idx) < 2, continue; end
                SE_band_out(p, ci, b) = compute_SE(pxx(idx));
            end
        end % channel loop

        fprintf('    Phase %-4s done\n', phase_labels{p});

    end 

end % function compute_subject_spen

% ---- ----

function se = compute_SE(pxx)
    pxx = pxx(:);
    pxx(pxx < 0) = 0;
    tot = sum(pxx);
    if tot == 0 || numel(pxx) < 2, se = 0; return; end
    p = pxx / tot;
    p(p <= 0) = [];
    se = -sum(p .* log2(p)) / log2(numel(pxx));
end

% ---- ----

function fill_between(x, y_lo, y_hi, col)
    x_fill = [x, fliplr(x)];
    y_fill = [y_lo, fliplr(y_hi)];
    fill(x_fill, y_fill, col, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
end

%% =========================================================
%  LOCAL FUNCTION — compute_PE
%  ---------------------------------------------------------
%  Permutation entropy of a time series x.
%  m   = embedding dimension (ordinal pattern length)
%  tau = time delay in samples
%% =========================================================

function pe = compute_PE(x, m, tau)
    x   = x(:)';
    n   = length(x);
    n_patterns = factorial(m);
    max_i = n - (m-1)*tau;

    if max_i < m
        pe = NaN;
        return;
    end

    % Build embedded matrix [max_i x m]
    idx = bsxfun(@plus, (1:max_i)', (0:m-1)*tau);  % [max_i x m]
    embedded = x(idx);                               % [max_i x m]

    % Get ordinal pattern for each row — argsort gives rank order
    [~, ranks] = sort(embedded, 2);   % [max_i x m], each row is a permutation

    % Encode each permutation as a unique integer using factorial number system
    % This avoids building a string key and is fast
    pattern_ids = zeros(max_i, 1);
    for k = 1:m
        pattern_ids = pattern_ids * (m - k + 1) + (ranks(:,k) - 1);
    end
    % Note: this gives a unique integer per ordinal pattern (Lehmer code)

    % Count pattern frequencies
    counts = histcounts(pattern_ids, 0:n_patterns);
    p = counts / sum(counts);
    p(p == 0) = [];

    pe = -sum(p .* log2(p)) / log2(n_patterns);
end
