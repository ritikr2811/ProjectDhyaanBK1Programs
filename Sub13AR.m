

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

bands.delta = [0.5  4];
bands.theta = [4    8];
bands.alpha = [8   13];
bands.beta  = [13  30];
bands.gamma = [30 100];
band_names  = fieldnames(bands);
n_bands     = numel(band_names);
n_phases    = numel(phase_labels);

% Welch parameters
window_sec  = 1;
overlap_pct = 50;
nfft        = 1024;

% Permutation entropy parameters
PE_m   = 5;
PE_tau = 1;

% File matching patterns (case-insensitive prefix)
phase_patterns = {
    {'EO1','E01'}, {'EC1'}, {'G1'}, {'M1'}, ...
    {'G2'},        {'EO2','E02'}, {'EC2'}, {'M2'}
};

group_colors = [0.53 0.29 0.72;   % purple — meditator
                0.40 0.40 0.40];  % gray   — control

%% =========================================================
%  SECTION 2 — SELECT FOLDERS
%% =========================================================

fprintf('Select meditator folder...\n');
med_folder = uigetdir('', 'Select MEDITATOR folder');
if isequal(med_folder, 0), error('No folder selected.'); end

fprintf('Select control folder...\n');
ctrl_folder = uigetdir('', 'Select CONTROL folder');
if isequal(ctrl_folder, 0), error('No folder selected.'); end

fprintf('\nMeditator: %s\n', med_folder);
fprintf('Control:   %s\n\n', ctrl_folder);

%% =========================================================
%  SECTION 3 — PROBE FIRST SUBJECT FOR DIMENSIONS
%% =========================================================

tmp = probe_subject(med_folder, phase_patterns);
n_ch       = tmp.n_ch;
all_labels = tmp.labels;
fs         = tmp.fs;

window_samp  = round(window_sec * fs);
overlap_samp = round(window_samp * overlap_pct / 100);

fprintf('Channels: %d  |  fs: %d Hz  |  Welch window: %d samples\n\n', ...
    n_ch, fs, window_samp);

%% =========================================================
%  SECTION 4 — COMPUTE SpEn AND PE FOR BOTH SUBJECTS
%% =========================================================

fprintf('=== Computing meditator ===\n');
[SE_med, SE_band_med, bad_med, PE_med] = compute_spen_pe( ...
    med_folder, phase_patterns, phase_labels, ...
    bands, band_names, window_samp, overlap_samp, nfft, fs, n_ch, PE_m, PE_tau);

fprintf('\n=== Computing control ===\n');
[SE_ctrl, SE_band_ctrl, bad_ctrl, PE_ctrl] = compute_spen_pe( ...
    ctrl_folder, phase_patterns, phase_labels, ...
    bands, band_names, window_samp, overlap_samp, nfft, fs, n_ch, PE_m, PE_tau);

% NaN bad channels
SE_med(:, bad_med)          = NaN;
SE_band_med(:, bad_med, :)  = NaN;
PE_med(:, bad_med)          = NaN;
SE_ctrl(:, bad_ctrl)        = NaN;
SE_band_ctrl(:, bad_ctrl,:) = NaN;
PE_ctrl(:, bad_ctrl)        = NaN;

fprintf('\nBad channels — Meditator: [%s]  Control: [%s]\n', ...
    num2str(bad_med), num2str(bad_ctrl));

%% =========================================================
%  SECTION 5 — NORMALISE (per-subject EO1 z-score)
%  z(p, ch) = (X(p,ch) - X(EO1,ch)) / std_across_channels(EO1)

EO1_idx = 1;

% SpEn normalisation
EO1_SE_med  = SE_med(EO1_idx, :);    % [1 x n_ch]
EO1_SE_ctrl = SE_ctrl(EO1_idx, :);

sig_med_SE  = nanstd(EO1_SE_med);    % scalar std across channels
sig_ctrl_SE = nanstd(EO1_SE_ctrl);
sig_med_SE(sig_med_SE   < 1e-10) = NaN;
sig_ctrl_SE(sig_ctrl_SE < 1e-10) = NaN;

SE_med_z  = (SE_med  - EO1_SE_med)  / sig_med_SE;   % [n_phases x n_ch]
SE_ctrl_z = (SE_ctrl - EO1_SE_ctrl) / sig_ctrl_SE;

% Band SpEn normalisation
SE_band_med_z  = nan(size(SE_band_med));
SE_band_ctrl_z = nan(size(SE_band_ctrl));
for b = 1:n_bands
    EO1_b_med  = SE_band_med(EO1_idx, :, b);
    EO1_b_ctrl = SE_band_ctrl(EO1_idx, :, b);
    sig_b_med  = nanstd(EO1_b_med);
    sig_b_ctrl = nanstd(EO1_b_ctrl);
    if sig_b_med  > 1e-10
        SE_band_med_z(:,:,b)  = (SE_band_med(:,:,b)  - EO1_b_med)  / sig_b_med;
    end
    if sig_b_ctrl > 1e-10
        SE_band_ctrl_z(:,:,b) = (SE_band_ctrl(:,:,b) - EO1_b_ctrl) / sig_b_ctrl;
    end
end

% PE normalisation
EO1_PE_med  = PE_med(EO1_idx, :);
EO1_PE_ctrl = PE_ctrl(EO1_idx, :);
sig_med_PE  = nanstd(EO1_PE_med);
sig_ctrl_PE = nanstd(EO1_PE_ctrl);
sig_med_PE(sig_med_PE   < 1e-10) = NaN;
sig_ctrl_PE(sig_ctrl_PE < 1e-10) = NaN;

PE_med_z  = (PE_med  - EO1_PE_med)  / sig_med_PE;
PE_ctrl_z = (PE_ctrl - EO1_PE_ctrl) / sig_ctrl_PE;

%% =========================================================
%  SECTION 6 — FIGURE 1: SpEn TRAJECTORY
%  Channel mean +/- SEM across channels
%% =========================================================

% Channel-average per phase [n_phases x 1]
mn_SE_med  = nanmean(SE_med_z,  2);
sem_SE_med = nanstd(SE_med_z,  0,2) / sqrt(sum(~isnan(SE_med_z(1,:))));
mn_SE_ctrl = nanmean(SE_ctrl_z, 2);
sem_SE_ctrl= nanstd(SE_ctrl_z, 0,2) / sqrt(sum(~isnan(SE_ctrl_z(1,:))));

x_pos = 1:n_phases;

fig1 = figure('Name','SpEn Trajectory','NumberTitle','off','Position',[50 50 900 400]);
hold on;

fill_between(x_pos, mn_SE_med'  - sem_SE_med',  mn_SE_med'  + sem_SE_med',  group_colors(1,:));
fill_between(x_pos, mn_SE_ctrl' - sem_SE_ctrl', mn_SE_ctrl' + sem_SE_ctrl', group_colors(2,:));

plot(x_pos, mn_SE_med,  '-o', 'Color', group_colors(1,:), 'LineWidth', 2.5, ...
    'MarkerSize', 7, 'MarkerFaceColor', group_colors(1,:));
plot(x_pos, mn_SE_ctrl, '-o', 'Color', group_colors(2,:), 'LineWidth', 2.5, ...
    'MarkerSize', 7, 'MarkerFaceColor', group_colors(2,:));

yline(0, '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.9);

set(gca, 'XTick', x_pos, 'XTickLabel', phase_labels, ...
    'FontSize', 11, 'Box', 'off', 'TickDir', 'out');
xlim([0.6  n_phases+0.4]);
ylabel('SpEn z-score (EO1 ref)', 'FontSize', 11);
title('Spectral entropy across phases  (shading = SEM across channels)', ...
    'FontWeight', 'normal', 'FontSize', 12);
legend('Meditator', 'Control', 'Location', 'best', 'FontSize', 10, 'Box', 'off');
grid on;
hold off;

%% =========================================================
%  SECTION 7 — FIGURE 2: PE TRAJECTORY
%% =========================================================

mn_PE_med   = nanmean(PE_med_z,  2);
sem_PE_med  = nanstd(PE_med_z,  0,2) / sqrt(sum(~isnan(PE_med_z(1,:))));
mn_PE_ctrl  = nanmean(PE_ctrl_z, 2);
sem_PE_ctrl = nanstd(PE_ctrl_z, 0,2) / sqrt(sum(~isnan(PE_ctrl_z(1,:))));

fig2 = figure('Name','PE Trajectory','NumberTitle','off','Position',[50 480 900 400]);
hold on;

fill_between(x_pos, mn_PE_med'  - sem_PE_med',  mn_PE_med'  + sem_PE_med',  group_colors(1,:));
fill_between(x_pos, mn_PE_ctrl' - sem_PE_ctrl', mn_PE_ctrl' + sem_PE_ctrl', group_colors(2,:));

plot(x_pos, mn_PE_med,  '-o', 'Color', group_colors(1,:), 'LineWidth', 2.5, ...
    'MarkerSize', 7, 'MarkerFaceColor', group_colors(1,:));
plot(x_pos, mn_PE_ctrl, '-o', 'Color', group_colors(2,:), 'LineWidth', 2.5, ...
    'MarkerSize', 7, 'MarkerFaceColor', group_colors(2,:));

yline(0, '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.9);

set(gca, 'XTick', x_pos, 'XTickLabel', phase_labels, ...
    'FontSize', 11, 'Box', 'off', 'TickDir', 'out');
xlim([0.6  n_phases+0.4]);
ylabel('PE z-score (EO1 ref)', 'FontSize', 11);
title('Permutation entropy across phases  (shading = SEM across channels)', ...
    'FontWeight', 'normal', 'FontSize', 12);
legend('Meditator', 'Control', 'Location', 'best', 'FontSize', 10, 'Box', 'off');
grid on;
hold off;

%% =========================================================
%  SECTION 8 — FIGURE 3: SpEn TOPOMAPS
%  Row 1 = meditator z-score per phase
%  Row 2 = control   z-score per phase
%  Row 3 = difference (med - ctrl)
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

[~, data_idx, topo_idx] = intersect(all_labels, topo_ch_names, 'stable');
ch_xy = topo_xy(topo_idx, :);

fprintf('\nTopomap diagnostics:\n');
fprintf('  Channels in data:         %d\n', numel(all_labels));
fprintf('  Channels in topo layout:  %d\n', numel(topo_ch_names));
fprintf('  Matched channels:         %d\n', numel(data_idx));

if numel(data_idx) < 3
    fprintf('\n  WARNING: Fewer than 3 channels matched.\n');
    fprintf('  First 5 data labels:  %s\n', strjoin(all_labels(1:min(5,end)), ', '));
    fprintf('  First 5 topo labels:  %s\n', strjoin(topo_ch_names(1:5), ', '));
    fprintf('  Trying case-insensitive match...\n');
    % Case-insensitive fallback
    [~, data_idx, topo_idx] = intersect(lower(all_labels), lower(topo_ch_names), 'stable');
    ch_xy = topo_xy(topo_idx, :);
    fprintf('  Matched after case fix: %d\n', numel(data_idx));
end

if numel(data_idx) < 3
    warning('Not enough matched channels for topomap interpolation. Skipping Figure 3.');
    fprintf('  Your labels: %s\n', strjoin(all_labels, ', '));
else
    fprintf('  Sample matched labels: %s\n', strjoin(all_labels(data_idx(1:min(5,end))), ', '));
end

% Check actual values going into topomap
test_vals = SE_med_z(2, data_idx);   % EC1 phase, matched channels
fprintf('  SE_med_z at EC1 — min: %.4f  max: %.4f  NaN count: %d / %d\n', ...
    nanmin(test_vals), nanmax(test_vals), sum(isnan(test_vals)), numel(test_vals));
fprintf('  clim_topo: [%.4f  %.4f]\n', clim_topo(1), clim_topo(2));

res       = 300;
[gx, gy]  = meshgrid(linspace(-1.1,1.1,res), linspace(-1.1,1.1,res));
head_mask = (gx.^2 + gy.^2) > 1.02^2;
th_head   = linspace(0, 2*pi, 300);

% Colour limits
all_vals = [SE_med_z(:); SE_ctrl_z(:)];
all_vals = all_vals(isfinite(all_vals));
clim_val = prctile(abs(all_vals), 95);
if clim_val < 1e-6, clim_val = 1; end
clim_topo = [-clim_val  clim_val];

diff_z    = SE_med_z - SE_ctrl_z;   % [n_phases x n_ch]
diff_vals = diff_z(isfinite(diff_z(:)));
if isempty(diff_vals) || max(abs(diff_vals)) == 0
    clim_diff = [-0.5  0.5];
else
    clim_diff = prctile(abs(diff_vals), 95) * [-1  1];
end

% Diverging colourmap
n_c      = 256;
cmap_div = [linspace(0.22,1,n_c/2)', linspace(0.54,1,n_c/2)', linspace(0.85,1,n_c/2)'; ...
            linspace(1,0.85,n_c/2)', linspace(1,0.35,n_c/2)', linspace(1,0.19,n_c/2)'];

% Skip EO1 (all zeros by definition)
display_phases = 2:n_phases;
n_disp = numel(display_phases);

row_data   = {SE_med_z, SE_ctrl_z, diff_z};
row_titles = {'Meditator (SpEn z)', 'Control (SpEn z)', 'Difference (med - ctrl)'};
row_clims  = {clim_topo, clim_topo, clim_diff};
row_cblabs = {'z-score', 'z-score', 'Delta z'};

if numel(data_idx) < 3
    warning('Skipping topomap figure — insufficient channel matches.');
else

fig3 = figure('Name','SpEn Topomaps','NumberTitle','off', ...
    'Position',[980 50 180*n_disp 680]);

tile_h  = 0.28;
tile_w  = 1/n_disp;
margins = struct('left',0.01,'right',0.09,'top',0.07,'bottom',0.03,'gap_h',0.04);

for row = 1:3
    for col = 1:n_disp
        p_idx = display_phases(col);

        ax_l = margins.left + (col-1)*tile_w;
        ax_b = 1 - (row*(tile_h+margins.gap_h)) - margins.top + margins.gap_h;
        ax_w = tile_w - 0.005;
        ax_h = tile_h;
        ax   = axes('Position',[ax_l ax_b ax_w ax_h]); %#ok<LAXES>

        vals = row_data{row}(p_idx, data_idx);
        gz   = griddata(ch_xy(:,1), ch_xy(:,2), vals(:), gx, gy, 'v4');
        gz(head_mask) = NaN;

        imagesc(linspace(-1.1,1.1,res), linspace(-1.1,1.1,res), gz);
        set(ax,'YDir','normal');
        colormap(ax, cmap_div);
        clim(row_clims{row});
        axis equal off;
        hold(ax,'on');
        plot(ax, cos(th_head), sin(th_head), 'k-', 'LineWidth', 1.4);
        plot(ax, [-0.08 0.00 0.08], [1.00 1.10 1.00], 'k-', 'LineWidth', 1.4);
        plot(ax, [-1.00 -1.06 -1.08 -1.06 -1.00], [ 0.08  0.04  0.00 -0.04 -0.08], 'k-', 'LineWidth', 1.4);
        plot(ax, [ 1.00  1.06  1.08  1.06  1.00], [ 0.08  0.04  0.00 -0.04 -0.08], 'k-', 'LineWidth', 1.4);
        scatter(ax, ch_xy(:,1), ch_xy(:,2), 8, 'k', 'filled', 'MarkerFaceAlpha', 0.35);
        hold(ax,'off');

        if row == 1
            title(ax, phase_labels{p_idx}, 'FontSize', 9, 'FontWeight', 'normal', ...
                'Color', phase_colors(p_idx,:)*0.7);
        end
        if col == 1
            ylabel(ax, row_titles{row}, 'FontSize', 8, 'FontWeight', 'normal');
        end
        if col == n_disp
            cb = colorbar(ax, 'Location', 'eastoutside');
            cb.Label.String  = row_cblabs{row};
            cb.Label.FontSize = 8;
            cb.FontSize = 8;
        end
    end
end

annotation(fig3,'textbox',[0 0.96 1 0.04], ...
    'String','SpEn z-scored to EO1 (std across channels)', ...
    'HorizontalAlignment','center','EdgeColor','none','FontSize',12,'FontWeight','normal');

end % topomap guard

%% =========================================================
%  SECTION 9 — FIGURE 4: M2 z-SCORED TO G2 BAR CHART
%  z = (SpEn(M2,ch) - SpEn(G2,ch)) / std_G2_across_channels
%  Channel-averaged, one bar per band + broadband
%% =========================================================

M2_idx = find(strcmp(phase_labels,'M2'));
G2_idx = find(strcmp(phase_labels,'G2'));

panel_labels = [{'Broadband'}, cellfun(@(b) sprintf('%s\n%.1f-%d Hz', ...
    upper(b), bands.(b)(1), bands.(b)(2)), band_names', 'UniformOutput', false)];
n_panels = n_bands + 1;

% Compute z for broadband
G2_SE_med  = SE_med(G2_idx, :);
G2_SE_ctrl = SE_ctrl(G2_idx, :);
M2_SE_med  = SE_med(M2_idx, :);
M2_SE_ctrl = SE_ctrl(M2_idx, :);

% std of G2 across channels (per subject) as denominator
sig_G2_med  = nanstd(G2_SE_med);
sig_G2_ctrl = nanstd(G2_SE_ctrl);
sig_G2_med(sig_G2_med   < 1e-10) = NaN;
sig_G2_ctrl(sig_G2_ctrl < 1e-10) = NaN;

z_M2_med_broad  = nanmean((M2_SE_med  - G2_SE_med)  / sig_G2_med);
z_M2_ctrl_broad = nanmean((M2_SE_ctrl - G2_SE_ctrl) / sig_G2_ctrl);

% Band-level
z_M2_med_bands  = nan(1, n_bands);
z_M2_ctrl_bands = nan(1, n_bands);
for b = 1:n_bands
    G2_b_med  = SE_band_med(G2_idx, :, b);
    G2_b_ctrl = SE_band_ctrl(G2_idx, :, b);
    M2_b_med  = SE_band_med(M2_idx, :, b);
    M2_b_ctrl = SE_band_ctrl(M2_idx, :, b);

    sig_b_med  = nanstd(G2_b_med);
    sig_b_ctrl = nanstd(G2_b_ctrl);
    sig_b_med(sig_b_med   < 1e-10) = NaN;
    sig_b_ctrl(sig_b_ctrl < 1e-10) = NaN;

    z_M2_med_bands(b)  = nanmean((M2_b_med  - G2_b_med)  / sig_b_med);
    z_M2_ctrl_bands(b) = nanmean((M2_b_ctrl - G2_b_ctrl) / sig_b_ctrl);
end

% All values: [broadband, delta, theta, alpha, beta, gamma]
z_med_all  = [z_M2_med_broad,  z_M2_med_bands];
z_ctrl_all = [z_M2_ctrl_broad, z_M2_ctrl_bands];

fig4 = figure('Name','M2 z-scored to G2','NumberTitle','off','Position',[50 50 1100 420]);
tl4  = tiledlayout(1, n_panels, 'TileSpacing','compact','Padding','compact');

bar_width = 0.35;

for col = 1:n_panels
    nexttile(tl4);
    hold on;

    % Bars side by side
    b1 = bar(1, z_med_all(col),  bar_width, 'FaceColor', group_colors(1,:), ...
        'EdgeColor', group_colors(1,:)*0.7, 'LineWidth', 1.2);
    b2 = bar(2, z_ctrl_all(col), bar_width, 'FaceColor', group_colors(2,:), ...
        'EdgeColor', group_colors(2,:)*0.7, 'LineWidth', 1.2);

    yline(0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);

    set(gca, 'XTick', [1 2], 'XTickLabel', {'Med','Ctrl'}, ...
        'FontSize', 10, 'Box', 'off', 'TickDir', 'out');
    xlim([0.4 2.6]);
    title(panel_labels{col}, 'FontSize', 9, 'FontWeight', 'normal', ...
        'Interpreter', 'none');
    if col == 1
        ylabel('z-score (M2 - G2)', 'FontSize', 10);
    end
    grid on;
    hold off;
end

annotation(fig4,'textbox',[0 0.94 1 0.06], ...
    'String','SpEn: M2 z-scored to G2  |  channel-averaged', ...
    'HorizontalAlignment','center','EdgeColor','none','FontSize',12,'FontWeight','normal');

% Legend
ax4 = nexttile(tl4, 1);
hold(ax4,'on');
plot(ax4, nan, nan, 's', 'Color', group_colors(1,:), ...
    'MarkerFaceColor', group_colors(1,:), 'MarkerSize', 10, 'DisplayName', 'Meditator');
plot(ax4, nan, nan, 's', 'Color', group_colors(2,:), ...
    'MarkerFaceColor', group_colors(2,:), 'MarkerSize', 10, 'DisplayName', 'Control');
hold(ax4,'off');
legend(ax4, 'Location', 'best', 'FontSize', 9, 'Box', 'off');

%% =========================================================
%  LOCAL FUNCTIONS
%% =========================================================

function info = probe_subject(folder, phase_patterns)
    all_files = dir(fullfile(folder, '*.mat'));
    all_names = {all_files.name};
    fname = '';
    for k = 1:numel(phase_patterns{1})
        pat  = phase_patterns{1}{k};
        hits = all_names(strncmpi(all_names, pat, numel(pat)));
        if ~isempty(hits), fname = hits{1}; break; end
    end
    if isempty(fname)
        error('probe_subject: could not find EO1 file in %s', folder);
    end
    S = load(fullfile(folder, fname));
    fns = fieldnames(S);
    d = [];
    for fi = 1:numel(fns)
        if isstruct(S.(fns{fi})) && isfield(S.(fns{fi}), 'trial')
            d = S.(fns{fi}); break;
        end
    end
    if isempty(d)
        error('probe_subject: no FieldTrip struct with .trial in %s', fname);
    end
    info.n_ch   = size(d.trial{1}, 1);
    info.labels = cellstr(d.label);
    info.fs     = double(d.fsample);
end

% ----

function [SE_broad, SE_band_out, bad_elecs, PE_broad] = compute_spen_pe( ...
        folder, phase_patterns, phase_labels, bands, band_names, ...
        window_samp, overlap_samp, nfft, fs, n_ch, PE_m, PE_tau)

    n_phases    = numel(phase_labels);
    n_bands     = numel(band_names);
    SE_broad    = nan(n_phases, n_ch);
    SE_band_out = nan(n_phases, n_ch, n_bands);
    PE_broad    = nan(n_phases, n_ch);
    bad_elecs   = [];

    all_files = dir(fullfile(folder, '*.mat'));
    all_names = {all_files.name};

    % Read bad electrodes from EO1
    eo1_fname = '';
    for k = 1:numel(phase_patterns{1})
        pat  = phase_patterns{1}{k};
        hits = all_names(strncmpi(all_names, pat, numel(pat)));
        if ~isempty(hits), eo1_fname = hits{1}; break; end
    end
    if ~isempty(eo1_fname)
        S_eo1 = load(fullfile(folder, eo1_fname));
        fns_e = fieldnames(S_eo1);
        for fi = 1:numel(fns_e)
            c = S_eo1.(fns_e{fi});
            if isstruct(c) && isfield(c, 'badElecs') && ~isempty(c.badElecs)
                bad_elecs = double(c.badElecs(:)');
                break;
            end
        end
    end

    for p = 1:n_phases
        fname = '';
        for k = 1:numel(phase_patterns{p})
            pat  = phase_patterns{p}{k};
            hits = all_names(strncmpi(all_names, pat, numel(pat)));
            if ~isempty(hits), fname = hits{1}; break; end
        end
        if isempty(fname)
            warning('No file for phase %s — filled with NaN', phase_labels{p});
            continue;
        end

        S = load(fullfile(folder, fname));
        d = [];
        fns = fieldnames(S);
        for fi = 1:numel(fns)
            if isstruct(S.(fns{fi})) && isfield(S.(fns{fi}), 'trial')
                d = S.(fns{fi}); break;
            end
        end
        if isempty(d)
            warning('No FieldTrip struct in %s — skipping', fname);
            continue;
        end

        trials    = d.trial;
        n_t       = numel(trials);
        trial_len = size(trials{1}, 2);
        eff_win   = min(window_samp, trial_len);
        eff_win   = max(eff_win, 2);
        eff_lap   = round(eff_win * overlap_samp / max(window_samp,1));
        eff_nfft  = max(nfft, eff_win);
        good      = setdiff(1:n_ch, bad_elecs);

        for ci = good
            pxx_acc   = zeros(eff_nfft/2+1, 1);
            pe_trials = nan(1, n_t);

            for t = 1:n_t
                seg            = double(trials{t}(ci,:));
                [pxx_t, freqs] = pwelch(seg, hann(eff_win), eff_lap, eff_nfft, fs);
                pxx_acc        = pxx_acc + pxx_t;
                pe_trials(t)   = compute_PE(seg, PE_m, PE_tau);
            end

            pxx = pxx_acc / n_t;
            SE_broad(p, ci)  = compute_SE(pxx);
            PE_broad(p, ci)  = nanmean(pe_trials);

            for b = 1:n_bands
                brange = bands.(band_names{b});
                idx    = freqs >= brange(1) & freqs <= brange(2);
                if sum(idx) < 2, continue; end
                SE_band_out(p, ci, b) = compute_SE(pxx(idx));
            end
        end

        fprintf('  Phase %-4s done\n', phase_labels{p});
    end
end

% ----

function se = compute_SE(pxx)
    pxx = pxx(:);
    pxx(pxx < 0) = 0;
    tot = sum(pxx);
    if tot == 0 || numel(pxx) < 2, se = 0; return; end
    p = pxx / tot;
    p(p <= 0) = [];
    se = -sum(p .* log2(p)) / log2(numel(pxx));
end

% ----

function pe = compute_PE(x, m, tau)
    x          = x(:)';
    n          = length(x);
    n_patterns = factorial(m);
    max_i      = n - (m-1)*tau;
    if max_i < m, pe = NaN; return; end
    idx        = bsxfun(@plus, (1:max_i)', (0:m-1)*tau);
    embedded   = x(idx);
    [~, ranks] = sort(embedded, 2);
    pattern_ids = zeros(max_i, 1);
    for k = 1:m
        pattern_ids = pattern_ids * (m-k+1) + (ranks(:,k)-1);
    end
    counts = histcounts(pattern_ids, 0:n_patterns);
    p      = counts / sum(counts);
    p(p == 0) = [];
    pe = -sum(p .* log2(p)) / log2(n_patterns);
end

% ----

function fill_between(x, y_lo, y_hi, col)
    x_fill = [x, fliplr(x)];
    y_fill = [y_lo, fliplr(y_hi)];
    fill(x_fill, y_fill, col, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
end