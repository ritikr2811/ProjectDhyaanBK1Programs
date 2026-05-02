
%% =========================================================
%  LOAD PRE-COMPUTED DATA
%% =========================================================

[data_file, data_path] = uigetfile('*.mat', ...
    'Select SpEn_group_computed.mat');
if isequal(data_file,0), error('No file selected.'); end
load(fullfile(data_path, data_file));
fprintf('Loaded: %s\n', data_file);

% Rebuild bands struct from loaded band_names
band_ranges = [0.5 4; 4 8; 8 13; 13 30; 30 100];
bands = struct();
for b = 1:numel(band_names)
    bands.(band_names{b}) = band_ranges(b,:);
end

% Restore plotting constants
phase_colors = [
    0.36 0.79 0.65; 0.94 0.62 0.15; 0.22 0.54 0.85; 0.85 0.35 0.19;
    0.22 0.54 0.85; 0.36 0.79 0.65; 0.94 0.62 0.15; 0.85 0.35 0.19];
group_colors = [0.53 0.29 0.72; 0.40 0.40 0.40];
pair_color   = [0.85 0.35 0.19];
n_bands      = numel(band_names);
n_phases     = numel(phase_labels);

fprintf('\nData loaded: %d pairs, %d phases, %d channels\n', ...
    n_pairs, n_phases, n_ch);
fprintf('Channel coverage — min: %d/%d  median: %.0f/%d subjects\n\n', ...
    min(n_good_per_ch), 2*n_pairs, median(n_good_per_ch), 2*n_pairs);

% ----------------------------------------------------------
% PLOT MODE — change this one variable to switch:
%   'raw'        absolute SpEn values (no normalisation)
%   'normalised' SpEn divided by each subject's EO1 baseline
%   'meanref'    group z-score referenced to EO1 baseline
%                z = (SpEn - mean_EO1_group) / std_EO1_group
% ----------------------------------------------------------
plot_mode = 'meanref';   % <-- change to 'raw' or 'normalised'

if strcmp(plot_mode, 'raw')
    SE_med_plot      = SE_med_raw;
    SE_ctrl_plot     = SE_ctrl_raw;
    SE_bmed_plot     = SE_band_med_raw;
    SE_bctrl_plot    = SE_band_ctrl_raw;
    y_ref            = NaN;
    y_label_traj     = 'Spectral Entropy (0-1)';
    y_label_diff     = 'SpEn diff (med - ctrl)';
    topo_title       = 'Group spectral entropy topomaps (absolute SpEn)';
    cb_label         = 'SpEn';
    fprintf('Plot mode: RAW (absolute SpEn)\n\n');

elseif strcmp(plot_mode, 'meanref')
    % Group z-score using EO1 as the reference distribution.
    %
    % For each channel, compute the mean and std of EO1 SpEn
    % ACROSS ALL SUBJECTS (both groups pooled).
    % Then z-score every subject's value at every phase using
    % those group-level EO1 statistics:
    %
    %   z(s, p, ch) = (SpEn(s,p,ch) - mean_EO1(ch)) / std_EO1(ch)
    %
    % Interpretation:
    %   z = 0   → that value equals the group's average resting baseline
    %   z = +1  → one SD above the group's resting baseline
    %   z = -1  → one SD below
    %
    
    EO1_idx = 1;

    
    EO1_all = [squeeze(SE_med_raw(:, EO1_idx, :)); ...
               squeeze(SE_ctrl_raw(:, EO1_idx, :))];   % [60 x n_ch]

    % Group std from EO1 — stable because it pools 60 subjects
    group_std_EO1 = nanstd(EO1_all, 0, 1);             % [1 x n_ch]
    group_std_EO1(group_std_EO1 < 1e-10) = NaN;

    % Per-subject EO1 values as the mean reference
    % [n_pairs x 1 x n_ch]
    EO1_med  = SE_med_raw(:,  EO1_idx, :);
    EO1_ctrl = SE_ctrl_raw(:, EO1_idx, :);

    sig = reshape(group_std_EO1, [1 1 n_ch]);  % broadcast-ready

    SE_med_plot  = (SE_med_raw  - EO1_med)  ./ sig;
    SE_ctrl_plot = (SE_ctrl_raw - EO1_ctrl) ./ sig;

   
    n_bands_local = numel(band_names);

    EO1_band_all = [squeeze(SE_band_med_raw(:, EO1_idx, :, :)); ...
                    squeeze(SE_band_ctrl_raw(:, EO1_idx, :, :))]; % [60 x n_ch x n_bands]

    group_std_bEO1 = nanstd(EO1_band_all, 0, 1);       % [1 x n_ch x n_bands]
    group_std_bEO1(group_std_bEO1 < 1e-10) = NaN;
    sig_b = reshape(group_std_bEO1, [1 1 n_ch n_bands_local]);

    % Per-subject EO1 band values [n_pairs x 1 x n_ch x n_bands]
    EO1_bmed  = SE_band_med_raw(:,  EO1_idx, :, :);
    EO1_bctrl = SE_band_ctrl_raw(:, EO1_idx, :, :);

    SE_bmed_plot  = (SE_band_med_raw  - EO1_bmed)  ./ sig_b;
    SE_bctrl_plot = (SE_band_ctrl_raw - EO1_bctrl) ./ sig_b;

    fprintf('EO1 group std (pooled %d subjects) — broadband mean: %.4f\n', ...
        size(EO1_all,1), nanmean(group_std_EO1));

    y_ref        = 0.0;
    y_label_traj = 'SpEn z (own EO1 baseline)';
    y_label_diff = 'z-score diff (med - ctrl)';
    topo_title   = 'Group SpEn topomaps (z relative to own EO1, group std)';
    cb_label     = 'z (own EO1)';
    fprintf('Plot mode: PER-SUBJECT EO1 Z-SCORE\n\n');

else
    % EO1 normalisation (original)
    SE_med_plot      = SE_med_norm;
    SE_ctrl_plot     = SE_ctrl_norm;
    SE_bmed_plot     = SE_band_med_norm;
    SE_bctrl_plot    = SE_band_ctrl_norm;
    y_ref            = 1.0;
    y_label_traj     = 'SpEn / EO1';
    y_label_diff     = 'SpEn diff (med - ctrl)';
    topo_title       = 'Group normalised spectral entropy topomaps (SpEn / EO1 baseline)';
    cb_label         = 'SpEn / EO1';
    fprintf('Plot mode: EO1-NORMALISED\n\n');
end

% Build group means and trajectories from whichever plot arrays are active
diff_plot      = SE_med_plot - SE_ctrl_plot;
gm_med         = squeeze(nanmean(SE_med_plot,   1));   % [n_phases x n_ch]
gm_ctrl        = squeeze(nanmean(SE_ctrl_plot,  1));
gm_diff        = squeeze(nanmean(diff_plot,     1));
gm_diff_sem    = squeeze(nanstd(diff_plot, 0,1)) / sqrt(n_pairs);
traj_med       = squeeze(nanmean(SE_med_plot,   3));   % [n_pairs x n_phases]
traj_ctrl      = squeeze(nanmean(SE_ctrl_plot,  3));
traj_diff      = squeeze(nanmean(diff_plot,     3));
diff_band_plot = SE_bmed_plot - SE_bctrl_plot;
traj_band_med  = squeeze(nanmean(SE_bmed_plot,   3));  % [n_pairs x n_phases x n_bands]
traj_band_ctrl = squeeze(nanmean(SE_bctrl_plot,  3));
traj_band_diff = squeeze(nanmean(diff_band_plot, 3));

%% =========================================================
%  TOPOGRAPHIC MAP SETUP
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
ch_xy     = topo_xy(topo_idx, :);
res       = 300;
[gx, gy]  = meshgrid(linspace(-1.1,1.1,res), linspace(-1.1,1.1,res));
head_mask = (gx.^2 + gy.^2) > 1.02^2;
th_head   = linspace(0, 2*pi, 300);


winsor_pct = [5  95];   

% SE_med_plot / SE_ctrl_plot already set correctly for whichever mode is active
all_subj_vals = [SE_med_plot(:); SE_ctrl_plot(:)];
all_subj_vals = all_subj_vals(isfinite(all_subj_vals));
if numel(all_subj_vals) > 10
    w_lo = prctile(all_subj_vals, winsor_pct(1));
    w_hi = prctile(all_subj_vals, winsor_pct(2));
else
    % Sensible fallbacks per mode
    if strcmp(plot_mode, 'meanref')
        w_lo = -2.0;  w_hi = 2.0;   % z-score: ±2 SD covers ~95% of data
    elseif strcmp(plot_mode, 'raw')
        w_lo = 0.5;   w_hi = 1.0;
    else
        w_lo = 0.8;   w_hi = 1.2;
    end
end
fprintf('Winsorise limits: [%.3f  %.3f]  (%d-%d pct)\n', w_lo, w_hi, winsor_pct(1), winsor_pct(2));

% Clamp group mean maps for display (saved data untouched)
gm_med_disp  = min(max(gm_med,  w_lo), w_hi);
gm_ctrl_disp = min(max(gm_ctrl, w_lo), w_hi);
gm_diff_disp = gm_med_disp - gm_ctrl_disp;

% Colour limits from winsorised group means
all_norm_vals = [gm_med_disp(:); gm_ctrl_disp(:)];
all_norm_vals = all_norm_vals(isfinite(all_norm_vals));
if numel(all_norm_vals) < 2 || range(all_norm_vals) == 0
    clim_topo = [0.5  1.0];
else
    lo = prctile(all_norm_vals, 2);
    hi = prctile(all_norm_vals, 98);
    if lo >= hi, lo = hi - 0.01; end
    clim_topo = [lo  hi];
end

diff_vals = gm_diff_disp(isfinite(gm_diff_disp(:)));
if isempty(diff_vals) || max(abs(diff_vals)) == 0
    if strcmp(plot_mode, 'meanref')
        clim_diff = [-1.0  1.0];   % z-score difference fallback
    else
        clim_diff = [-0.05  0.05];
    end
else
    clim_diff = max(abs(diff_vals)) * [-1  1];
end

% Diverging colourmap for difference row
n_c      = 256;
cmap_div = [linspace(0.22,1.00,n_c/2)', linspace(0.54,1.00,n_c/2)', linspace(0.85,1.00,n_c/2)'; ...
            linspace(1.00,0.85,n_c/2)', linspace(1.00,0.35,n_c/2)', linspace(1.00,0.19,n_c/2)'];

%% =========================================================
%  FIGURE 1: TOPO GRID
%  3 rows × 7 columns (EO1 skipped — it is the baseline)
%  Row 1 = meditators group mean normalised SpEn
%  Row 2 = controls group mean normalised SpEn
%  Row 3 = difference (med - ctrl), diverging colourmap
%% =========================================================

display_phases = 2:n_phases;   % skip EO1 (all ≈ 1.0 by definition)
n_disp = numel(display_phases);

row_data   = {gm_med_disp, gm_ctrl_disp, gm_diff_disp};  % winsorised display copies
row_titles = {'Meditators (norm. SpEn)', 'Controls (norm. SpEn)', 'Paired diff (med - ctrl)'};
row_cmaps  = {hot, hot, cmap_div};
row_clims  = {clim_topo, clim_topo, clim_diff};
row_cblabs = {'SpEn / EO1', 'SpEn / EO1', 'Delta SpEn'};

fig1 = figure('Name','Normalised SpEn Topomaps — All Phases', ...
    'NumberTitle','off', 'Position',[20 20 200*n_disp 680]);

tile_h  = 0.28;
tile_w  = 1/n_disp;
margins = struct('left',0.01,'right',0.08,'top',0.06,'bottom',0.04,'gap_h',0.04);

for row = 1:3
    data_row = row_data{row};   % [n_phases x n_good]

    for col = 1:n_disp
        p_idx = display_phases(col);

        ax_l = margins.left + (col-1)*tile_w;
        ax_b = 1 - (row * (tile_h + margins.gap_h)) - margins.top + margins.gap_h;
        ax_w = tile_w - 0.005;
        ax_h = tile_h;
        ax   = axes('Position',[ax_l  ax_b  ax_w  ax_h]); %#ok<LAXES>

        vals = data_row(p_idx, data_idx);
        gz   = griddata(ch_xy(:,1), ch_xy(:,2), vals(:), gx, gy, 'v4');
        gz(head_mask) = NaN;

        imagesc(linspace(-1.1,1.1,res), linspace(-1.1,1.1,res), gz);
        set(ax,'YDir','normal');
        colormap(ax, row_cmaps{row});
        clim(row_clims{row});
        axis equal off;
        hold(ax,'on');

        % Head outline
        plot(ax, cos(th_head), sin(th_head), 'k-', 'LineWidth',1.4);
        % Nose
        plot(ax,[-0.08  0.00  0.08],[1.00  1.10  1.00],'k-','LineWidth',1.4);
        % Ears
        plot(ax,[-1.00 -1.06 -1.08 -1.06 -1.00],[ 0.08  0.04  0.00 -0.04 -0.08],'k-','LineWidth',1.4);
        plot(ax,[ 1.00  1.06  1.08  1.06  1.00],[ 0.08  0.04  0.00 -0.04 -0.08],'k-','LineWidth',1.4);
        % Electrode dots
        scatter(ax, ch_xy(:,1), ch_xy(:,2), 8, 'k', 'filled', 'MarkerFaceAlpha',0.35);

        hold(ax,'off');

        if row == 1
            title(ax, phase_labels{p_idx}, 'FontSize',10, 'FontWeight','normal', ...
                'Color', phase_colors(p_idx,:)*0.7);
        end
        if col == 1
            ylabel(ax, row_titles{row}, 'FontSize',9, 'FontWeight','normal');
        end
        if col == n_disp
            cb = colorbar(ax, 'Location','eastoutside');
            cb.Label.String  = strrep(row_cblabs{row}, 'SpEn / EO1', cb_label);
            cb.Label.FontSize = 8;
            cb.FontSize = 8;
        end
    end
end

annotation(fig1,'textbox',[0 0.97 1 0.03], ...
    'String',topo_title, ...
    'HorizontalAlignment','center','EdgeColor','none','FontSize',13,'FontWeight','normal');

%% =========================================================
%  FIGURE 2: PHASE TRAJECTORY
%  Row 1 = meditators vs controls (normalised SpEn)
%  Row 2 = paired difference (med - ctrl)
%  One panel per band + broadband
%% =========================================================

n_panels = n_bands + 1;
x_pos    = 1:n_phases;

fig2 = figure('Name','Normalised SpEn Phase Trajectory', ...
    'NumberTitle','off', 'Position',[20 740 1300 520]);
tl = tiledlayout(fig2, 2, n_panels, 'TileSpacing','compact', 'Padding','compact');

for row_t = 1:2
    for pan = 1:n_panels

        nexttile(tl);
        hold on;

        if pan <= n_bands
            b = pan;
            y_med   = squeeze(traj_band_med(:,:,b));
            y_ctrl  = squeeze(traj_band_ctrl(:,:,b));
            y_diff  = squeeze(traj_band_diff(:,:,b));
            % Fix label for delta (0.5 Hz shows as 5e-01 with %d)
            blo = bands.(band_names{b})(1);
            bhi = bands.(band_names{b})(2);
            if blo < 1
                panel_title = sprintf('%s\n%.1f-%d Hz', upper(band_names{b}), blo, bhi);
            else
                panel_title = sprintf('%s\n%d-%d Hz', upper(band_names{b}), blo, bhi);
            end
        else
            y_med   = traj_med;
            y_ctrl  = traj_ctrl;
            y_diff  = traj_diff;
            panel_title = 'Broadband';
        end

        % Winsorise trajectory data for display (same limits as topo)
        y_med  = min(max(y_med,  w_lo), w_hi);
        y_ctrl = min(max(y_ctrl, w_lo), w_hi);
        y_diff = min(max(y_diff, -(w_hi-1)*2), (w_hi-1)*2);

        if row_t == 1
            % Spaghetti
            for s = 1:n_pairs
                plot(x_pos, y_med(s,:),  '-', 'Color', [group_colors(1,:) 0.20], 'LineWidth',0.7);
                plot(x_pos, y_ctrl(s,:), '-', 'Color', [group_colors(2,:) 0.20], 'LineWidth',0.7);
            end
            % Group means + SEM bands
            mn_med   = nanmean(y_med,  1);
            sem_med  = nanstd(y_med,  0,1) / sqrt(sum(~isnan(y_med(:,1))));
            mn_ctrl  = nanmean(y_ctrl, 1);
            sem_ctrl = nanstd(y_ctrl, 0,1) / sqrt(sum(~isnan(y_ctrl(:,1))));
            fill_between(x_pos, mn_med  - sem_med,  mn_med  + sem_med,  group_colors(1,:));
            fill_between(x_pos, mn_ctrl - sem_ctrl, mn_ctrl + sem_ctrl, group_colors(2,:));
            plot(x_pos, mn_med,  '-o', 'Color', group_colors(1,:), ...
                'LineWidth',2.5,'MarkerSize',6,'MarkerFaceColor',group_colors(1,:));
            plot(x_pos, mn_ctrl, '-o', 'Color', group_colors(2,:), ...
                'LineWidth',2.5,'MarkerSize',6,'MarkerFaceColor',group_colors(2,:));
            if ~isnan(y_ref), yline(y_ref,'--','Color',[0.75 0.75 0.75],'LineWidth',0.9); end
            ylabel(y_label_traj,'FontSize',9);
            title(panel_title,'FontWeight','normal','FontSize',10,'Interpreter','none');

        else
            % Paired difference spaghetti
            for s = 1:n_pairs
                plot(x_pos, y_diff(s,:), '-', 'Color', [pair_color 0.22], 'LineWidth',0.7);
            end
            mn_diff  = nanmean(y_diff, 1);
            sem_diff = nanstd(y_diff, 0,1) / sqrt(sum(~isnan(y_diff(:,1))));
            fill_between(x_pos, mn_diff - sem_diff, mn_diff + sem_diff, pair_color);
            plot(x_pos, mn_diff, '-o', 'Color', pair_color, ...
                'LineWidth',2.5,'MarkerSize',6,'MarkerFaceColor',pair_color);
            yline(0,'--','Color',[0.75 0.75 0.75],'LineWidth',0.9);
            ylabel(y_label_diff,'FontSize',9);
        end

        set(gca,'XTick',x_pos,'XTickLabel',phase_labels, ...
            'FontSize',9,'Box','off','TickDir','out');
        xlim([0.6  n_phases+0.4]);
        grid on;
        hold off;
    end
end

% Legend on first tile
ax_first = nexttile(tl, 1);
hold(ax_first,'on');
h_leg = [ ...
    plot(ax_first,nan,nan,'-o','Color',group_colors(1,:),'LineWidth',2.5, ...
         'MarkerFaceColor',group_colors(1,:),'MarkerSize',6, ...
         'DisplayName',sprintf('Meditators (n=%d)', n_pairs)); ...
    plot(ax_first,nan,nan,'-o','Color',group_colors(2,:),'LineWidth',2.5, ...
         'MarkerFaceColor',group_colors(2,:),'MarkerSize',6, ...
         'DisplayName',sprintf('Controls (n=%d)', n_pairs)); ...
    plot(ax_first,nan,nan,'-o','Color',pair_color,'LineWidth',2.5, ...
         'MarkerFaceColor',pair_color,'MarkerSize',6, ...
         'DisplayName','Paired diff (med - ctrl)')];
hold(ax_first,'off');
legend(ax_first, h_leg, 'Location','northwest','FontSize',9,'Box','off');

annotation(fig2,'textbox',[0 0.97 1 0.03], ...
    'String','Normalised spectral entropy across phases (shading = SEM)', ...
    'HorizontalAlignment','center','EdgeColor','none', ...
    'FontSize',13,'FontWeight','normal');

%% =========================================================
%  FIGURE 3: M2 vs G2 TOPOMAP
%  z = (SpEn(M2) - SpEn(G2)) / std_G2_group(ch)
%  One row per group + difference row
%  One column per band + broadband
%% =========================================================

M2_idx = find(strcmp(phase_labels, 'M2'));
G2_idx = find(strcmp(phase_labels, 'G2'));

% --- Broadband z-score: M2 referenced to G2 ---
% Pool G2 across both groups to get stable std per channel
G2_all_broad = [squeeze(SE_med_raw(:, G2_idx, :)); ...
                squeeze(SE_ctrl_raw(:, G2_idx, :))];   % [60 x n_ch]
std_G2_broad = nanstd(G2_all_broad, 0, 1);             % [1 x n_ch]
std_G2_broad(std_G2_broad < 1e-10) = NaN;

% Per-subject G2 subtraction, group std scaling
G2_med_broad  = squeeze(SE_med_raw(:,  G2_idx, :));   % [n_pairs x n_ch]
G2_ctrl_broad = squeeze(SE_ctrl_raw(:, G2_idx, :));
M2_med_broad  = squeeze(SE_med_raw(:,  M2_idx, :));
M2_ctrl_broad = squeeze(SE_ctrl_raw(:, M2_idx, :));

z_M2_med_broad  = (M2_med_broad  - G2_med_broad)  ./ std_G2_broad;  % [n_pairs x n_ch]
z_M2_ctrl_broad = (M2_ctrl_broad - G2_ctrl_broad) ./ std_G2_broad;
z_M2_diff_broad = z_M2_med_broad - z_M2_ctrl_broad;

% --- Band-level z-score ---
n_bands_local = numel(band_names);
z_M2_med_band  = nan(n_pairs, n_ch, n_bands_local);
z_M2_ctrl_band = nan(n_pairs, n_ch, n_bands_local);
z_M2_diff_band = nan(n_pairs, n_ch, n_bands_local);

for b = 1:n_bands_local
    G2_band_all = [squeeze(SE_band_med_raw(:,  G2_idx, :, b)); ...
                   squeeze(SE_band_ctrl_raw(:, G2_idx, :, b))];  % [60 x n_ch]
    std_G2_b = nanstd(G2_band_all, 0, 1);
    std_G2_b(std_G2_b < 1e-10) = NaN;

    G2_med_b  = squeeze(SE_band_med_raw(:,  G2_idx, :, b));
    G2_ctrl_b = squeeze(SE_band_ctrl_raw(:, G2_idx, :, b));
    M2_med_b  = squeeze(SE_band_med_raw(:,  M2_idx, :, b));
    M2_ctrl_b = squeeze(SE_band_ctrl_raw(:, M2_idx, :, b));

    z_M2_med_band(:,:,b)  = (M2_med_b  - G2_med_b)  ./ std_G2_b;
    z_M2_ctrl_band(:,:,b) = (M2_ctrl_b - G2_ctrl_b) ./ std_G2_b;
    z_M2_diff_band(:,:,b) = z_M2_med_band(:,:,b) - z_M2_ctrl_band(:,:,b);
end

% --- Group means [n_ch] for topomap ---
gm_M2_med  = [nanmean(z_M2_med_broad,  1); ...   % [n_bands+1 x n_ch]
              squeeze(nanmean(z_M2_med_band,  1))'];
gm_M2_ctrl = [nanmean(z_M2_ctrl_broad, 1); ...
              squeeze(nanmean(z_M2_ctrl_band, 1))'];
gm_M2_diff = [nanmean(z_M2_diff_broad, 1); ...
              squeeze(nanmean(z_M2_diff_band, 1))'];

% Panel labels
panel_labels = [{'Broadband'}, cellfun(@(b) sprintf('%s\n%.1f-%d Hz', ...
    upper(b), bands.(b)(1), bands.(b)(2)), band_names', 'UniformOutput', false)];
n_panels_f3 = n_bands_local + 1;

% --- Draw figure: bar + dot plot, one panel per band + broadband ---
fig3 = figure('Name','M2 z-scored to G2', ...
    'NumberTitle','off', 'Position',[100 100 1200 420]);

tl3 = tiledlayout(1, n_panels_f3, 'TileSpacing','compact', 'Padding','compact');

for col = 1:n_panels_f3
    nexttile(tl3);
    hold on;

    % Channel-averaged z-score per subject [n_pairs x 1]
    if col == 1
        z_med_subj  = nanmean(z_M2_med_broad,  2);   % broadband
        z_ctrl_subj = nanmean(z_M2_ctrl_broad, 2);
    else
        z_med_subj  = nanmean(z_M2_med_band(:,:,col-1),  2);
        z_ctrl_subj = nanmean(z_M2_ctrl_band(:,:,col-1), 2);
    end

    mn_med  = nanmean(z_med_subj);
    sem_med = nanstd(z_med_subj)  / sqrt(sum(~isnan(z_med_subj)));
    mn_ctrl = nanmean(z_ctrl_subj);
    sem_ctrl= nanstd(z_ctrl_subj) / sqrt(sum(~isnan(z_ctrl_subj)));

    % Individual subject dots with jitter
    jitter = (rand(n_pairs,1) - 0.5) * 0.15;
    scatter(ones(n_pairs,1) + jitter,      z_med_subj,  28, group_colors(1,:), ...
        'filled', 'MarkerFaceAlpha', 0.45);
    scatter(2*ones(n_pairs,1) + jitter,    z_ctrl_subj, 28, group_colors(2,:), ...
        'filled', 'MarkerFaceAlpha', 0.45);

    % Paired lines connecting each meditator to their control
    for s = 1:n_pairs
        plot([1 2] + jitter(s), [z_med_subj(s) z_ctrl_subj(s)], ...
            '-', 'Color', [0.6 0.6 0.6 0.25], 'LineWidth', 0.7);
    end

    % Mean ± SEM bars
    errorbar(1, mn_med,  sem_med,  'o', 'Color', group_colors(1,:), ...
        'MarkerFaceColor', group_colors(1,:), 'MarkerSize', 8, 'LineWidth', 2.0, ...
        'CapSize', 8);
    errorbar(2, mn_ctrl, sem_ctrl, 'o', 'Color', group_colors(2,:), ...
        'MarkerFaceColor', group_colors(2,:), 'MarkerSize', 8, 'LineWidth', 2.0, ...
        'CapSize', 8);

    yline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.9);

    set(gca, 'XTick', [1 2], 'XTickLabel', {'Med','Ctrl'}, ...
        'FontSize', 10, 'Box', 'off', 'TickDir', 'out');
    xlim([0.5 2.5]);
    title(panel_labels{col}, 'FontSize', 10, 'FontWeight', 'normal', ...
        'Interpreter', 'none');
    if col == 1
        ylabel('z-score (M2 − G2)', 'FontSize', 10);
    end
    grid on;
    hold off;
end

annotation(fig3,'textbox',[0 0.95 1 0.05], ...
    'String', sprintf('M2 z-scored to G2  |  channel-averaged  |  n=%d pairs', n_pairs), ...
    'HorizontalAlignment','center','EdgeColor','none','FontSize',13,'FontWeight','normal');

% Legend
ax3_first = nexttile(tl3, 1);
hold(ax3_first, 'on');
h3 = [plot(ax3_first, nan, nan, 'o', 'Color', group_colors(1,:), ...
           'MarkerFaceColor', group_colors(1,:), 'MarkerSize', 7, ...
           'DisplayName', sprintf('Meditators (n=%d)', n_pairs)); ...
      plot(ax3_first, nan, nan, 'o', 'Color', group_colors(2,:), ...
           'MarkerFaceColor', group_colors(2,:), 'MarkerSize', 7, ...
           'DisplayName', sprintf('Controls (n=%d)', n_pairs))];
hold(ax3_first, 'off');
legend(ax3_first, h3, 'Location', 'best', 'FontSize', 9, 'Box', 'off');


%% =========================================================
%  FIGURE 7: PERMUTATION ENTROPY TRAJECTORY
%  Same layout as Figure 2 top row but for PE instead of SpEn
%  Z-scored using EO1 group baseline (same method as SpEn)
%% =========================================================

% Per-subject EO1 z-score for PE — matches SpEn normalisation:
%   z(s,p,ch) = (PE(s,p,ch) - PE(s,EO1,ch)) / std_EO1_group(ch)
EO1_PE_all  = [squeeze(PE_med_raw(:,1,:)); squeeze(PE_ctrl_raw(:,1,:))];  % [60 x n_ch]
PE_sig_EO1  = reshape(nanstd(EO1_PE_all, 0, 1), [1 1 n_ch]);  % group std only
PE_sig_EO1(PE_sig_EO1 < 1e-10) = NaN;
EO1_PE_med  = PE_med_raw(:,  1, :);   % each subject's own EO1 [n_pairs x 1 x n_ch]
EO1_PE_ctrl = PE_ctrl_raw(:, 1, :);

PE_med_z    = (PE_med_raw  - EO1_PE_med)  ./ PE_sig_EO1;  % [n_pairs x n_phases x n_ch]
PE_ctrl_z   = (PE_ctrl_raw - EO1_PE_ctrl) ./ PE_sig_EO1;
PE_diff_z   = PE_med_z - PE_ctrl_z;

% Channel-average trajectories [n_pairs x n_phases]
traj_PE_med  = squeeze(nanmean(PE_med_z,  3));
traj_PE_ctrl = squeeze(nanmean(PE_ctrl_z, 3));
traj_PE_diff = squeeze(nanmean(PE_diff_z, 3));

% -- Dot plot: one panel per phase, meditators vs controls --
fig7 = figure('Name','Permutation Entropy — EO1 z-score by phase', ...
    'NumberTitle','off', 'Position',[100 100 1300 440]);

tl7 = tiledlayout(1, n_phases, 'TileSpacing','compact', 'Padding','compact');

for p = 1:n_phases
    nexttile(tl7);
    hold on;

    v_med  = traj_PE_med(:,  p);   % [n_pairs x 1]
    v_ctrl = traj_PE_ctrl(:, p);

    mn_med   = nanmean(v_med);
    sem_med  = nanstd(v_med)  / sqrt(sum(~isnan(v_med)));
    mn_ctrl  = nanmean(v_ctrl);
    sem_ctrl = nanstd(v_ctrl) / sqrt(sum(~isnan(v_ctrl)));

    % Paired lines
    jitter = (rand(n_pairs,1) - 0.5) * 0.15;
    for s = 1:n_pairs
        plot([1 2] + jitter(s), [v_med(s) v_ctrl(s)], ...
            '-', 'Color', [0.6 0.6 0.6 0.25], 'LineWidth', 0.7);
    end

    % Individual dots
    scatter(ones(n_pairs,1)   + jitter, v_med,  26, group_colors(1,:), ...
        'filled', 'MarkerFaceAlpha', 0.5);
    scatter(2*ones(n_pairs,1) + jitter, v_ctrl, 26, group_colors(2,:), ...
        'filled', 'MarkerFaceAlpha', 0.5);

    % Mean +/- SEM
    errorbar(1, mn_med,  sem_med,  'o', 'Color', group_colors(1,:), ...
        'MarkerFaceColor', group_colors(1,:), 'MarkerSize', 8, ...
        'LineWidth', 2, 'CapSize', 7);
    errorbar(2, mn_ctrl, sem_ctrl, 'o', 'Color', group_colors(2,:), ...
        'MarkerFaceColor', group_colors(2,:), 'MarkerSize', 8, ...
        'LineWidth', 2, 'CapSize', 7);

    yline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.9);

    set(gca, 'XTick', [1 2], 'XTickLabel', {'Med','Ctrl'}, ...
        'FontSize', 9, 'Box', 'off', 'TickDir', 'out');
    xlim([0.5 2.5]);
    title(phase_labels{p}, 'FontSize', 10, 'FontWeight', 'normal', ...
        'Color', phase_colors(p,:) * 0.7);
    if p == 1
        ylabel('PE z-score (EO1 ref)', 'FontSize', 10);
    end
    grid on;
    hold off;
end

annotation(fig7, 'textbox', [0 0.95 1 0.05], ...
    'String', sprintf('Permutation entropy z-scored to EO1  |  n=%d pairs', n_pairs), ...
    'HorizontalAlignment', 'center', 'EdgeColor', 'none', ...
    'FontSize', 13, 'FontWeight', 'normal');

% Legend on first tile
ax7 = nexttile(tl7, 1);
hold(ax7, 'on');
h7 = [plot(ax7, nan, nan, 'o', 'Color', group_colors(1,:), ...
           'MarkerFaceColor', group_colors(1,:), 'MarkerSize', 7, ...
           'DisplayName', sprintf('Meditators (n=%d)', n_pairs)); ...
      plot(ax7, nan, nan, 'o', 'Color', group_colors(2,:), ...
           'MarkerFaceColor', group_colors(2,:), 'MarkerSize', 7, ...
           'DisplayName', sprintf('Controls (n=%d)', n_pairs))];
hold(ax7, 'off');
legend(ax7, h7, 'Location', 'best', 'FontSize', 9, 'Box', 'off');

%% =========================================================
%  FIGURE 8: SpEn vs PE SCATTER — M2 only
%  Each dot = one subject, colour = group
%  Shows whether SpEn and PE agree or dissociate at M2
%% =========================================================

M2_idx_pe = find(strcmp(phase_labels,'M2'));

% Channel-averaged M2 z-scores
SpEn_M2_med  = nanmean(squeeze(SE_med_plot(:, M2_idx_pe, :)),  2);  % [n_pairs x 1]
SpEn_M2_ctrl = nanmean(squeeze(SE_ctrl_plot(:, M2_idx_pe, :)), 2);
PE_M2_med    = nanmean(squeeze(PE_med_z(:,  M2_idx_pe, :)), 2);
PE_M2_ctrl   = nanmean(squeeze(PE_ctrl_z(:, M2_idx_pe, :)), 2);

fig8 = figure('Name','SpEn vs PE at M2','NumberTitle','off','Position',[100 620 520 450]);
hold on;
scatter(SpEn_M2_med,  PE_M2_med,  55, group_colors(1,:), 'filled','MarkerFaceAlpha',0.7);
scatter(SpEn_M2_ctrl, PE_M2_ctrl, 55, group_colors(2,:), 'filled','MarkerFaceAlpha',0.7);

% Regression lines per group
all_x = [SpEn_M2_med;  SpEn_M2_ctrl];
all_y = [PE_M2_med;    PE_M2_ctrl];
finite_mask = isfinite(all_x) & isfinite(all_y);
if sum(finite_mask) > 3
    p_fit = polyfit(all_x(finite_mask), all_y(finite_mask), 1);
    x_fit = linspace(min(all_x(finite_mask)), max(all_x(finite_mask)), 100);
    plot(x_fit, polyval(p_fit, x_fit), 'k--', 'LineWidth', 1.2);
    [r, pval] = corr(all_x(finite_mask), all_y(finite_mask), 'Rows','complete');
    text(0.05, 0.95, sprintf('r = %.2f,  p = %.3f', r, pval), ...
        'Units','normalized','FontSize',9,'VerticalAlignment','top');
end

xline(0,'--','Color',[0.8 0.8 0.8],'LineWidth',0.8);
yline(0,'--','Color',[0.8 0.8 0.8],'LineWidth',0.8);
xlabel('SpEn z-score at M2','FontSize',11);
ylabel('PE z-score at M2','FontSize',11);
title('Spectral entropy vs permutation entropy at M2',...
    'FontWeight','normal','FontSize',12);
legend(sprintf('Meditators (n=%d)',n_pairs), sprintf('Controls (n=%d)',n_pairs),...
    'Location','best','FontSize',9,'Box','off');
set(gca,'FontSize',10,'Box','off','TickDir','out'); grid on;
hold off;

%% =========================================================
%  LOCAL FUNCTIONS
%% =========================================================

function fill_between(x, y_lo, y_hi, col)
    x_fill = [x, fliplr(x)];
    y_fill = [y_lo, fliplr(y_hi)];
    fill(x_fill, y_fill, col, 'FaceAlpha',0.15, 'EdgeColor','none');
end