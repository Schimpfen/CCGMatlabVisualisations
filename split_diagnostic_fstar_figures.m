% diagnostic_fstar_figures.m
% Pooled and per-condition separation S(f) for log10|Z| and phase.
% - Uses ALL *.xlsx in cwd (exports are auto-skipped)
% - One figure per electrode (E): pooled + each-condition overlaid in same axes
% - One collector figure with all electrodes tiled (rows) for both metrics (2 cols)
%
% Filename pattern (preferred): E<digits>_<condition>_*.xlsx
% Data columns: [Re(Z) Ohm, Im(Z) Ohm, (optional) Frequency Hz]
% If no frequency column, script reads "Start Frequency"/"Stop Frequency" to synthesize.

clear; clc;

%% ---------------- User settings ----------------
dataDir         = pwd;
CONDS           = [];                % [] = accept all discovered conditions
BASELINE        = 'cleaned';         % preferred baseline label (case-insensitive)
F_TRIM_MIN_HZ   = 1e4;               % [] disables
F_TRIM_MAX_HZ   = 1.5e5;             % [] disables
NGRID           = 300;

% Plot style
fontName      = 'Helvetica';
fontSizeAx    = 9;  fontSizeLbl = 10;  fontSizeTtl = 10;
exportDPI     = 300;  doPDFexport = true;
pubWidthOne   = 7.4;  pubHeightOne = 3.2;   % per-electrode figure (two axes)
pubWidthGrid  = 7.4;  rowHeight    = 2.4;   % collector

figOutPrefix  = fullfile(dataDir, ['fstar_figs_' datestr(now,'yyyymmdd_HHMMSS')]);

%% ---------------- Discover data files (XLSX; natural sort) --------------
files = dir(fullfile(dataDir,'*.xlsx'));
if isempty(files), error('No .xlsx files found in %s.', dataDir); end

names = {files.name};
maskSkip = startsWith(names,'sim_results_') | startsWith(names,'fit_params_') | ...
           startsWith(names,'deltaZ_results_') | startsWith(names,'fstar_separation_') | ...
           startsWith(names,'fstar_figs_');
files = files(~maskSkip);
if isempty(files), error('No .xlsx files left after skipping exports.'); end

keys  = natpad_keys({files.name});
[~,ord] = sort(keys); files = files(ord);

%% ---------------- Load & group by electrode and condition ----------------
H = struct();         % H.(electrode).(cond) -> cell array of recs
allRecs = {};

allowedConds = [];
if ~isempty(CONDS), allowedConds = lower(string(CONDS)); end

for k = 1:numel(files)
    fname = files(k).name; fpath = fullfile(files(k).folder, fname);

    tok = regexp(fname, '^(E\d+)_([A-Za-z0-9]+)_', 'tokens', 'once');
    if isempty(tok)
        % fallback: put into EALL and infer cond from first token
        [~,stem,~] = fileparts(fname);
        us = strfind(stem,'_');
        if isempty(us), electrode='EALL'; cond=lower(stem);
        else, electrode='EALL'; cond=lower(stem(1:us(1)-1));
        end
    else
        electrode = tok{1};
        cond      = lower(tok{2});
    end

    if ~isempty(allowedConds) && ~ismember(string(cond), allowedConds)
        fprintf('Skip (cond not in CONDS): %s (cond=%s)\n', fname, cond);
        continue;
    end

    T = readmatrix(fpath);
    if size(T,2) < 2
        fprintf('Skip (need Re/Im cols): %s\n', fname); continue;
    end
    Re = T(:,1); Im = T(:,2);
    ok = isfinite(Re) & isfinite(Im);
    Re = Re(ok); Im = Im(ok);
    if numel(Re) < 4, fprintf('Skip (too few points): %s\n', fname); continue; end

    % Frequency from col3 or Start/Stop
    f = [];
    if size(T,2) >= 3
        fc = T(ok,3);
        g2 = isfinite(fc) & (fc > 0);
        if nnz(g2) >= 4
            f  = fc(g2); Re = Re(g2); Im = Im(g2);
        end
    end
    if isempty(f)
        try
            [fStart, fStop] = read_start_stop_xlsx(fpath);
            f = logspace(log10(fStart), log10(fStop), numel(Re)).';
        catch ME
            fprintf('Skip (no usable freq; need col3 or Start/Stop): %s ; %s\n', fname, ME.message);
            continue;
        end
    end

    keep = isfinite(f) & (f > 0);
    if ~isempty(F_TRIM_MIN_HZ), keep = keep & (f >= F_TRIM_MIN_HZ); end
    if ~isempty(F_TRIM_MAX_HZ), keep = keep & (f <= F_TRIM_MAX_HZ); end
    f = f(keep); Re = Re(keep); Im = Im(keep);
    if numel(f) < 5
        fprintf('Skip (too few usable points after trim): %s\n', fname); continue;
    end

    rec.f    = f(:);
    rec.mag  = hypot(Re, Im);
    rec.lmz  = log10(max(rec.mag, realmin));
    rec.phu  = (180/pi) * unwrap(atan2(Im, Re));
    rec.name = fname;

    if ~isfield(H, electrode), H.(electrode) = struct(); end
    if ~isfield(H.(electrode), cond), H.(electrode).(cond) = {}; end
    H.(electrode).(cond){end+1} = rec; %#ok<AGROW>
    allRecs{end+1} = rec; %#ok<AGROW>
end

% Inventory
electrodes = fieldnames(H);
fprintf('Found electrodes: %s\n', strjoin(electrodes, ', '));
for e = 1:numel(electrodes)
    E = electrodes{e};
    fprintf('  %s -> conditions: %s\n', E, strjoin(fieldnames(H.(E)), ', '));
end
if isempty(electrodes), error('No usable datasets loaded.'); end

%% ---------------- Baseline detection (robust) ----------------------------
baseLower = lower(BASELINE);
condsPresent = {};
for e = 1:numel(electrodes)
    condsPresent = union(condsPresent, fieldnames(H.(electrodes{e})));
end
condsPresent = lower(condsPresent(:)).';

if ~ismember(baseLower, condsPresent)
    idxClean = find(contains(condsPresent, 'clean', 'IgnoreCase', true), 1, 'first');
    if ~isempty(idxClean)
        fprintf('Baseline "%s" not found; using fuzzy "%s".\n', BASELINE, condsPresent{idxClean});
        baseLower = condsPresent{idxClean};
    else
        [fallbackLabel, nOcc] = most_common_condition(H, electrodes);
        if ~isempty(fallbackLabel)
            fprintf('Baseline "%s" not found; using most common "%s" (n=%d).\n', BASELINE, fallbackLabel, nOcc);
            baseLower = fallbackLabel;
        else
            error('No usable baseline found (no conditions present).');
        end
    end
end

hasAnyBaseline = any(arrayfun(@(i) isfield(H.(electrodes{i}), baseLower) && ...
    ~isempty(H.(electrodes{i}).(baseLower)), 1:numel(electrodes)));
if ~hasAnyBaseline
    error('Baseline "%s" not present on any electrode.', baseLower);
end

%% ---------------- Per-electrode analysis & plotting ----------------------
% Also assemble a collector figure with all electrodes.
numE = numel(electrodes);
figGrid = figure('Color','w','Units','inches','Position',[0.5 0.5 pubWidthGrid max(rowHeight*numE, rowHeight)], ...
                 'PaperPositionMode','auto','Name','All Electrodes — S(f)');
tlo = tiledlayout(figGrid, numE, 2, 'TileSpacing','compact','Padding','compact');

% Colors reused for per-condition curves
colorBank = lines(16);

% To export CSVs per electrode
CSV_rows = {};

for e = 1:numE
    E = electrodes{e};
    condsE = lower(fieldnames(H.(E)));
    if isempty(condsE), fprintf('No conditions for %s; skip.\n', E); continue; end
    if ~ismember(baseLower, condsE) || isempty(H.(E).(baseLower))
        fprintf('Electrode %s lacks baseline "%s"; skip.\n', E, baseLower);
        continue;
    end

    % Overlap grid for this electrode
    recsE = {};
    for ci = 1:numel(condsE), recsE = [recsE, H.(E).(condsE{ci})]; end %#ok<AGROW>
    [fLoE, fHiE] = overlap_span(recsE);
    if ~(isfinite(fLoE) && isfinite(fHiE) && fHiE > fLoE)
        fprintf('No frequency overlap on %s; skip.\n', E); continue;
    end
    fgrid_e = logspace(log10(fLoE), log10(fHiE), NGRID);

    % Build matrices for this electrode
    Glmz_e = struct(); Gph_e = struct();
    for ci = 1:numel(condsE) 
        c = condsE{ci};
        recs_c = H.(E).(c);
        Mm = nan(numel(fgrid_e), numel(recs_c));
        Mp = nan(numel(fgrid_e), numel(recs_c));
        for r = 1:numel(recs_c)
            [fu, iu] = unique(recs_c{r}.f, 'stable');
            Mm(:,r) = interp1(fu, reallog10(recs_c{r}.mag(iu)), fgrid_e, 'pchip', 'extrap');
            Mp(:,r) = interp1(fu,        recs_c{r}.phu(iu),     fgrid_e, 'pchip', 'extrap');
        end
        Glmz_e.(c) = Mm;  Gph_e.(c) = Mp;
    end

    otherE = condsE(~strcmpi(condsE, baseLower));
    otherE = cellstr(otherE);

    % --- Pooled (within electrode; robust for n=1) ---
    [sep_mag_pooled_e,   fstar_mag_e]   = pick_fstar_pooled(Glmz_e, fgrid_e, baseLower, otherE);
    [sep_phase_pooled_e, fstar_phase_e] = pick_fstar_pooled(Gph_e,  fgrid_e, baseLower, otherE);

    % --- Per-condition separations (robust for n=1) ---
    [S_mag_per_e,  fstar_mag_per_e] = per_condition_separation(Glmz_e, fgrid_e, baseLower, otherE);
    [S_ph_per_e,   fstar_ph_per_e]  = per_condition_separation(Gph_e,  fgrid_e, baseLower, otherE);

    % ===================== Per-electrode figure (two axes) =====================
    figE = figure('Color','w','Units','inches','Position',[0.6 0.6 pubWidthOne pubHeightOne], ...
                  'PaperPositionMode','auto','Name',sprintf('%s — S(f) pooled + per-cond',E));
    tlE = tiledlayout(figE, 1, 2, 'TileSpacing','compact','Padding','compact');

    % ---- Axis 1: log10|Z| ----
    axM = nexttile(tlE, 1); hold(axM,'on'); grid(axM,'on'); box(axM,'on'); set(axM,'XScale','log');
    set(axM,'FontName',fontName,'FontSize',fontSizeAx);
    title(axM, sprintf('%s — S_{log_{10}|Z|}(f)', E), 'FontName',fontName,'FontSize',fontSizeTtl,'FontWeight','bold');
    xlabel(axM,'Frequency (Hz)','FontName',fontName,'FontSize',fontSizeLbl);
    ylabel(axM,'S_{log_{10}|Z|}(f)','FontName',fontName,'FontSize',fontSizeLbl,'Interpreter','tex');

    legEntriesM = {}; hLegM = [];

    % pooled (bold black)
    if any(isfinite(sep_mag_pooled_e))
        [hp,~] = plot_finite_segments(axM, fgrid_e, sep_mag_pooled_e, ...
                    '-', 'LineWidth',1.7, 'Color',[0 0 0], ...
                    'DisplayName','Pooled (others vs baseline)');
        if ~isempty(hp), hLegM(end+1) = hp(1); legEntriesM{end+1} = hp(1).DisplayName; end %#ok<AGROW>
        xline(axM, fstar_mag_e, '--k', sprintf('f* pooled = %.3g Hz', fstar_mag_e), ...
              'HandleVisibility','off','LabelVerticalAlignment','bottom','LabelOrientation','aligned');
        [sepMax, ~] = max(sep_mag_pooled_e(isfinite(sep_mag_pooled_e)));
        if isfinite(sepMax)
            ii = find(isfinite(sep_mag_pooled_e)); imax = ii(argmax(sep_mag_pooled_e(ii)));
            plot(axM, fgrid_e(imax), sepMax, 'ko', 'MarkerFaceColor','k','MarkerSize',4.5, 'HandleVisibility','off');
        end
    end

    % each condition (colored)
    for cidx = 1:numel(otherE)
        c = otherE{cidx};
        if ~isfield(S_mag_per_e,c), continue; end
        S = S_mag_per_e.(c);
        if ~any(isfinite(S)), continue; end
        clr = colorBank( 1+mod(cidx-1,size(colorBank,1)) , : );
        [hSeg,~] = plot_finite_segments(axM, fgrid_e, S, '-', 'LineWidth',1.2, 'Color', clr, ...
                                        'DisplayName', sprintf('%s vs baseline', c));
        if ~isempty(hSeg)
            hLegM(end+1) = hSeg(1); legEntriesM{end+1} = hSeg(1).DisplayName; %#ok<AGROW>
            [Smax, idx] = max(S(isfinite(S)));
            ii = find(isfinite(S));
            plot(axM, fgrid_e(ii(idx)), Smax, 'o', 'MarkerSize',4.0, ...
                 'MarkerFaceColor',clr,'MarkerEdgeColor','w','LineWidth',0.7,'HandleVisibility','off');
        end
    end

    % y-limits from stack
    SstackAllM = sep_mag_pooled_e(:);
    for cidx = 1:numel(otherE)
        c = otherE{cidx};
        if isfield(S_mag_per_e,c), SstackAllM = [SstackAllM; S_mag_per_e.(c)(:)]; end %#ok<AGROW>
    end
    if any(isfinite(SstackAllM)), ylim(axM, padlims(SstackAllM, 0.08, false)); end
    if ~isempty(hLegM), legend(axM, hLegM, legEntriesM, 'Location','best','Interpreter','tex','AutoUpdate','off'); end

    % ---- Axis 2: Phase ----
    axP = nexttile(tlE, 2); hold(axP,'on'); grid(axP,'on'); box(axP,'on'); set(axP,'XScale','log');
    set(axP,'FontName',fontName,'FontSize',fontSizeAx);
    title(axP, sprintf('%s — S_{phase}(f)', E), 'FontName',fontName,'FontSize',fontSizeTtl,'FontWeight','bold');
    xlabel(axP,'Frequency (Hz)','FontName',fontName,'FontSize',fontSizeLbl);
    ylabel(axP,'S_{phase}(f)','FontName',fontName,'FontSize',fontSizeLbl,'Interpreter','tex');

    legEntriesP = {}; hLegP = [];

    if any(isfinite(sep_phase_pooled_e))
        [hp,~] = plot_finite_segments(axP, fgrid_e, sep_phase_pooled_e, ...
                    '-', 'LineWidth',1.7, 'Color',[0 0 0], ...
                    'DisplayName','Pooled (others vs baseline)');
        if ~isempty(hp), hLegP(end+1) = hp(1); legEntriesP{end+1} = hp(1).DisplayName; end %#ok<AGROW>
        xline(axP, fstar_phase_e, '--k', sprintf('f* pooled = %.3g Hz', fstar_phase_e), ...
              'HandleVisibility','off','LabelVerticalAlignment','bottom','LabelOrientation','aligned');
        [sepMax, ~] = max(sep_phase_pooled_e(isfinite(sep_phase_pooled_e)));
        if isfinite(sepMax)
            ii = find(isfinite(sep_phase_pooled_e)); imax = ii(argmax(sep_phase_pooled_e(ii)));
            plot(axP, fgrid_e(imax), sepMax, 'ko', 'MarkerFaceColor','k','MarkerSize',4.5, 'HandleVisibility','off');
        end
    end

    for cidx = 1:numel(otherE)
        c = otherE{cidx};
        if ~isfield(S_ph_per_e,c), continue; end
        S = S_ph_per_e.(c); if ~any(isfinite(S)), continue; end
        clr = colorBank( 1+mod(cidx-1,size(colorBank,1)) , : );
        [hSeg,~] = plot_finite_segments(axP, fgrid_e, S, '-', 'LineWidth',1.2, 'Color', clr, ...
                                        'DisplayName', sprintf('%s vs baseline', c));
        if ~isempty(hSeg)
            hLegP(end+1) = hSeg(1); legEntriesP{end+1} = hSeg(1).DisplayName; %#ok<AGROW>
            [Smax, idx] = max(S(isfinite(S)));
            ii = find(isfinite(S));
            plot(axP, fgrid_e(ii(idx)), Smax, 'o', 'MarkerSize',4.0, ...
                 'MarkerFaceColor',clr,'MarkerEdgeColor','w','LineWidth',0.7,'HandleVisibility','off');
        end
    end

    SstackAllP = sep_phase_pooled_e(:);
    for cidx = 1:numel(otherE)
        c = otherE{cidx};
        if isfield(S_ph_per_e,c), SstackAllP = [SstackAllP; S_ph_per_e.(c)(:)]; end %#ok<AGROW>
    end
    if any(isfinite(SstackAllP)), ylim(axP, padlims(SstackAllP, 0.08, false)); end
    if ~isempty(hLegP), legend(axP, hLegP, legEntriesP, 'Location','best','Interpreter','tex','AutoUpdate','off'); end

    % Export per-electrode
    try
        fnE = sprintf('%s_%s_pooled_plus_percond', figOutPrefix, E);
        print(figE,[fnE '.png'],'-dpng',sprintf('-r%d',exportDPI));
        if doPDFexport, print(figE,[fnE '.pdf'],'-dpdf','-painters'); end
    catch ME
        warning('[%s] Export failed: %s', E, ME.message);
    end

% -------- Left cell: log10|Z| (WITH LEGEND + per-condition f* xlines) --------
axM2 = nexttile(tlo, (e-1)*2 + 1); 
hold(axM2,'on'); grid(axM2,'on'); box(axM2,'on'); set(axM2,'XScale','log');
set(axM2,'FontName',fontName,'FontSize',fontSizeAx);
if e==1, title(axM2,'S_{log_{10}|Z|}(f)','Interpreter','tex'); end
ylabel(axM2, E, 'FontWeight','bold','Rotation',0,'HorizontalAlignment','right','FontName',fontName);

hLegM2 = gobjects(0,1); legStrM2 = {};

% pooled (black) — visible in legend; pooled f* guide hidden
if any(isfinite(sep_mag_pooled_e))
    [hp,~] = plot_finite_segments(axM2, fgrid_e, sep_mag_pooled_e, '-', ...
        'LineWidth',1.5, 'Color',[0 0 0], 'DisplayName','Pooled (others vs baseline)');
    if ~isempty(hp), hLegM2(end+1) = hp(1); legStrM2{end+1} = hp(1).DisplayName; end %#ok<AGROW>
    if isfinite(fstar_mag_e), xline(axM2, fstar_mag_e, '--k', 'HandleVisibility','off'); end
end

% each condition (colored) — visible in legend; per-condition f*_c guides hidden
for cidx = 1:numel(otherE)
    c = otherE{cidx}; 
    if ~isfield(S_mag_per_e,c), continue; end
    S = S_mag_per_e.(c); 
    if ~any(isfinite(S)), continue; end

    clr = colorBank(1+mod(cidx-1,size(colorBank,1)), :);
    [hSeg,~] = plot_finite_segments(axM2, fgrid_e, S, '-', ...
        'LineWidth',1.0, 'Color', clr, 'DisplayName', sprintf('%s vs baseline', c));
    if ~isempty(hSeg), hLegM2(end+1) = hSeg(1); legStrM2{end+1} = hSeg(1).DisplayName; end %#ok<AGROW>

    if isfield(fstar_mag_per_e,c)
        fstar_c = fstar_mag_per_e.(c);
        if isfinite(fstar_c)
            xline(axM2, fstar_c, '--', 'Color', clr, 'HandleVisibility','off');
        end
    end
end

if any(isfinite(SstackAllM)), ylim(axM2, padlims(SstackAllM, 0.10, false)); end
if ~isempty(hLegM2)
    legend(axM2, hLegM2, legStrM2, 'Location','best','Interpreter','tex','AutoUpdate','off');
end

% -------- Right cell: phase (WITH LEGEND + per-condition f* xlines) --------
axP2 = nexttile(tlo, (e-1)*2 + 2); 
hold(axP2,'on'); grid(axP2,'on'); box(axP2,'on'); set(axP2,'XScale','log');
set(axP2,'FontName',fontName,'FontSize',fontSizeAx);
if e==1, title(axP2,'S_{phase}(f)','Interpreter','tex'); end

hLegP2 = gobjects(0,1); legStrP2 = {};

% pooled (black) — visible in legend; pooled f* guide hidden
if any(isfinite(sep_phase_pooled_e))
    [hp,~] = plot_finite_segments(axP2, fgrid_e, sep_phase_pooled_e, '-', ...
        'LineWidth',1.5, 'Color',[0 0 0], 'DisplayName','Pooled (others vs baseline)');
    if ~isempty(hp), hLegP2(end+1) = hp(1); legStrP2{end+1} = hp(1).DisplayName; end %#ok<AGROW>
    if isfinite(fstar_phase_e), xline(axP2, fstar_phase_e, '--k', 'HandleVisibility','off'); end
end

% each condition (colored) — visible in legend; per-condition f*_c guides hidden
for cidx = 1:numel(otherE)
    c = otherE{cidx}; 
    if ~isfield(S_ph_per_e,c), continue; end
    S = S_ph_per_e.(c); 
    if ~any(isfinite(S)), continue; end

    clr = colorBank(1+mod(cidx-1,size(colorBank,1)), :);
    [hSeg,~] = plot_finite_segments(axP2, fgrid_e, S, '-', ...
        'LineWidth',1.0, 'Color', clr, 'DisplayName', sprintf('%s vs baseline', c));
    if ~isempty(hSeg), hLegP2(end+1) = hSeg(1); legStrP2{end+1} = hSeg(1).DisplayName; end %#ok<AGROW>

    if isfield(fstar_ph_per_e,c)
        fstar_c = fstar_ph_per_e.(c);
        if isfinite(fstar_c)
            xline(axP2, fstar_c, '--', 'Color', clr, 'HandleVisibility','off');
        end
    end
end

if any(isfinite(SstackAllP)), ylim(axP2, padlims(SstackAllP, 0.10, false)); end
if ~isempty(hLegP2)
    legend(axP2, hLegP2, legStrP2, 'Location','best','Interpreter','tex','AutoUpdate','off');
end



    % ------- accumulate CSV rows (per electrode) -------
    CSV_rows{end+1,1} = E; %#ok<AGROW>
    CSV_rows{end,2}   = fgrid_e(:);
    CSV_rows{end,3}   = sep_mag_pooled_e(:);
    CSV_rows{end,4}   = sep_phase_pooled_e(:);
    CSV_rows{end,5}   = S_mag_per_e;
    CSV_rows{end,6}   = S_ph_per_e;
end

%% ---------------- CSV exports (pooled + per-cond, all E) ----------------
% Pooled (stack across electrodes with NaN-pad to longest)
try
    % Stack pooled by electrode (variable-length)
    maxN = 0; for i=1:size(CSV_rows,1), maxN = max(maxN, numel(CSV_rows{i,2})); end
    freq_col = nan(maxN,1);
    Smag_cols = nan(maxN, size(CSV_rows,1));
    Sph_cols  = nan(maxN, size(CSV_rows,1));
    elabels   = strings(1, size(CSV_rows,1));
    for i=1:size(CSV_rows,1)
        fi = CSV_rows{i,2}; n = numel(fi);
        freq_col(1:n) = fi;
        Smag_cols(1:n,i) = CSV_rows{i,3};
        Sph_cols(1:n,i)  = CSV_rows{i,4};
        elabels(i) = string(CSV_rows{i,1});
    end
    T_pooled = table(freq_col, Smag_cols, Sph_cols, 'VariableNames', ...
        {'freq_Hz','S_log10AbsZ_pooled_byE','S_phase_pooled_byE'});
    writetable(T_pooled, fullfile(dataDir, ['fstar_separation_pooled_byE_' datestr(now,'yyyymmdd_HHMMSS') '.csv']));
catch ME
    warning('Export pooled-by-electrode CSV failed: %s', ME.message);
end

% Per-condition wide export (stacking across electrodes inside each cond struct is complex;
% here we export per-electrode per-condition bundles in a simple cell-based .mat for fidelity)
try
    MAT = struct('perElectrode', cell(1, size(CSV_rows,1)));
    for i=1:size(CSV_rows,1)
        MAT.perElectrode(i).electrode = CSV_rows{i,1};
        MAT.perElectrode(i).freq_Hz   = CSV_rows{i,2};
        MAT.perElectrode(i).S_mag     = CSV_rows{i,5};
        MAT.perElectrode(i).S_phase   = CSV_rows{i,6};
    end
    save(fullfile(dataDir, ['fstar_per_condition_byE_' datestr(now,'yyyymmdd_HHMMSS') '.mat']), '-struct','MAT');
catch ME
    warning('Export per-condition MAT failed: %s', ME.message);
end

%% ---------------- Helper functions --------------------------------------
function [fStart, fStop] = read_start_stop_xlsx(fpath)
    raw = readcell(fpath);
    [r_sf,c_sf] = find(strcmp(raw, 'Start Frequency'), 1);
    [r_sp,c_sp] = find(strcmp(raw, 'Stop Frequency'), 1);
    assert(~isempty(r_sf)&&~isempty(r_sp), 'Start/Stop markers not found.');
    fStart = raw{r_sf,c_sf+1};
    fStop  = raw{r_sp,c_sp+1};
    validateattributes(fStart,{'numeric'},{'scalar','positive','finite'});
    validateattributes(fStop, {'numeric'},{'scalar','positive','finite'});
    assert(fStop > fStart, 'Stop Frequency must be > Start Frequency.');
end

function [fLo,fHi] = overlap_span(recs)
    fLo = -inf; fHi = +inf;
    for r=1:numel(recs)
        fi = recs{r}.f; if isempty(fi)||~isvector(fi), continue; end
        fLo = max(fLo, min(fi)); fHi = min(fHi, max(fi));
    end
end

function y = reallog10(x), y = log10(max(x, realmin)); end
function i = argmax(v), [~,i] = max(v); end

% ---------- Robust pooled separation (works for n=1) ----------
function [sep,fstar] = pick_fstar_pooled(M, fgrid, baseline, otherConds)
    sep = nan(numel(fgrid),1); 
    fstar = fgrid(1);
    if ~isfield(M, baseline) || isempty(M.(baseline)), return; end
    xb = M.(baseline);

    % Stack "other" conditions columns
    Moth = [];
    for i=1:numel(otherConds)
        c = otherConds{i};
        if isfield(M, c) && ~isempty(M.(c)), Moth = [Moth, M.(c)]; end %#ok<AGROW>
    end
    if isempty(Moth), return; end

    % Global sigma at each frequency (across all conditions)
    sigma_glob = global_sigma_from_struct(M);

    for i=1:numel(fgrid)
        vb = xb(i,:); vb = vb(isfinite(vb)); nb = numel(vb);
        vo = Moth(i,:); vo = vo(isfinite(vo)); no = numel(vo);
        if nb==0 || no==0, continue; end

        mb = mean(vb,'omitnan'); mo = mean(vo,'omitnan');
        sb = std(vb,0,'omitnan'); so = std(vo,0,'omitnan');

        if ~(isfinite(sb)&&sb>0) || ~(isfinite(so)&&so>0)
            % Fallback stds within groups if possible
            sRob = robust_sigma([vb(:); vo(:)]);
            if ~(isfinite(sb)&&sb>0), sb = sRob; end
            if ~(isfinite(so)&&so>0), so = sRob; end
        end

        df = nb+no-2;
        if df > 0
            sp = sqrt(((nb-1)*sb^2 + (no-1)*so^2) / df);
        else
            % KEY: Use global per-frequency sigma, not pairwise-based
            sp = sigma_glob(i);
        end
        sp = max(sp, eps);
        sep(i) = abs(mb - mo) / sp;
    end

    [~,idx] = max(sep); 
    if isempty(idx) || ~isfinite(idx), idx = 1; end
    fstar = fgrid(idx);
end


% ---------- Robust per-condition separation (works for n=1) ----------
function [S_per, fstar_per] = per_condition_separation(M, fgrid, baseline, otherConds)
    S_per = struct(); fstar_per = struct();
    if ~isfield(M, baseline) || isempty(M.(baseline)), return; end
    Mb = M.(baseline);

    % Global per-frequency sigma across ALL conditions (once)
    sigma_glob = global_sigma_from_struct(M);

    if isstring(otherConds), otherConds = cellstr(otherConds); end
    for i=1:numel(otherConds)
        c = otherConds{i};
        if ~isfield(M, c) || isempty(M.(c)), continue; end
        Mc = M.(c);

        S = nan(numel(fgrid),1);
        for k=1:numel(fgrid)
            vb = Mb(k,:); vb = vb(isfinite(vb)); nb = numel(vb);
            vc = Mc(k,:); vc = vc(isfinite(vc)); nc = numel(vc);
            if nb==0 || nc==0, continue; end

            mb = mean(vb,'omitnan'); mc = mean(vc,'omitnan');
            sb = std(vb,0,'omitnan'); sc = std(vc,0,'omitnan');

            if ~(isfinite(sb)&&sb>0) || ~(isfinite(sc)&&sc>0)
                sRob = robust_sigma([vb(:); vc(:)]);
                if ~(isfinite(sb)&&sb>0), sb = sRob; end
                if ~(isfinite(sc)&&sc>0), sc = sRob; end
            end

            df = nb+nc-2;
            if df > 0
                sp = sqrt(((nb-1)*sb^2 + (nc-1)*sc^2) / df);
            else
                % KEY: Use global per-frequency sigma here
                sp = sigma_glob(k);
            end
            sp = max(sp, eps);
            S(k) = abs(mb - mc) / sp;
        end
        S_per.(c) = S;

        ii = find(isfinite(S));
        if isempty(ii), fstar_per.(c) = fgrid(1); 
        else, [~,m] = max(S(ii)); fstar_per.(c) = fgrid(ii(m)); 
        end
    end
end


function s = robust_sigma(x)
    x = x(isfinite(x));
    if isempty(x), s = 1; return; end
    madv = mad(x,1); s = 1.4826*madv;
    if ~(isfinite(s) && s>0)
        q1 = quantile(x,0.25); q3 = quantile(x,0.75);
        s = (q3-q1)/1.349;
    end
    if ~(isfinite(s) && s>0)
        s = max(eps, 0.01*max(abs(x)));
    end
end

function L = padlims(y, frac, symmetricToZero)
    if nargin<2||isempty(frac), frac=0.08; end
    if nargin<3, symmetricToZero=false; end
    y = y(:); y = y(isfinite(y)); if isempty(y), L=[-1 1]; return; end
    lo=min(y); hi=max(y); if lo==hi, lo=lo-1; hi=hi+1; end
    dy=hi-lo; lo=lo-frac*dy; hi=hi+frac*dy;
    if symmetricToZero, lo=min(lo,0); hi=max(hi,0); end
    L=[lo hi];
end

function keys = natpad_keys(names)
    padWidth = 12; if isstring(names), names = cellstr(names); end
    if ischar(names), names = cellstr(names); end
    keys = cell(size(names));
    for i = 1:numel(names)
        s = names{i}; if ~ischar(s), s = char(s); end
        toks = regexp(s, '\d+|\D+', 'match'); out = '';
        for t = 1:numel(toks)
            tok = toks{t};
            if ~isempty(tok) && all(isstrprop(tok,'digit'))
                n = str2double(tok); if ~isfinite(n), n = 0; end
                out = [out, sprintf(['%0',num2str(padWidth),'d'], n)]; %#ok<AGROW>
            else
                out = [out, lower(tok)]; %#ok<AGROW>
            end
        end
        keys{i} = out;
    end
end

function [label, countMax] = most_common_condition(H, electrodes)
    counts = containers.Map('KeyType','char','ValueType','double');
    for e = 1:numel(electrodes)
        E = electrodes{e};
        condsE = fieldnames(H.(E));
        for i = 1:numel(condsE)
            c = lower(condsE{i});
            if ~isKey(counts, c), counts(c) = 0; end
            counts(c) = counts(c) + numel(H.(E).(c));
        end
    end
    if counts.Count == 0, label = ''; countMax = 0; return; end
    ks = counts.keys; vs = counts.values;
    [countMax, idx] = max(cell2mat(vs)); label = ks{idx};
end

% --------- plot only finite contiguous segments (prevents broken lines) ---------
function [hSeg, segIdx] = plot_finite_segments(ax, x, y, varargin)
    x = x(:); y = y(:);
    fin = isfinite(x) & isfinite(y);
    hSeg = gobjects(0); segIdx = {};
    if ~any(fin), return; end
    % find contiguous runs of 'true' in fin
    d = diff([false; fin; false]);
    starts = find(d==1);
    ends   = find(d==-1) - 1;
    for i = 1:numel(starts)
        ii = starts(i):ends(i);
        if numel(ii) < 2, continue; end
        if isempty(hSeg)
            % first segment: keep legend/display name
            h = plot(ax, x(ii), y(ii), varargin{:});
        else
            % subsequent segments: same style but hide from legend
            args = varargin;
            for k = 1:2:numel(args)
                if strcmpi(args{k}, 'DisplayName')
                    args{k+1} = ''; % remove name for extra segments
                end
            end
            h = plot(ax, x(ii), y(ii), args{:}, 'HandleVisibility','off');
        end
        hSeg(end+1) = h; %#ok<AGROW>
        segIdx{end+1} = ii; %#ok<AGROW>
    end
end
function sigma_glob = global_sigma_from_struct(M)
% sigma_glob(k): robust sigma across ALL conditions & replicates at each frequency k
    % Collect all matrices column-wise
    fn = fieldnames(M);
    if isempty(fn)
        sigma_glob = [];
        return;
    end
    % Determine number of frequency points from first non-empty field
    nF = [];
    for i = 1:numel(fn)
        if isfield(M, fn{i}) && ~isempty(M.(fn{i}))
            nF = size(M.(fn{i}),1); break;
        end
    end
    if isempty(nF), sigma_glob = []; return; end

    sigma_glob = nan(nF,1);
    for k = 1:nF
        pool_k = [];
        for i = 1:numel(fn)
            if ~isempty(M.(fn{i}))
                row = M.(fn{i})(k,:); row = row(isfinite(row));
                pool_k = [pool_k; row(:)]; %#ok<AGROW>
            end
        end
        if isempty(pool_k)
            sigma_glob(k) = NaN;
        else
            sigma_glob(k) = robust_sigma(pool_k);
            if ~(isfinite(sigma_glob(k)) && sigma_glob(k) > 0)
                sigma_glob(k) = max(eps, 0.01*max(abs(pool_k)));
            end
        end
    end
end
