function [fstar_mag, fstar_phase, Tout] = single_freq_delta_Z(dataDir, varargin)
% SINGLE_FREQ_DELTA_Z  (Electrode-specific deltas; stats on log10|Z|)
%   Computes and plots single-frequency summaries for EIS datasets with
%   electrode-specific baselines (same electrode 'cleaned').
%   f* (magnitude) is selected by maximizing separation on log10|Z|.
%
% Usage:
%   [fstar_mag, fstar_phase, T] = single_freq_delta_Z
%   [fstar_mag, fstar_phase, T] = single_freq_delta_Z('C:\path', 'FIXED_FREQ_HZ', 30000)
%
% Expected filename pattern:
%   E<digit+>_<condition>_*.xlsx  (conds default: {'cleaned','e4','e6','e9'})
%
% Columns:
%   col1: Re(Z) [Ohm]
%   col2: Im(Z) [Ohm]
%   col3 (optional): Frequency [Hz]
%
% Output table main columns:
%   file, electrode, condition,
%   fstar_mag_Hz,  AbsZ_fstar_mag,  Delta_AbsZ_vs_cleaned_sameE, Baseline_AbsZ_sameE,
%                  log10AbsZ_fstar_mag, Delta_log10AbsZ_vs_cleaned_sameE, Baseline_log10AbsZ_sameE,
%   fstar_phase_Hz, phase_deg_fstar_phase, Delta_phase_deg_vs_cleaned_sameE, Baseline_phase_deg_sameE
%
% Robustness (optional):
%   Delta_AbsZ_fixedHz_vs_cleaned_sameE, Delta_log10AbsZ_fixedHz_vs_cleaned_sameE
%   Delta_AbsZ_bandMed_vs_cleaned_sameE, Delta_log10AbsZ_bandMed_vs_cleaned_sameE
%
% Notes for Methods:
%   "Magnitude f* was selected by maximizing the standardized separation on log10|Z|.
%    Statistical analyses (ANOVA/pairs) were performed on log10|Z|; results are
%    reported as Δ|Z| (Ω) for interpretability."

    if nargin < 1 || isempty(dataDir), dataDir = pwd; end

    % ---------------- User settings ----------------
    p = inputParser;
    addParameter(p,'CONDS',{'cleaned','e4','e6','e9'});
    addParameter(p,'BASELINE','cleaned');
    addParameter(p,'F_TRIM_MIN_HZ',1e4,@(x)isnumeric(x)&&isscalar(x));
    addParameter(p,'NGRID',300,@(x)isnumeric(x)&&isscalar(x));
    addParameter(p,'FIXED_FREQ_HZ',[],@(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0));
    addParameter(p,'BAND',[1e4 1e5],@(x)isnumeric(x)&&numel(x)==2&&x(1)>0&&x(2)>x(1));
    parse(p,varargin{:});
    condsWanted   = string(p.Results.CONDS);
    baselineCond  = string(p.Results.BASELINE);
    F_TRIM_MIN_HZ = p.Results.F_TRIM_MIN_HZ;
    nGridPick     = p.Results.NGRID;
    FIXED_FREQ_HZ = p.Results.FIXED_FREQ_HZ;
    BAND          = p.Results.BAND;

    % ------- Figure/Export settings -------
    clr = [  0 114 178;   213 94  0;   0 158 115;   204 121 167;  230 159 0]/255;
    fontName     = 'Helvetica'; fontSizeAx=9; fontSizeLbl=10; fontSizeTtl=10;
    dotSize      = 28; boxWidth=0.32; exportDPI=300; doPDFexport=true;
    pubWidth=3.6; pubHeight=3.0;
    figOutPrefix = fullfile(dataDir, ['deltaZ_figs_' datestr(now,'yyyymmdd_HHMMSS')]);

    % ---------------- Discover files ----------------
    files = dir(fullfile(dataDir, '*.xlsx'));
    maskSkip = startsWith({files.name}, 'sim_results_') | startsWith({files.name}, 'fit_params_') ...
             | startsWith({files.name}, 'deltaZ_results_');
    files = files(~maskSkip);
    assert(~isempty(files), 'No .xlsx files found (after skipping exports).');

    % ---------------- Load & group: by electrode and condition ----------------
    H = struct();
    allRecs = {};
    recMeta = {};  % {fileName, electrode, condition}

    for k = 1:numel(files)
        fname = files(k).name; fpath = fullfile(files(k).folder, fname);
        tok = regexp(fname, '^(E\d+)_([A-Za-z0-9]+)_', 'tokens', 'once');
        if isempty(tok), fprintf('Skip (name pattern): %s\n', fname); continue; end
        electrode = tok{1}; cond = lower(tok{2});
        if ~ismember(string(cond), condsWanted), continue; end

        M = readmatrix(fpath);
        if size(M,2) < 2, fprintf('Skip (need Re/Im): %s\n', fname); continue; end
        Re = M(:,1); Im = M(:,2);
        ok = isfinite(Re) & isfinite(Im); Re=Re(ok); Im=Im(ok);

        % Frequency parsing
        if size(M,2) >= 3
            fc = M(ok,3); g2 = isfinite(fc) & (fc > 0);
            if nnz(g2) >= 4
                f  = fc(g2); Re=Re(g2); Im=Im(g2);
            else
                [fStart, fStop] = read_start_stop(fpath);
                f = logspace(log10(fStart), log10(fStop), numel(Re)).';
            end
        else
            [fStart, fStop] = read_start_stop(fpath);
            f = logspace(log10(fStart), log10(fStop), numel(Re)).';
        end

        keep = isfinite(f) & (f > 0);
        if ~isempty(F_TRIM_MIN_HZ), keep = keep & (f >= F_TRIM_MIN_HZ); end
        f = f(keep); Re = Re(keep); Im = Im(keep);
        if numel(f) < 5, fprintf('Too few usable points after trim: %s\n', fname); continue; end

        rec.f    = f(:);
        rec.mag  = hypot(Re, Im);                       % |Z|
        rec.lmz  = log10(max(rec.mag, realmin));        % log10|Z|
        rec.phu  = rad2deg(unwrap(atan2(Im, Re)));      % unwrapped phase (deg)
        rec.name = fname;

        if ~isfield(H, electrode), H.(electrode) = struct(); end
        if ~isfield(H.(electrode), cond), H.(electrode).(cond) = {}; end
        H.(electrode).(cond){end+1} = rec; %#ok<AGROW>

        allRecs{end+1} = rec; %#ok<AGROW>
        recMeta(end+1,:) = {fname, electrode, lower(cond)}; %#ok<AGROW>
    end

    % Ensure baseline exists somewhere
    electrodes = fieldnames(H);
    hasAnyBaseline = any(arrayfun(@(i) isfield(H.(electrodes{i}), baselineCond) && ~isempty(H.(electrodes{i}).(baselineCond)), 1:numel(electrodes)));
    assert(hasAnyBaseline, 'No electrode has a "%s" baseline.', baselineCond);

    % ---------------- Global overlap for f* selection ----------------
    [fLoAll, fHiAll] = overlap_span(allRecs);
    assert(isfinite(fLoAll) && isfinite(fHiAll) && fHiAll > fLoAll, 'No overlapping frequency span across all records.');
    fgrid = logspace(log10(fLoAll), log10(fHiAll), nGridPick);

    % Build condition-wise matrices on fgrid for f* selection (log10|Z|, phase)
    Glmz = struct(); Gph = struct();
    for ci = 1:numel(condsWanted)
        cond = char(condsWanted(ci));
        recs_c = {};
        for e = 1:numel(electrodes)
            E = electrodes{e};
            if isfield(H.(E), cond), recs_c = [recs_c, H.(E).(cond)]; end %#ok<AGROW>
        end
        if isempty(recs_c), Glmz.(cond)=nan(numel(fgrid),0); Gph.(cond)=nan(numel(fgrid),0); continue; end
        Mm  = nan(numel(fgrid), numel(recs_c));
        Mp  = nan(numel(fgrid), numel(recs_c));
        for r = 1:numel(recs_c)
            [fu, iu] = unique(recs_c{r}.f, 'stable');
            Mm(:,r) = interp1(fu, reallog10(recs_c{r}.mag(iu)), fgrid, 'pchip', 'extrap');  % log10|Z|
            Mp(:,r) = interp1(fu, recs_c{r}.phu(iu),            fgrid, 'pchip', 'extrap');  % phase
        end
        Glmz.(cond) = Mm; Gph.(cond) = Mp;
    end

    % Choose f* (same for all electrodes)
    [~, fstar_mag]   = pick_fstar(Glmz, fgrid, char(baselineCond), cellstr(condsWanted)); % sep on log10|Z|
    [~, fstar_phase] = pick_fstar(Gph,  fgrid, char(baselineCond), cellstr(condsWanted)); % sep on phase

    fprintf('Chosen f*_mag   = %.4g Hz (global separation on log10|Z|)\n', fstar_mag);
    fprintf('Chosen f*_phase = %.4g Hz (global separation on phase)\n',    fstar_phase);

    % ---------------- Precompute baselines @ f* ----------------
    cleanedMean_mag  = containers.Map; % |Z| @ f*mag
    cleanedMean_lmz  = containers.Map; % log10|Z| @ f*mag
    cleanedMean_ph   = containers.Map; % phase @ f*phase
    if ~isempty(FIXED_FREQ_HZ)
        cleanedMean_mag_fix = containers.Map;
        cleanedMean_lmz_fix = containers.Map;
    end

    for e = 1:numel(electrodes)
        E = electrodes{e};
        if isfield(H.(E), char(baselineCond)) && ~isempty(H.(E).(char(baselineCond)))
            recsB = H.(E).(char(baselineCond));
            vB_mag = interp_per_files(recsB, fstar_mag,   'mag');
            vB_lmz = interp_per_files(recsB, fstar_mag,   'lmz');
            vB_ph  = interp_per_files(recsB, fstar_phase, 'phase');
            cleanedMean_mag(E) = mean(vB_mag, 'omitnan');
            cleanedMean_lmz(E) = mean(vB_lmz, 'omitnan');
            cleanedMean_ph(E)  = mean(vB_ph,  'omitnan');
            if ~isempty(FIXED_FREQ_HZ)
                vB_magF = interp_per_files(recsB, FIXED_FREQ_HZ, 'mag');
                vB_lmzF = interp_per_files(recsB, FIXED_FREQ_HZ, 'lmz');
                cleanedMean_mag_fix(E) = mean(vB_magF,'omitnan');
                cleanedMean_lmz_fix(E) = mean(vB_lmzF,'omitnan');
            end
        end
    end

    % ---------------- Per-file values & deltas ----------------
    nFiles = numel(allRecs);
    absZ_star = nan(nFiles,1); dAbsZ = nan(nFiles,1);
    lmz_star  = nan(nFiles,1); dlmz  = nan(nFiles,1);
    ph_star   = nan(nFiles,1); dph   = nan(nFiles,1);
    baseUsed_mag = nan(nFiles,1); baseUsed_lmz = nan(nFiles,1); baseUsed_ph = nan(nFiles,1);

    % robustness holders
    dAbsZ_fix = nan(nFiles,1); dlmz_fix = nan(nFiles,1);
    dAbsZ_band = nan(nFiles,1); dlmz_band = nan(nFiles,1);

    for i = 1:nFiles
        rec = allRecs{i}; electrode = recMeta{i,2};

        % at f*_mag
        absZ_star(i) = interp_one(rec, fstar_mag, 'mag');
        lmz_star(i)  = interp_one(rec, fstar_mag, 'lmz');
        if isKey(cleanedMean_mag, electrode)
            muZ  = cleanedMean_mag(electrode);
            muL  = cleanedMean_lmz(electrode);
            dAbsZ(i) = absZ_star(i) - muZ;
            dlmz(i)  = lmz_star(i)  - muL;
            baseUsed_mag(i) = muZ; baseUsed_lmz(i) = muL;
        end

        % phase @ f*_phase
        ph_star(i) = interp_one(rec, fstar_phase, 'phase');
        if isKey(cleanedMean_ph, electrode)
            muP = cleanedMean_ph(electrode);
            dph(i) = ph_star(i) - muP;
            baseUsed_ph(i) = muP;
        end

        % fixed frequency robustness
        if ~isempty(FIXED_FREQ_HZ) && isKey(cleanedMean_mag_fix, electrode)
            zF  = interp_one(rec, FIXED_FREQ_HZ, 'mag');
            lF  = interp_one(rec, FIXED_FREQ_HZ, 'lmz');
            dAbsZ_fix(i) = zF - cleanedMean_mag_fix(electrode);
            dlmz_fix(i)  = lF - cleanedMean_lmz_fix(electrode);
        end

        % band-median robustness (median over BAND on each record)
        inBand = rec.f >= BAND(1) & rec.f <= BAND(2);
        if any(inBand)
            zMed  = median(rec.mag(inBand),'omitnan');
            lMed  = median(reallog10(rec.mag(inBand)),'omitnan');
            % need baseline band medians for the same electrode:
            if isfield(H.(electrode), char(baselineCond))
                zB=[]; lB=[];
                for rr = 1:numel(H.(electrode).(char(baselineCond)))
                    rB = H.(electrode).(char(baselineCond)){rr};
                    vb = rB.f >= BAND(1) & rB.f <= BAND(2);
                    if any(vb)
                        zB(end+1) = median(rB.mag(vb),'omitnan'); %#ok<AGROW>
                        lB(end+1) = median(reallog10(rB.mag(vb)),'omitnan'); %#ok<AGROW>
                    end
                end
                if ~isempty(zB)
                    dAbsZ_band(i) = zMed - mean(zB,'omitnan');
                    dlmz_band(i)  = lMed - mean(lB,'omitnan');
                end
            end
        end
    end

    Tout = table( ...
        string(recMeta(:,1)), string(recMeta(:,2)), string(recMeta(:,3)), ...
        repmat(fstar_mag,   nFiles,1), absZ_star, dAbsZ, baseUsed_mag, lmz_star, dlmz, baseUsed_lmz, ...
        repmat(fstar_phase, nFiles,1), ph_star,   dph,    baseUsed_ph, ...
        'VariableNames', {'file','electrode','condition', ...
        'fstar_mag_Hz','AbsZ_fstar_mag','Delta_AbsZ_vs_cleaned_sameE','Baseline_AbsZ_sameE', ...
        'log10AbsZ_fstar_mag','Delta_log10AbsZ_vs_cleaned_sameE','Baseline_log10AbsZ_sameE', ...
        'fstar_phase_Hz','phase_deg_fstar_phase','Delta_phase_deg_vs_cleaned_sameE','Baseline_phase_deg_sameE'} ...
    );

    % add robustness columns if available
    if ~isempty(FIXED_FREQ_HZ)
        Tout.Delta_AbsZ_fixedHz_vs_cleaned_sameE    = dAbsZ_fix;
        Tout.Delta_log10AbsZ_fixedHz_vs_cleaned_sameE = dlmz_fix;
        Tout.fixedHz_Hz = repmat(FIXED_FREQ_HZ, nFiles, 1);
    end
    Tout.Delta_AbsZ_bandMed_vs_cleaned_sameE     = dAbsZ_band;
    Tout.Delta_log10AbsZ_bandMed_vs_cleaned_sameE  = dlmz_band;
    Tout.bandLo_Hz = repmat(BAND(1), nFiles, 1);
    Tout.bandHi_Hz = repmat(BAND(2), nFiles, 1);

    % ---------------- Export CSV ----------------
    csvName = fullfile(dataDir, ['deltaZ_results_' datestr(now,'yyyymmdd_HHMMSS') '.csv']);
    try, writetable(Tout, csvName); fprintf('Results exported to: %s\n', csvName);
    catch ME, warning('CSV export failed: %s', ME.message); end

    % ---------------- Plots: Δ|Z| and Δphase (visual only; stats will run on log in compare script) ----
    condsAvail = unique(Tout.condition, 'stable');
    [~, orderIdx] = ismember(condsWanted, condsAvail); orderIdx(orderIdx==0)=[];
    condsPlot = condsAvail(orderIdx); nC = numel(condsPlot);
    cmap = containers.Map; for ci=1:nC, cmap(char(condsPlot(ci))) = clr(1+mod(ci-1,size(clr,1)),:); end; xs = 1:nC;

    % Δ|Z|
    fig1 = figure('Color','w','Name','Δ|Z| (Ω) vs cleaned (same electrode)', ...
        'Units','inches','Position',[1 1 pubWidth pubHeight],'PaperPositionMode','auto','Renderer','painters');
    ax1 = axes(fig1); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on'); yline(ax1,0,'k:','LineWidth',0.8); 
    yAll = cell(1,nC); nPer = zeros(1,nC);
    for c=1:nC
        cond = condsPlot(c);
        rows = strcmp(Tout.condition, cond);
        d = Tout.Delta_AbsZ_vs_cleaned_sameE(rows); d = d(isfinite(d));
        yAll{c} = d; nPer(c)=numel(d); if isempty(d), continue; end
        thisColor = cmap(char(cond));
        if exist('swarmchart','file')
            swarmchart(ax1, c*ones(size(d)), d, dotSize, 'filled', 'MarkerFaceAlpha',0.85,'MarkerEdgeColor','none','CData',thisColor);
        else
            jit = 0.10*(rand(size(d))-0.5);
            scatter(ax1, c+jit, d, dotSize, 'filled','MarkerFaceColor',thisColor,'MarkerEdgeColor','none');
        end
        if exist('boxchart','file')
            boxchart(ax1, c*ones(size(d)), d, 'BoxWidth',0.32,'BoxFaceColor',thisColor,'BoxFaceAlpha',0.25,'WhiskerLineColor','k','MarkerStyle','none','LineWidth',1.0);
            med = median(d,'omitnan'); plot(ax1,[c-0.15 c+0.15],[med med],'k-','LineWidth',1.6);
        end
        m = mean(d,'omitnan'); s = std(d,0,'omitnan');
        plot(ax1,[c c],[m-s m+s],'-','Color',thisColor*0.4,'LineWidth',1.0);
        plot(ax1,c,m,'o','MarkerSize',4.5,'MarkerFaceColor','k','MarkerEdgeColor','w','LineWidth',0.7);
    end
    set(ax1,'XTick',xs,'XTickLabel',cellstr(condsPlot),'FontName',fontName,'FontSize',fontSizeAx);
    xlabel(ax1,'Condition','FontName',fontName,'FontSize',fontSizeLbl);
    ylabel(ax1,sprintf('Δ|Z| (Ω) @ f*_{mag}=%.3g Hz (baseline: same electrode "cleaned")', fstar_mag),'FontName',fontName,'FontSize',fontSizeLbl,'Interpreter','tex');
    xt = get(ax1,'XTick'); yl = ylim(ax1);
    text(ax1, xt, repmat(yl(1)+0.03*range(yl),1,numel(xt)), compose('n=%d', nPer), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontName',fontName,'FontSize',fontSizeAx,'Color',[0 0 0]*0.6);
    title(ax1,'Δ|Z| (visualization; stats on log10|Z|)','FontName',fontName,'FontSize',fontSizeTtl,'FontWeight','bold');
    ylim(ax1, padlims(cell2mat(yAll), 0.08, true));

    % Δ phase
    fig2 = figure('Color','w','Name','Δ phase (deg) vs cleaned (same electrode)', ...
        'Units','inches','Position',[1 1 pubWidth pubHeight],'PaperPositionMode','auto','Renderer','painters');
    ax2 = axes(fig2); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on'); yline(ax2,0,'k:','LineWidth',0.8);
    yAll2 = cell(1,nC); nPer2 = zeros(1,nC);
    for c=1:nC
        cond = condsPlot(c);
        rows = strcmp(Tout.condition, cond);
        d = Tout.Delta_phase_deg_vs_cleaned_sameE(rows); d = d(isfinite(d));
        yAll2{c} = d; nPer2(c)=numel(d); if isempty(d), continue; end
        thisColor = cmap(char(cond));
        if exist('swarmchart','file')
            swarmchart(ax2, c*ones(size(d)), d, dotSize, 'filled','MarkerFaceAlpha',0.85,'MarkerEdgeColor','none','CData',thisColor);
        else
            jit = 0.10*(rand(size(d))-0.5);
            scatter(ax2, c+jit, d, dotSize, 'filled','MarkerFaceColor',thisColor,'MarkerEdgeColor','none');
        end
        if exist('boxchart','file')
            boxchart(ax2, c*ones(size(d)), d, 'BoxWidth',0.32,'BoxFaceColor',thisColor,'BoxFaceAlpha',0.25,'WhiskerLineColor','k','MarkerStyle','none','LineWidth',1.0);
            med = median(d,'omitnan'); plot(ax2,[c-0.15 c+0.15],[med med],'k-','LineWidth',1.6);
        end
        m = mean(d,'omitnan'); s = std(d,0,'omitnan');
        plot(ax2,[c c],[m-s m+s],'-','Color',thisColor*0.4,'LineWidth',1.0);
        plot(ax2,c,m,'o','MarkerSize',4.5,'MarkerFaceColor','k','MarkerEdgeColor','w','LineWidth',0.7);
    end
    set(ax2,'XTick',xs,'XTickLabel',cellstr(condsPlot),'FontName',fontName,'FontSize',fontSizeAx);
    xlabel(ax2,'Condition','FontName',fontName,'FontSize',fontSizeLbl);
    ylabel(ax2,sprintf('Δφ (deg) @ f*_{phase}=%.3g Hz (baseline: same electrode "cleaned")', fstar_phase),'FontName',fontName,'FontSize',fontSizeLbl,'Interpreter','tex');
    xt = get(ax2,'XTick'); yl2 = ylim(ax2);
    text(ax2, xt, repmat(yl2(1)+0.03*range(yl2),1,numel(xt)), compose('n=%d', nPer2), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontName',fontName,'FontSize',fontSizeAx,'Color',[0 0 0]*0.6);
    title(ax2,'Δ phase (deg)','FontName',fontName,'FontSize',fontSizeTtl,'FontWeight','bold');
    ylim(ax2, padlims(cell2mat(yAll2), 0.08, true));

    % Export figs
    try
        fn1 = [figOutPrefix '_dAbsZ']; fn2 = [figOutPrefix '_dphase'];
        print(fig1,[fn1 '.png'],'-dpng',sprintf('-r%d',exportDPI));
        print(fig2,[fn2 '.png'],'-dpng',sprintf('-r%d',exportDPI));
        if doPDFexport, print(fig1,[fn1 '.pdf'],'-dpdf','-painters'); print(fig2,[fn2 '.pdf'],'-dpdf','-painters'); end
        fprintf('Exported figures:\n  %s.[png/pdf]\n  %s.[png/pdf]\n', fn1, fn2);
    catch ME, warning('Figure export failed: %s', ME.message); end
end % ================== end main ==================

% ------------------------- local helpers -------------------------
function [fStart, fStop] = read_start_stop(fpath)
    raw = readcell(fpath);
    [r_sf,c_sf] = find(strcmp(raw, 'Start Frequency'), 1);
    [r_sp,c_sp] = find(strcmp(raw, 'Stop Frequency'), 1);
    assert(~isempty(r_sf)&&~isempty(r_sp), 'No Frequency col and missing Start/Stop Frequency in %s', fpath);
    fStart = raw{r_sf,c_sf+1}; fStop = raw{r_sp,c_sp+1};
    validateattributes(fStart,{'numeric'},{'scalar','positive','finite'});
    validateattributes(fStop, {'numeric'},{'scalar','positive','finite','>',fStart});
end

function [fLo,fHi] = overlap_span(recs)
    fLo = -inf; fHi = +inf;
    for r=1:numel(recs)
        fi = recs{r}.f; if isempty(fi)||~isvector(fi), continue; end
        fLo = max(fLo, min(fi)); fHi = min(fHi, max(fi));
    end
end

function v = interp_one(rec, fstar, mode)
    [fu, iu] = unique(rec.f,'stable');
    switch mode
        case 'mag',  mu = rec.mag(iu); v = interp1(fu, mu, fstar,'pchip','extrap');
        case 'lmz',  mu = reallog10(rec.mag(iu)); v = interp1(fu, mu, fstar,'pchip','extrap');
        case 'phase',pu = rec.phu(iu); v = interp1(fu, pu, fstar,'pchip','extrap');
        otherwise, error('Unknown mode: %s', mode);
    end
end

function vv = interp_per_files(recs, fstar, mode)
    vv = nan(1,numel(recs));
    for r=1:numel(recs), vv(r) = interp_one(recs{r}, fstar, mode); end
end

function [sep,fstar] = pick_fstar(M, fgrid, baseline, condsOrder)
    xb = M.(baseline);
    otherConds = condsOrder(~strcmp(condsOrder, baseline));
    Moth = [];
    for i=1:numel(otherConds), Moth = [Moth, M.(otherConds{i})]; end %#ok<AGROW>
    sep = nan(numel(fgrid),1);
    for i=1:numel(fgrid)
        vb = xb(i,:); vb = vb(isfinite(vb));
        vo = Moth(i,:); vo = vo(isfinite(vo));
        if numel(vb)>=2 && numel(vo)>=2
            mb = mean(vb,'omitnan'); mo = mean(vo,'omitnan');
            sb = std(vb,0,'omitnan'); so = std(vo,0,'omitnan');
            sp = sqrt(((numel(vb)-1)*sb^2 + (numel(vo)-1)*so^2) / max(numel(vb)+numel(vo)-2,1));
            if sp>0, sep(i) = abs(mb-mo)/sp; end
        end
    end
    [~,idx] = max(sep); fstar = fgrid(idx);
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

function y = reallog10(x)
    x = max(x, realmin); y = log10(x);
end
