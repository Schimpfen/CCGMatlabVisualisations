function Split_fitnybodeL_pub
% Exact figure layout & style to match fitnybodeL_groups4_ONEFIG
% - Files: ONLY *.xlsx in current folder (all are pulled)
% - Layout: rows = groups of 4 datasets; cols = [Bode|Z|, Phase, Nyquist]
% - Style: dashed grids, box on, markers only for data, per-row legend on Phase
% - Bode/Phase x-lims = [1e3, 1e5] Hz
% - Nyquist limits fixed to [0..5000] x [0..10000]
% - Exports PNG (300 dpi) and optional PDF (vector)

clear; close all; clc;

%% ---------------- User settings ----------------
GROUP_SIZE     = 4;        % datasets per group (row)
% Outlier removal (union of Re/Im)
DO_OUTLIER     = true;     
OUT_METHOD     = 'gesd'; OUT_ALPHA = 1;
OUT_MAX_FRAC   = 0.12; OUT_TF = 3.0; OUT_MOV_WIN_FR = 0.10;

% Frequency window on Bode/Phase
XLIM_HZ        = [1e3, 1e5];

% Optional pre-trim before plotting
F_TRIM_MIN_HZ  = [];           % [] disables
F_TRIM_MAX_HZ  = [];           % [] disables

% Aesthetics (match look)
FONT_NAME      = 'Helvetica';
rowHeightIn    = 2.6; 
figWidthIn     = 12.0;
markerSize     = 5.0; 
dataLineWidth  = 0.9; 
titleFS        = 10; 
labelFS        = 9;

% Fixed Nyquist window
NY_XMAX        = 5000;
NY_YMAX        = 10000;

% Export
PNG_DPI        = 300;                  
PNG_BASENAME   = 'Groups4_XLSX_PUB';   
DO_EXPORT_PDF  = true;

%% ---------------- Discover files (natural sort) ----------------
files = dir('*.xlsx');
if isempty(files), error('No .xlsx files found in the current folder.'); end

names = {files.name};
keys  = natpad_keys(names);
[~, ord] = sort(keys);
files = files(ord);

shortNames      = erase({files.name}',{'.xlsx'});
shortNamesLaTeX = cellfun(@latexify_filename_label, shortNames,'UniformOutput',false);
numFiles        = numel(files);
colors          = lines(max(numFiles,8));

%% ---------------- Prepare figure grid ----------------
nGroups     = ceil(numFiles / GROUP_SIZE);
figHeightIn = max(1, nGroups) * rowHeightIn;

fig = figure('Units','inches','Position',[0.5 0.5 figWidthIn figHeightIn], ...
             'Color','w','PaperPositionMode','auto','Name','Groups of 4 (XLSX, ONE FIG)');

tlo = tiledlayout(fig, nGroups, 3, 'TileSpacing','compact','Padding','compact');

%% ---------------- Iterate groups (rows) ----------------
for g = 1:nGroups
    idx = (g-1)*GROUP_SIZE + 1 : min(g*GROUP_SIZE, numFiles);

    % --- Axes for this row (exact column order) ---
    axMag = nexttile(tlo, (g-1)*3 + 1); hold(axMag,'on'); set(axMag,'XScale','log','YScale','log','FontName',FONT_NAME);
    xlabel(axMag,'f [Hz]','Interpreter','latex','FontSize',labelFS);
    ylabel(axMag,'$|Z|\,[\Omega]$','Interpreter','latex','FontSize',labelFS);
    title(axMag,sprintf('Bode |Z| — Group %d',g),'FontWeight','bold','FontSize',titleFS);
    grid(axMag,'on'); axMag.GridLineStyle='--'; axMag.Box='on'; xlim(axMag, XLIM_HZ);

    axPh  = nexttile(tlo, (g-1)*3 + 2); hold(axPh,'on'); set(axPh,'XScale','log','FontName',FONT_NAME);
    xlabel(axPh,'f [Hz]','Interpreter','latex','FontSize',labelFS);
    ylabel(axPh,'$\angle Z\,[^\circ]$','Interpreter','latex','FontSize',labelFS);
    title(axPh,sprintf('Phase — Group %d',g),'FontWeight','bold','FontSize',titleFS);
    grid(axPh,'on'); axPh.GridLineStyle='--'; axPh.Box='on'; xlim(axPh, XLIM_HZ);

    axNy  = nexttile(tlo, (g-1)*3 + 3); hold(axNy,'on'); set(axNy,'FontName',FONT_NAME);
    xlabel(axNy,'$Re(Z)\;[\Omega]$','Interpreter','latex','FontSize',labelFS);
    ylabel(axNy,'$-Im(Z)\;[\Omega]$','Interpreter','latex','FontSize',labelFS);
    title(axNy,sprintf('Nyquist — Group %d',g),'FontWeight','bold','FontSize',titleFS);
    grid(axNy,'on'); axNy.GridLineStyle='--'; axNy.Box='on';

    plotted_names = {};

    % ---- Files in this group (row) ----
    for k = idx
        fname = files(k).name;  c = colors(k,:);
        try
            [f, mag, ph_deg, ~] = read_xlsx_positional_or_texty(fname);
        catch ME
            warning('Skipping "%s": %s', fname, ME.message);
            continue;
        end

        % Clean & sanity
        good = isfinite(f) & isfinite(mag) & isfinite(ph_deg) & f>0 & mag>0;
        f = f(good); mag = mag(good); ph_deg = ph_deg(good);
        if numel(f) < 3, continue; end

        % Optional global trim
        keep = true(size(f));
        if ~isempty(F_TRIM_MIN_HZ), keep = keep & (f >= F_TRIM_MIN_HZ); end
        if ~isempty(F_TRIM_MAX_HZ), keep = keep & (f <= F_TRIM_MAX_HZ); end
        f = f(keep); mag = mag(keep); ph_deg = ph_deg(keep);
        if numel(f) < 3, continue; end

        % Complex reconstruction
        ph_rad = deg2rad(ph_deg);
        Re = mag .* cos(ph_rad);
        Im = mag .* sin(ph_rad);

        % Outlier removal (union of Re/Im)
        if DO_OUTLIER
            nPts = numel(Re); maxNum = max(1, floor(OUT_MAX_FRAC*nPts));
            wMov = max(5, round(OUT_MOV_WIN_FR*nPts)); if mod(wMov,2)==0, wMov=wMov+1; end
            idxR = robust_outliers(Re,OUT_METHOD,OUT_ALPHA,maxNum,OUT_TF,wMov);
            idxI = robust_outliers(Im,OUT_METHOD,OUT_ALPHA,maxNum,OUT_TF,wMov);
            keep = ~(idxR | idxI);
            f=f(keep); mag=mag(keep); ph_deg=ph_deg(keep); Re=Re(keep); Im=Im(keep);
        end
        if numel(f) < 3, continue; end

        % Restrict plotted Bode/Phase to window
        inWin = (f >= XLIM_HZ(1)) & (f <= XLIM_HZ(2));
        f_bp  = f(inWin); Re_bp = Re(inWin); Im_bp = Im(inWin); mag_bp = mag(inWin); ph_bp = ph_deg(inWin);
        if numel(f_bp) < 2
            f_bp=f; Re_bp=Re; Im_bp=Im; mag_bp=mag; ph_bp=ph_deg;
        end

        % Plots: markers only for measured data (exact look)
        plot(axMag, f_bp, mag_bp, 'o', 'Color',c, 'MarkerSize',markerSize, 'LineWidth',dataLineWidth, 'DisplayName',shortNamesLaTeX{k});
        plot(axPh,  f_bp, ph_bp,  'o', 'Color',c, 'MarkerSize',markerSize, 'LineWidth',dataLineWidth, 'DisplayName',shortNamesLaTeX{k});
        plot(axNy,  Re,   -Im,    'o', 'Color',c, 'MarkerSize',markerSize, 'LineWidth',dataLineWidth, 'DisplayName',shortNamesLaTeX{k});

        plotted_names{end+1} = shortNamesLaTeX{k}; %#ok<SAGROW>
    end

    % Legend per row on Phase axis
    if ~isempty(plotted_names)
        legend(axPh, plotted_names, 'Interpreter','latex', ...
           'Location','southoutside', 'NumColumns', max(2, ceil(numel(plotted_names)/2)), 'Box','off');
    end

    % Fixed Nyquist window (exact)
    set(axNy, 'XLim',[0 NY_XMAX], 'YLim',[0 NY_YMAX], ...
              'XLimMode','manual','YLimMode','manual');
end

% Re-apply fixed Nyquist to any Nyquist axes (guard vs autoscale)
axAll = findobj(fig, 'Type','axes');
for ax = axAll'
    t = get(get(ax,'Title'),'String');
    if (ischar(t)   && contains(t,'Nyquist')) || ...
       (isstring(t) && contains(t,"Nyquist"))
        set(ax,'XLim',[0 NY_XMAX], 'YLim',[0 NY_YMAX], ...
               'XLimMode','manual','YLimMode','manual');
    end
end

%% ---------------- Export ----------------
ts = datestr(now,'yyyymmdd_HHMMSS');
pngFile = sprintf('%s_%s.png', PNG_BASENAME, ts);
try
    exportgraphics(fig, pngFile, 'Resolution', PNG_DPI);
    fprintf('PNG exported: %s (dpi=%d)\n', pngFile, PNG_DPI);
catch ME
    warning('PNG export failed: %s', ME.message);
end

if DO_EXPORT_PDF
    pdfFile = sprintf('%s_%s.pdf', PNG_BASENAME, ts);
    try
        exportgraphics(fig, pdfFile, 'ContentType','vector');
        fprintf('PDF exported: %s (vector)\n', pdfFile);
    catch ME
        warning('PDF export failed: %s', ME.message);
    end
end

end % ===== main =====


%% ============================ Robust Reader ============================
function [f, mag, ph_deg, T] = read_xlsx_positional_or_texty(fname)
% Reads XLSX where columns are real OR each row is one text cell "ts,f,|Z|,phase,temp"
% Returns column vectors. Phase auto-detected deg/rad.

    % First try: normal table
    try
        Ttab = readtable(fname, 'PreserveVariableNames', true);
    catch
        Ttab = readtable(fname, 'FileType','spreadsheet','PreserveVariableNames', true);
    end

    if width(Ttab) >= 4
        f   = to_numeric(Ttab{:, 2});
        mag = to_numeric(Ttab{:, 3});
        ph  = to_numeric(Ttab{:, 4});
        if width(Ttab) >= 5, T = to_numeric(Ttab{:, 5}); else, T = nan(size(f)); end
    else
        raw = readcell(fname);
        if isempty(raw), error('Empty file: %s', fname); end
        if size(raw,2) >= 1, rows = raw(:,1); else, rows = raw(:); end

        fList=[]; magList=[]; phList=[]; Tlist=[];
        for i = 1:numel(rows)
            s = rows{i};
            if ismissing(s), continue; end
            if iscell(s), s = s{1}; end
            if isnumeric(s)
                if isempty(s) || any(~isfinite(s),'all'), continue; end
                s = string(s);
            elseif ischar(s)
                s = string(s);
            elseif ~isstring(s)
                continue;
            end
            if numel(s)~=1 || strlength(s)==0, continue; end

            parts = regexp(char(s), '[,\;\t]+', 'split');
            if numel(parts) < 4, continue; end
            fList(end+1,1)   = str2double(parts{2}); %#ok<AGROW>
            magList(end+1,1) = str2double(parts{3}); %#ok<AGROW>
            phList(end+1,1)  = str2double(parts{4}); %#ok<AGROW>
            if numel(parts) >= 5
                Tlist(end+1,1) = str2double(parts{5}); %#ok<AGROW>
            else
                Tlist(end+1,1) = NaN; %#ok<AGROW>
            end
        end

        if isempty(fList)
            error('Could not parse rows as comma-separated values in "%s".', fname);
        end

        f = fList; mag = magList; ph = phList; T = Tlist;
    end

    % Shape & clean
    f   = f(:); mag = mag(:); ph = ph(:);
    if isempty(T), T = nan(size(f)); else, T = T(:); end

    % Auto-detect radians vs degrees
    ph_deg = ph;
    if ~isempty(ph)
        if all(isfinite(ph)) && median(abs(ph),'omitnan') < 3.5
            ph_deg = rad2deg(ph);
        end
    end
end

%% ============================ Utils ====================================
function x = to_numeric(col)
    if isnumeric(col)
        x = double(col);
    elseif iscell(col)
        try
            x = cellfun(@str2double, col);
        catch
            x = nan(size(col));
            for i = 1:numel(col)
                try, x(i) = str2double(string(col{i}));
                catch, x(i) = nan; end
            end
        end
    else
        x = str2double(string(col));
    end
end

function keys = natpad_keys(names)
    padWidth = 12;
    if isstring(names), names = cellstr(names); end
    if ischar(names),   names = cellstr(names); end
    keys = cell(size(names));
    for i = 1:numel(names)
        s = names{i}; if ~ischar(s), s = char(s); end
        toks = regexp(s, '\d+|\D+', 'match');
        out = '';
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

function idx = robust_outliers(x, method, alpha, maxNum, tf, w)
    x = x(:);
    idx = false(size(x));
    n = numel(x);
    if n < 3, return; end
    switch lower(method)
        case 'gesd'
            idx = gesd_rosner(x, maxNum, alpha);
        case 'median'
            idx = isoutlier_fallback(x,'median',tf,maxNum);
        case 'quartiles'
            idx = isoutlier_fallback(x,'quartiles',tf,maxNum);
        case 'movmedian'
            try
                idx = isoutlier(x,'movmedian',w,'ThresholdFactor',tf,'MaxNumOutliers',maxNum);
            catch
                xm = movmedian(x,w,'omitnan');
                xs = movmad(x,w,1); xs(xs==0) = median(xs(xs>0));
                z = abs(x - xm) ./ (1.4826*xs);
                idx = capmax(z > tf, maxNum, z);
            end
        otherwise
            idx = gesd_rosner(x, maxNum, alpha);
    end
    function idx2 = isoutlier_fallback(x,method,tf,maxNum)
        try
            idx2 = isoutlier(x,method,'ThresholdFactor',tf,'MaxNumOutliers',maxNum);
        catch
            switch method
                case 'median'
                    med = median(x,'omitnan'); madv = mad(x,1); if madv==0, idx2=false(size(x)); return; end
                    z = abs(x - med)/(1.4826*madv);
                    idx2 = capmax(z > tf, maxNum, z);
                case 'quartiles'
                    q1 = quantile(x,0.25); q3 = quantile(x,0.75); iqrV = q3-q1;
                    lf = q1 - tf*iqrV; uf = q3 + tf*iqrV;
                    d = zeros(n,1);
                    mask = (x<lf) | (x>uf);
                    d(x<lf) = lf - x(x<lf);
                    d(x>uf) = x(x>uf) - uf;
                    idx2 = capmax(mask, maxNum, d);
            end
        end
    end
    function idx3 = capmax(mask, mmax, score)
        if nnz(mask) <= mmax, idx3 = mask; return; end
        [~,ord] = sort(score(mask),'descend');
        idx3 = false(size(mask));
        inds = find(mask);
        idx3(inds(ord(1:mmax))) = true;
    end
end

function idxOut = gesd_rosner(x, maxOut, alpha)
    x = x(:); n = numel(x);
    maxOut = min(maxOut, n-3);
    idxOut = false(n,1);
    if maxOut <= 0, return; end
    workIdx = true(n,1);
    R = zeros(maxOut,1); lambda = zeros(maxOut,1); rmOrder = zeros(maxOut,1);
    for i = 1:maxOut
        xi = x(workIdx); ni = numel(xi);
        if ni <= 2, break; end
        mu = mean(xi,'omitnan'); sd = std(xi,0,'omitnan'); if sd==0, break; end
        r = abs(xi - mu)/sd;
        [R(i), localIdx] = max(r);
        globalIdxCandidates = find(workIdx);
        rmOrder(i) = globalIdxCandidates(localIdx);
        p = 1 - alpha/(2*ni);
        try, tcrit = tinv(p, ni-2); catch, tcrit = sqrt(2) * erfcinv(alpha/ni); end
        lambda(i) = ((ni-1)*tcrit)/sqrt((ni-2 + tcrit^2)*(ni-1));
        workIdx(rmOrder(i)) = false;
    end
    m = 0; for i=1:maxOut, if R(i)>lambda(i), m=i; end, end
    if m>0, idxOut(rmOrder(1:m)) = true; end
end

function s = latexify_filename_label(s0)
    s = char(s0);
    pat = '(^|\s)(\d+)\^([\-+]?\d+)';
    [starts, ends, toks] = regexp(s, pat, 'start','end','tokens');
    if ~isempty(starts)
        out = ''; last = 1;
        for i = 1:numel(starts)
            if starts(i) > last, out = [out, s(last:starts(i)-1)]; end %#ok<AGROW>
            tok = toks{i}; if iscell(tok) && numel(tok)==1 && iscell(tok{1}), tok = tok{1}; end
            lead = tok{1}; base = tok{2}; expo = tok{3};
            out = [out, lead, '$', base, '^{', expo, '}$']; %#ok<AGROW>
            last = ends(i) + 1;
        end
        out = [out, s(last:end)];
        s = out;
    end
    s = strrep(s,'\','\textbackslash{}');
    s = strrep(s,'_','\_');
    s = strrep(s,'%','\%');
    s = strrep(s,'&','\&');
    s = strrep(s,'#','\#');
end
