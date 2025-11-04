function summary = kk_test_fullband(dataDir, varargin)
% KK_TEST_FULLBAND  Full-band Kramers–Kronig check with band-validity map (>=100 Hz).
% Matches your loader: col1=Re(Ω), col2=Im(Ω), col3(optional)=f(Hz) or Start/Stop labels.
%
% Options (name/value):
%   'F_MIN_HZ'     (default 100)   % analyze f >= this
%   'hfFraction'   (default 0.20)  % top fraction of points to estimate Rs/L/Cdl
%   'winFrac'      (default 0.15)  % rolling-window fraction for local NRMSE
%   'thrGood'      (default 10)    % local NRMSE < thrGood  => Valid
%   'thrMarg'      (default 20)    % thrGood–thrMarg       => Marginal; else Invalid
%   'makePlots'    (default true)
%   'exportPDF'    (default true)

if nargin < 1 || isempty(dataDir), dataDir = pwd; end
p = inputParser;
addParameter(p,'F_MIN_HZ',100,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'hfFraction',0.20,@(x)isnumeric(x)&&x>0&&x<=0.6);
addParameter(p,'winFrac',0.15,@(x)isnumeric(x)&&x>0&&x<=0.6);
addParameter(p,'thrGood',10,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'thrMarg',20,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'makePlots',true,@islogical);
addParameter(p,'exportPDF',true,@islogical);
parse(p,varargin{:});
F_MIN_HZ  = p.Results.F_MIN_HZ;
hfFraction= p.Results.hfFraction;
winFrac   = p.Results.winFrac;
thrGood   = p.Results.thrGood;
thrMarg   = p.Results.thrMarg;
makePlots = p.Results.makePlots;
exportPDF = p.Results.exportPDF;

files = dir(fullfile(dataDir,'*.xlsx'));
maskSkip = startsWith({files.name}, 'sim_results_') | ...
           startsWith({files.name}, 'fit_params_')  | ...
           startsWith({files.name}, 'deltaZ_results_') | ...
           startsWith({files.name}, 'kk_summary_') | ...
           startsWith({files.name}, 'kk1k_summary_');
files = files(~maskSkip);
assert(~isempty(files),'No .xlsx files found (after skipping exports).');

rows = {};
ts = datestr(now,'yyyymmdd_HHMMSS');
csvPath = fullfile(dataDir, ['kk_fullband_summary_' ts '.csv']);

for k = 1:numel(files)
    fpath = fullfile(files(k).folder, files(k).name);
    try
        [f, Re, Im] = read_eis_like_yours(fpath);
        keep = isfinite(f)&isfinite(Re)&isfinite(Im)&(f>0);
        if F_MIN_HZ>0, keep = keep & (f>=F_MIN_HZ); end
        f = f(keep); Re = Re(keep); Im = Im(keep);

        if numel(f) < 20
            warning('Too few usable points after trim: %s', files(k).name);
            continue;
        end

        % ---------- KK core on available band ----------
        out = kk_core_full(f, Re, Im, hfFraction, winFrac, thrGood, thrMarg, makePlots);

        % ---------- Plot export ----------
        if makePlots
            base = fullfile(dataDir, ['kk_full_' strip_ext(files(k).name) '_' ts]);
            print(gcf, [base '.png'], '-dpng','-r300');
            if exportPDF, print(gcf, [base '.pdf'], '-dpdf','-painters'); end
        end

        % ---------- Band segments to strings ----------
        sValid   = segs_to_str(out.freq, out.bandFlag==1);
        sMarg    = segs_to_str(out.freq, out.bandFlag==0);
        sInvalid = segs_to_str(out.freq, out.bandFlag==-1);

        rows(end+1,:) = { files(k).name, numel(f), min(f), max(f), ...
                          out.Rs, out.L, out.Cdl, ...
                          out.metrics.NRMSE_Re_pct, out.metrics.NRMSE_Im_pct, ...
                          sValid, sMarg, sInvalid }; %#ok<AGROW>
    catch ME
        warning('KK failed for %s: %s', files(k).name, ME.message);
    end
end

if isempty(rows)
    warning('No successful KK results to summarize.'); summary = table(); return;
end

summary = cell2table(rows, 'VariableNames', ...
    {'file','N','fmin_Hz','fmax_Hz','Rs_ohm','L_H','Cdl_F', ...
     'Global_NRMSE_Re_pct','Global_NRMSE_Im_pct', ...
     'ValidBands_fHz','MarginalBands_fHz','InvalidBands_fHz'});

try
    writetable(summary, csvPath);
    fprintf('Wrote KK full-band summary: %s\n', csvPath);
catch ME
    warning('Failed to write CSV: %s', ME.message);
end
end

%% ======= Core full-band KK with local validity map =======
function out = kk_core_full(f_Hz, ReZ, ImZ, hfFraction, winFrac, thrGood, thrMarg, makePlots)

% Clean/sort, merge duplicate f by median
[f, Re, Im] = clean_sort(f_Hz(:), ReZ(:), ImZ(:));
N = numel(f); w = 2*pi*f;

% ----- HF parasitics (fit Rs, L, Cdl on top fraction) -----
nHF   = max(12, round(hfFraction*N));
idxHF = (N-nHF+1):N;

% Rs: robust median of Re at HF
Rs_est = median(Re(idxHF), 'omitnan');

% Fit Im ≈ a + ωL − 1/(ω Cdl) on HF slice (a soaks residual fixture offset)
w_hf  = w(idxHF); im_hf = Im(idxHF);
A = [ones(numel(w_hf),1), w_hf, -1./w_hf];    % [a, L, 1/Cdl term]
theta = A \ im_hf;
a_est   = theta(1);
L_est   = max(theta(2), 0);                   % physical constraint
Cdl_est = 1 / max(theta(3), eps);

% Remove parasitics within band
Re_p = Re - Rs_est;
Im_p = Im - (a_est + w.*L_est - 1./(w.*Cdl_est));

% ----- KK transforms (principal value, trapezoid in ω) -----
xi = w(:);
dxi = trapz_steps(xi);
Xi  = xi.';  Wi = xi;
den = (Xi.^2) - (Wi.^2);
den(1:numel(xi)+1:end) = NaN;   % explicit diagonal removal (PV)

% Re' from Im'
Re_p_KK = (2/pi) * nansum( (Xi .* Im_p(:)')./den .* dxi(:).', 2 );
% Im' from Re'
Im_p_KK = -(2/pi) * Wi .* nansum( (Re_p(:)')./den .* dxi(:).', 2 );

% Re-add parasitics for comparison
Re_KK = Re_p_KK + Rs_est;
Im_KK = Im_p_KK + (a_est + w.*L_est - 1./(w.*Cdl_est));

% ----- Global metrics -----
res_Re = Re - Re_KK;
res_Im = Im - Im_KK;
rmse_Re  = sqrt(mean(res_Re.^2,'omitnan'));
rmse_Im  = sqrt(mean(res_Im.^2,'omitnan'));
rng_Re   = max(Re)-min(Re); if rng_Re < 10*eps, rng_Re = norm(Re,2)+eps; end
rng_Im   = max(Im)-min(Im); if rng_Im < 10*eps, rng_Im = norm(Im,2)+eps; end
nrmse_Re = 100*rmse_Re / rng_Re;
nrmse_Im = 100*rmse_Im / rng_Im;

% ----- Local (rolling) NRMSE vs frequency -----
winN = max(7, round(winFrac*N));   % window length in points on log-grid
locNrmse = local_nrmse(log10(f), [Re Im], [Re_KK Im_KK], winN);  % returns [locNrmseRe, locNrmseIm]
locNrmseRe = locNrmse(:,1);
locNrmseIm = locNrmse(:,2);

% Band flags per frequency (more conservative: take max of Re/Im local NRMSE)
locMax = max(locNrmseRe, locNrmseIm);
bandFlag = zeros(N,1);                       % 0 = marginal
bandFlag(locMax <  thrGood) = 1;             % valid
bandFlag(locMax >  thrMarg) = -1;            % invalid

out = struct('freq',f,'w',w,'ReZ',Re,'ImZ',Im, ...
             'ReZ_KK',Re_KK,'ImZ_KK',Im_KK, ...
             'res_Re',res_Re,'res_Im',res_Im, ...
             'Rs',Rs_est,'L',L_est,'Cdl',Cdl_est, ...
             'a_offset_Im',a_est, ...
             'metrics',struct('RMSE_Re',rmse_Re,'NRMSE_Re_pct',nrmse_Re, ...
                              'RMSE_Im',rmse_Im,'NRMSE_Im_pct',nrmse_Im), ...
             'locNrmseRe',locNrmseRe,'locNrmseIm',locNrmseIm,'bandFlag',bandFlag);

% ----- Plots -----
if makePlots
    figure('Color','w','Units','inches','Position',[1 1 12 8]);

    % A: Re overlay
    subplot(3,2,1); hold on; box on; grid on;
    plot(f, Re, 'o','MarkerSize',4);
    plot(f, Re_KK, '-', 'LineWidth',1.5);
    set(gca,'XScale','log'); xlabel('Frequency (Hz)'); ylabel('Re(Z) (\Omega)');
    title(sprintf('Re: obs vs KK  (Global NRMSE=%.2f%%)', nrmse_Re));
    legend('Obs','KK','Location','best');

    % B: Im overlay
    subplot(3,2,2); hold on; box on; grid on;
    plot(f, Im, 'o','MarkerSize',4);
    plot(f, Im_KK, '-', 'LineWidth',1.5);
    set(gca,'XScale','log'); xlabel('Frequency (Hz)'); ylabel('Im(Z) (\Omega)');
    title(sprintf('Im: obs vs KK  (Global NRMSE=%.2f%%)', nrmse_Im));
    legend('Obs','KK','Location','best');

    % C: Re residuals
    subplot(3,2,3); hold on; box on; grid on;
    semilogx(f, res_Re, '-'); yline(0,'k:');
    set(gca,'XScale','log')

    xlabel('Frequency (Hz)'); ylabel('Residual Re (\Omega)');
    title('Residual Re(obs) - Re(KK)');

    % D: Im residuals
    subplot(3,2,4); hold on; box on; grid on;
    semilogx(f, res_Im, '-'); yline(0,'k:');
    set(gca,'XScale','log')

    xlabel('Frequency (Hz)'); ylabel('Residual Im (\Omega)');
    title('Residual Im(obs) - Im(KK)');

    % E: Local NRMSE bands (Re)
    subplot(3,2,5); hold on; box on; grid on;
    semilogx(f, locNrmseRe, '-','LineWidth',1.2);
    set(gca,'XScale','log')

    yline(thrGood,'g--','<10% valid'); yline(thrMarg,'r--','>20% invalid');
    xlabel('Frequency (Hz)'); ylabel('Local NRMSE Re (%)');
    title(sprintf('Local NRMSE (Re) | window=%d pts (%.0f%%)', winN, 100*winFrac));

    % F: Validity map (max of Re/Im)
    subplot(3,2,6); hold on; box on; grid on;
    semilogx(f, locMax, '-','LineWidth',1.2);
    set(gca,'XScale','log')

    yline(thrGood,'g--'); yline(thrMarg,'r--');
    xlabel('Frequency (Hz)'); ylabel('Local NRMSE max(Re,Im) (%)');
    title('Band validity: <10% Valid, 10–20% Marginal, >20% Invalid');

    % Shade regions for readability
    yl = ylim; 
    x = f(:);
    maskValid   = out.bandFlag==1;
    maskInvalid = out.bandFlag==-1;
    shade_regions(x, maskValid, [0.85 1.0 0.85], yl);     % light green
    shade_regions(x, maskInvalid, [1.0 0.85 0.85], yl);   % light red
    uistack(findobj(gca,'Type','line'),'top');
end
end

%% ======= Helpers =======
function [f, Re, Im] = read_eis_like_yours(fpath)
M = readmatrix(fpath);
assert(size(M,2) >= 2, 'Need at least 2 cols (Re, Im): %s', fpath);
Re = M(:,1); Im = M(:,2);
ok = isfinite(Re)&isfinite(Im);
Re = Re(ok); Im = Im(ok);

if size(M,2) >= 3
    fc = M(ok,3);
    g2 = isfinite(fc) & (fc>0);
    if nnz(g2) >= 4
        f  = fc(g2); Re = Re(g2); Im = Im(g2);
    else
        [fStart,fStop] = read_start_stop(fpath);
        f = logspace(log10(fStart), log10(fStop), numel(Re)).';
    end
else
    [fStart,fStop] = read_start_stop(fpath);
    f = logspace(log10(fStart), log10(fStop), numel(Re)).';
end
end

function [fStart, fStop] = read_start_stop(fpath)
raw = readcell(fpath);
[r_sf,c_sf] = find(strcmp(raw,'Start Frequency'),1);
[r_sp,c_sp] = find(strcmp(raw,'Stop Frequency'),1);
assert(~isempty(r_sf)&&~isempty(r_sp), ...
    'No Frequency col and missing Start/Stop Frequency labels in %s', fpath);
fStart = raw{r_sf, c_sf+1};
fStop  = raw{r_sp, c_sp+1};
validateattributes(fStart,{'numeric'},{'scalar','positive','finite'});
validateattributes(fStop, {'numeric'},{'scalar','positive','finite','>',fStart});
end

function [f, Re, Im] = clean_sort(f, Re, Im)
mask = isfinite(f)&isfinite(Re)&isfinite(Im)&(f>0);
f=f(mask); Re=Re(mask); Im=Im(mask);
[f, i] = sort(f,'ascend'); Re=Re(i); Im=Im(i);
[uf,~,ic] = unique(f);
if numel(uf) < numel(f)
    Re = accumarray(ic, Re, [], @median);
    Im = accumarray(ic, Im, [], @median);
    f  = uf;
end
end

function dxi = trapz_steps(xi)
xi = xi(:);
N = numel(xi);
if N<2, dxi = zeros(N,1); return; end
d = diff(xi);
dxi = zeros(N,1);
dxi(1)     = d(1)/2;
dxi(2:N-1) = (d(1:end-1)+d(2:end))/2;
dxi(N)     = d(end)/2;
end

function L = contiguous_segments(x, mask)
% Return Nx2 segments [x_start, x_end] where mask==true and contiguous in index
mask = mask(:) & isfinite(x(:));
if ~any(mask), L = []; return; end
idx = find(mask);
splits = [1; 1+find(diff(idx) > 1); numel(idx)+1];
L = [];
for s = 1:numel(splits)-1
    block = idx(splits(s):splits(s+1)-1);
    L(end+1,:) = [x(block(1)), x(block(end))]; %#ok<AGROW>
end
end

function S = segs_to_str(x, mask)
segs = contiguous_segments(x, mask);
if isempty(segs), S = ''; return; end
buf = strings(size(segs,1),1);
for i=1:size(segs,1)
    buf(i) = sprintf('[%.3g–%.3g]', segs(i,1), segs(i,2));
end
S = strjoin(buf, ', ');
end

function shade_regions(x, mask, rgb, yl)
% Shade contiguous regions defined by mask across x (semilog axis compatible)
segs = contiguous_segments(x, mask);
for i=1:size(segs,1)
    xx = segs(i,:);
    patch([xx(1) xx(2) xx(2) xx(1)], [yl(1) yl(1) yl(2) yl(2)], rgb, ...
          'EdgeColor','none','FaceAlpha',0.25,'HandleVisibility','off');
end
end

function loc = local_nrmse(logf, Yobs2, Ykk2, winN, varargin)
% LOCAL_NRMSE  Rolling NRMSE with global-range normalization, floor, and smoothing.
% Inputs:
%   logf   : log10(frequency) vector [N x 1] (only used for alignment)
%   Yobs2  : [N x 2] observed {Re, Im}
%   Ykk2   : [N x 2] KK-predicted {Re, Im}
%   winN   : window length in points (along index space)
% Options (name/value):
%   'floorFrac' : floor as fraction of global range (default 0.05 = 5%)
%   'smoothN'   : movmedian window length for smoothing (default = max(5,round(winN/2)))

p = inputParser;
addParameter(p,'floorFrac',0.05,@(x)isnumeric(x)&&x>0&&x<1);
addParameter(p,'smoothN',[],@(x)isnumeric(x)&&x>=1);
parse(p,varargin{:});
floorFrac = p.Results.floorFrac;
smoothN   = p.Results.smoothN;

N = size(Yobs2,1);
winN = max(3, min(round(winN), N));
if isempty(smoothN), smoothN = max(5, round(winN/2)); end

% Global ranges (per channel) with safe fallbacks
glob_rng = zeros(1,2);
for c = 1:2
    rg = max(Yobs2(:,c)) - min(Yobs2(:,c));
    if rg < 10*eps
        rg = norm(Yobs2(:,c),2) + eps;  % fallback if essentially flat globally
    end
    glob_rng(c) = rg;
end
floorRe = floorFrac * glob_rng(1);
floorIm = floorFrac * glob_rng(2);

loc = nan(N,2);
half = floor(winN/2);

for i = 1:N
    i0 = max(1, i-half);
    i1 = min(N, i+half);
    idx = i0:i1;

    for c = 1:2
        y  = Yobs2(idx,c);
        yk = Ykk2(idx,c);

        rmse = sqrt(mean((y - yk).^2, 'omitnan'));

        % Robust local spread with floor tied to global range
        rg_local = iqr(y);                % robust spread (less sensitive than max-min)
        if ~isfinite(rg_local) || rg_local < 10*eps
            rg_local = 0;                 % treat as flat window
        end
        denom = max(rg_local, (c==1)*floorRe + (c==2)*floorIm);

        loc(i,c) = 100 * rmse / denom;
    end
end

% Optional smoothing to reduce salt-and-pepper artifacts
for c = 1:2
    loc(:,c) = movmedian(loc(:,c), smoothN, 'omitnan');
end
end


function s = strip_ext(fn)
[~,s] = fileparts(fn);
end
