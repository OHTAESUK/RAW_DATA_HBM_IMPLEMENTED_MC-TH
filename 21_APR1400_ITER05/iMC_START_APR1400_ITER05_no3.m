% ============================================================
% APR1400: Power uncertainty comparison @ x=50 (YZ)
%   - APP STD  vs  HBM STD
%   - 2 panels vertical (2x1)
%   - Fuel-based mask so guide tubes are clean holes
% ============================================================

clear; clc;

%% ===== USER SETTINGS =====
Nx = 136; Ny = 136; Nz = 38;
x0 = 50;  x0 = max(1,min(Nx,x0));

posFig = [260 60 900 950];     % 세로 2장 보기 적당히

%% ===== PATHS =====
baseDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(baseDir,'ITER05_TRIAL01');

% fuel mean (mask용)  : "통상 MC" 결과 사용
fuel_mu_file = fullfile(dataDir,'APR1400_MESH_FUEL_TEMP_MC.out');

% power std files
pow_std_APP_file = fullfile(dataDir,'APR1400_MESH_POWER_MC_STD.out');
pow_std_HBM_file = fullfile(dataDir,'APR1400_HBM_MESH_POWER_MC_STD.out');

%% ===== LOAD =====
Fmu = load_mesh(fuel_mu_file,     Nx,Ny,Nz);     % mask source
Sapp = load_mesh(pow_std_APP_file, Nx,Ny,Nz);    % APP sigma
Shbm = load_mesh(pow_std_HBM_file, Nx,Ny,Nz);    % HBM sigma

% 0 -> NaN
Fmu(Fmu==0)=NaN;
Sapp(Sapp==0)=NaN;
Shbm(Shbm==0)=NaN;

%% ===== YZ @ x0 =====
F_yz = squeeze(Fmu(x0,:,:));      % (Ny,Nz)
fuelMask_yz = isfinite(F_yz);     % guide tube holes etc.

APP_yz = squeeze(Sapp(x0,:,:));   % (Ny,Nz)
HBM_yz = squeeze(Shbm(x0,:,:));   % (Ny,Nz)

APP_yz(~fuelMask_yz) = NaN;
HBM_yz(~fuelMask_yz) = NaN;

%% ===== PLOT (2x1 vertical) =====
figure('Name',sprintf('APR1400 Power σ comparison @ x=%d (YZ)',x0), ...
       'Position',posFig);
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

% --- APP ---
nexttile;
h = imagesc(APP_yz.');
axis image; axis xy;
set(h,'AlphaData',double(isfinite(APP_yz.')), 'AlphaDataMapping','none');
set(gca,'Color','w');
xlabel('Y index'); ylabel('Z index');
title(sprintf('Power uncertainty σ (APP)  —  YZ @ x=%d',x0));
cb = colorbar; cb.Label.String = 'σ';
clim([0 0.08]);

% --- HBM ---
nexttile;
h = imagesc(HBM_yz.');
axis image; axis xy;
set(h,'AlphaData',double(isfinite(HBM_yz.')), 'AlphaDataMapping','none');
set(gca,'Color','w');
xlabel('Y index'); ylabel('Z index');
title(sprintf('Power uncertainty σ (HBM)  —  YZ @ x=%d',x0));
cb = colorbar; cb.Label.String = 'σ';
clim([0 0.08]);

%% ============================================================
% HBM temperature uncertainty @ x = 50 (YZ)
%   - Fuel T STD (HBM)
%   - Coolant T STD (HBM)
%   - Fuel-based mask (same geometry)
%% ============================================================

% ---- FILE PATHS (HBM STD only) ----
fuel_std_HBM_file = fullfile(dataDir,'APR1400_HBM_MESH_FUEL_TEMP_MC_STD.out');
bulk_std_HBM_file = fullfile(dataDir,'APR1400_HBM_MESH_BULK_TEMP_MC_STD.out');

% ---- LOAD ----
Fstd = load_mesh(fuel_std_HBM_file, Nx,Ny,Nz);          % fuel σ
Bstd = load_mesh(bulk_std_HBM_file, Nx+1,Ny+1,Nz+1);   % bulk σ

Fstd(Fstd==0)=NaN;
Bstd(Bstd==0)=NaN;

% ---- YZ @ x0 ----
Fstd_yz = squeeze(Fstd(x0,:,:));          % (Ny,Nz)
Bstd_yz = squeeze(Bstd(x0,1:Ny,1:Nz));    % (Ny,Nz)

% 동일 fuel mask 적용
FOK = isfinite(F_yz);                     % 위에서 만든 fuelMask
Fstd_yz(~FOK) = NaN;
Bstd_yz(~FOK) = NaN;

% ---- PLOT (2x1 vertical) ----
figure('Name',sprintf('APR1400 HBM temperature σ @ x=%d (YZ)',x0), ...
       'Position',[300 60 900 950]);
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

% ---- Fuel temperature σ ----
nexttile;
h = imagesc(Fstd_yz.');
axis image; axis xy;
set(h,'AlphaData',double(isfinite(Fstd_yz.')), ...
      'AlphaDataMapping','none');
set(gca,'Color','w');
xlabel('Y index'); ylabel('Z index');
title(sprintf('Fuel temperature uncertainty σ (HBM) — YZ @ x=%d',x0));
cb = colorbar; cb.Label.String='[K]';
% clim([0 10]);

% ---- Coolant temperature σ ----
nexttile;
h = imagesc(Bstd_yz.');
axis image; axis xy;
set(h,'AlphaData',double(isfinite(Bstd_yz.')), ...
      'AlphaDataMapping','none');
set(gca,'Color','w');
xlabel('Y index'); ylabel('Z index');
title(sprintf('Coolant temperature uncertainty σ (HBM) — YZ @ x=%d',x0));
cb = colorbar; cb.Label.String='[K]';
% clim([0 0.2]);

%% ============================================================
% L2 relative uncertainty (STATISTICAL, not iteration diff)
%   - APP vs HBM
%   - whole domain
%% ============================================================

% --- LOAD MEAN (통상 MC) ---
power_mu = load_mesh( ...
    fullfile(dataDir,'APR1400_MESH_POWER_MC_AVG.out'), Nx,Ny,Nz);

% --- LOAD STD ---
power_std_APP = load_mesh( ...
    fullfile(dataDir,'APR1400_MESH_POWER_MC_STD.out'), Nx,Ny,Nz);

power_std_HBM = load_mesh( ...
    fullfile(dataDir,'APR1400_HBM_MESH_POWER_MC_STD.out'), Nx,Ny,Nz);

% --- CLEAN ---
power_mu(power_mu==0)         = NaN;
power_std_APP(power_std_APP==0)= NaN;
power_std_HBM(power_std_HBM==0)= NaN;

% --- CONSISTENT MASK (fuel-based geometry) ---
fuelMask3D = isfinite(Fmu);  % 이미 위에서 읽어둔 fuel mean
mask = fuelMask3D & isfinite(power_mu) & isfinite(power_std_APP) & isfinite(power_std_HBM);

mu  = power_mu(mask);
sgA = power_std_APP(mask);
sgH = power_std_HBM(mask);

% --- RELATIVE L2 UNCERTAINTY ---
Erel_APP = sqrt(sum(sgA.^2)) / sqrt(sum(mu.^2));
Erel_HBM = sqrt(sum(sgH.^2)) / sqrt(sum(mu.^2));

fprintf('\n=== Statistical relative L2 uncertainty (POWER) ===\n');
fprintf('APP : %.6e\n', Erel_APP);
fprintf('HBM : %.6e\n', Erel_HBM);
fprintf('HBM/APP ratio = %.3f\n', Erel_HBM/Erel_APP);

%% ============================================================
% ADD: Statistical relative L2 uncertainty (FUEL, COOLANT)
%   - avg: normal MC (.out)
%   - std: HBM only
%   - mask: fuel-mean based (same geometry holes)
%% ============================================================

% ---------- (A) FUEL ----------
fuel_mu = load_mesh(fullfile(dataDir,'APR1400_MESH_FUEL_TEMP_MC.out'), Nx,Ny,Nz);
fuel_std_HBM = load_mesh(fullfile(dataDir,'APR1400_HBM_MESH_FUEL_TEMP_MC_STD.out'), Nx,Ny,Nz);

fuel_mu(fuel_mu==0) = NaN;
fuel_std_HBM(fuel_std_HBM==0) = NaN;

maskF = isfinite(fuel_mu) & isfinite(fuel_std_HBM);
NF = nnz(maskF);
if NF==0
    error('FUEL mask has zero valid cells. Check fuel_mu / fuel_std_HBM files.');
end

muF = fuel_mu(maskF);
sgF = fuel_std_HBM(maskF);

denF = sqrt(sum(muF.^2));
if denF <= eps
    error('FUEL denominator ~0 (||mu||_2). This should not happen unless mu is all zeros/NaN.');
end

Erel_FUEL_HBM = sqrt(sum(sgF.^2)) / denF;

% ---------- (B) COOLANT ----------
bulk_mu = load_mesh(fullfile(dataDir,'APR1400_MESH_BULK_TEMP_MC.out'), Nx+1,Ny+1,Nz+1);
bulk_std_HBM = load_mesh(fullfile(dataDir,'APR1400_HBM_MESH_BULK_TEMP_MC_STD.out'), Nx+1,Ny+1,Nz+1);

bulk_mu(bulk_mu==0) = NaN;
bulk_std_HBM(bulk_std_HBM==0) = NaN;

% bulk를 fuel 격자와 같은 크기로 맞춤 (IMPORTANT)
muB3 = bulk_mu(1:Nx,1:Ny,1:Nz);
sgB3 = bulk_std_HBM(1:Nx,1:Ny,1:Nz);

% fuel mean(Fmu) 기반 geometry mask를 bulk에도 동일 적용
% (이미 위에서 Fmu 읽어둔 게 있으면 fuel_mu 대신 그걸 써도 됨)
fuelMask3D = isfinite(fuel_mu);        % (Nx,Ny,Nz)

muB3(~fuelMask3D) = NaN;
sgB3(~fuelMask3D) = NaN;

maskB = isfinite(muB3) & isfinite(sgB3);
NB = nnz(maskB);
if NB==0
    error('COOLANT mask has zero valid cells. Likely mask mismatch (bulk indexing or fuelMask3D).');
end

muB = muB3(maskB);
sgB = sgB3(maskB);

denB = sqrt(sum(muB.^2));
if denB <= eps
    error('COOLANT denominator ~0 (||mu||_2). If bulk_mu is delta-T or near-zero field, relative L2 is undefined.');
end

Erel_BULK_HBM = sqrt(sum(sgB.^2)) / denB;

fprintf('\n=== Statistical relative L2 uncertainty (HBM) ===\n');
fprintf('Fuel temperature   : %.6e   (N=%d)\n', Erel_FUEL_HBM, NF);
fprintf('Coolant temperature: %.6e   (N=%d)\n', Erel_BULK_HBM, NB);


%% ===== FUNCTION =====
function A = load_mesh(path,Nx,Ny,Nz)
    txt  = fileread(path);
    vals = sscanf(txt,'%f');
    nExp = Nx*Ny*Nz;
    assert(numel(vals)==nExp, 'Mesh size mismatch: %s (got %d, expected %d)', ...
           path, numel(vals), nExp);
    A = reshape(vals,[Ny,Nx,Nz]);
    A = permute(A,[2 1 3]);   % (Nx,Ny,Nz)
end
