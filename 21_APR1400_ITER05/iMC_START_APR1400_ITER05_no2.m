% ============================================================
% APR1400 iMC tally visualization (NO batch, direct .out files)
%  - 2D z-averaged maps:   Power / Fuel T / Coolant T   (1x3, horizontal)
%  - Y–Z slice @ x = 50:   Power / Fuel T / Coolant T   (3x1, vertical)
%  - Fuel-based mask applied to all (guide tube holes show cleanly)
% ============================================================

clear; clc;

%% ================= USER SETTINGS =================
Nx = 136; Ny = 136; Nz = 38;

x0 = 50;                              % Y–Z slice location (cell index)
x0 = max(1, min(Nx, x0));             % clamp

% window placements (easy to grab)
pos2D = [120 140 1800 560];           % 2D (wide)
posYZ = [260 10  900  1050];          % YZ (tall, but not crazy)

%% ================= PATHS ==========================
baseDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(baseDir, 'ITER05_TRIAL01');

powerFile = fullfile(dataDir, 'APR1400_MESH_POWER_MC_AVG.out');
fuelFile  = fullfile(dataDir, 'APR1400_MESH_FUEL_TEMP_MC.out');
bulkFile  = fullfile(dataDir, 'APR1400_MESH_BULK_TEMP_MC.out');

%% ================= LOAD ===========================
% Power & Fuel : (Nx,Ny,Nz)
power = load_mesh(powerFile, Nx, Ny, Nz);
fuel  = load_mesh(fuelFile,  Nx, Ny, Nz);

% Bulk : (Nx+1,Ny+1,Nz+1)
bulk  = load_mesh(bulkFile,  Nx+1, Ny+1, Nz+1);

% Convert "0" to NaN where appropriate
power(power==0) = NaN;
fuel(fuel==0)   = NaN;
bulk(bulk==0)   = NaN;

%% ================= FUEL-BASED GEOMETRY MASK =================
fuelMask3D = isfinite(fuel);          % (Nx,Ny,Nz)

% power는 fuel 비활성 영역 같이 빵구
power(~fuelMask3D) = NaN;

% bulk은 (Nx+1,Ny+1,Nz+1)라서 (Nx,Ny,Nz) mask를 bulk에 매핑
bulkMask3D = false(size(bulk));
bulkMask3D(1:Nx, 1:Ny, 1:Nz) = fuelMask3D;
bulk(~bulkMask3D) = NaN;

%% ============================================================
% (1) 2D MAPS  <z>  (Power/Fuel/Bulk)  1x3
%% ============================================================
mapP = mean(power, 3, 'omitnan');
mapF = mean(fuel,  3, 'omitnan');
mapB = mean(bulk(1:Nx,1:Ny,1:Nz), 3, 'omitnan');   % bulk을 Nx,Ny로 맞춤

% Normalize power map
mapP = mapP / mean(mapP,'all','omitnan');

% fuel 기반 2D 마스크
fuelMask2D = isfinite(mapF);
mapP(~fuelMask2D) = NaN;
mapB(~fuelMask2D) = NaN;

figure('Name','APR1400 2D z-avg maps (fuel-masked)', 'Position', pos2D);
tiledlayout(1,3,"Padding","compact","TileSpacing","compact");

% ---- Power <z> ----
nexttile;
h = imagesc(mapP);
axis image; axis xy;
set(h,'AlphaData', double(isfinite(mapP)), 'AlphaDataMapping','none');
set(gca,'Color','w');
title('Power  ⟨·⟩_z (normalized)');
xlabel('X index'); ylabel('Y index');
cb = colorbar; cb.Label.String = 'Normalized power';
clim([0.0 1.8]);

% ---- Fuel T <z> ----
nexttile;
h = imagesc(mapF);
axis image; axis xy;
set(h,'AlphaData', double(isfinite(mapF)), 'AlphaDataMapping','none');
set(gca,'Color','w');
title('Fuel T  ⟨·⟩_z');
xlabel('X index'); ylabel('Y index');
cb = colorbar; cb.Label.String = '[K]';
clim([580.0 740.0]);

% ---- Coolant T <z> ----
nexttile;
h = imagesc(mapB);
axis image; axis xy;
set(h,'AlphaData', double(isfinite(mapB)), 'AlphaDataMapping','none');
set(gca,'Color','w');
title('Coolant T  ⟨·⟩_z');
xlabel('X index'); ylabel('Y index');
cb = colorbar; cb.Label.String = '[K]';
clim([564.0 567.0]);

%% ============================================================
% (2) Y–Z @ x = x0  (Power/Fuel/Bulk)  3x1  (VERTICAL!!)
%% ============================================================
P_yz = squeeze(power(x0,:,:));                 % (Ny,Nz)
F_yz = squeeze(fuel(x0,:,:));                  % (Ny,Nz)
B_yz = squeeze(bulk(x0,1:Ny,1:Nz));            % (Ny,Nz)로 맞춤

% Normalize power slice (slice 내부 평균)
P_yz = P_yz / mean(P_yz(:),'omitnan');

% fuel mask 적용(2D)
fuelMask_yz = isfinite(F_yz); 
P_yz(~fuelMask_yz) = NaN;
B_yz(~fuelMask_yz) = NaN;

figure('Name',sprintf('APR1400 Y–Z @ x=%d (fuel-masked)',x0), 'Position', posYZ);
tiledlayout(3,1,"Padding","compact","TileSpacing","compact");

% ---- Power YZ ----
nexttile;
h = imagesc(P_yz.');
axis image; axis xy;
set(h,'AlphaData', double(isfinite(P_yz.')), 'AlphaDataMapping','none');
set(gca,'Color','w');
title(sprintf('Power (YZ @ x=%d, normalized)',x0));
xlabel('Y index'); ylabel('Z index');
cb = colorbar; cb.Label.String = 'Normalized power';
clim([0.0 2.2]);

% ---- Fuel YZ ----
nexttile;
h = imagesc(F_yz.');
axis image; axis xy;
set(h,'AlphaData', double(isfinite(F_yz.')), 'AlphaDataMapping','none');
set(gca,'Color','w');
title(sprintf('Fuel T (YZ @ x=%d)',x0));
xlabel('Y index'); ylabel('Z index');
cb = colorbar; cb.Label.String = '[K]';
clim([580.0 800.0]);

% ---- Coolant YZ ----
nexttile;
h = imagesc(B_yz.');
axis image; axis xy;
set(h,'AlphaData', double(isfinite(B_yz.')), 'AlphaDataMapping','none');
set(gca,'Color','w');
title(sprintf('Coolant T (YZ @ x=%d)',x0));
xlabel('Y index'); ylabel('Z index');
cb = colorbar; cb.Label.String = '[K]';
clim([564.0 569.0]);

%% ================= SAVE FUEL GEOMETRY MASK =================
% Purpose:
%   - Construct a clean, reference fuel geometry mask
%   - Guide tubes / water holes are FALSE
%   - Fuel regions are TRUE
%   - This mask is reused across TH / wo_TH / HBM ON-OFF scripts
%
% Definition:
%   Fuel region is inferred from fuel temperature mesh
%   (non-finite or zero-temperature cells are treated as non-fuel)
% 
% fuelMask3D_ref = isfinite(fuel);
% 
% % Optional sanity check: enforce strict binary
% fuelMask3D_ref = logical(fuelMask3D_ref);
% 
% % Save once (overwrite is intended)
% maskFile = fullfile(baseDir,'APR1400_fuelMask3D_ref.mat');
% save(maskFile,'fuelMask3D_ref');

%% ================= FUNCTIONS =====================
function A = load_mesh(path, Nx, Ny, Nz)
    txt  = fileread(path);
    vals = sscanf(txt, '%f');

    nExp = Nx*Ny*Nz;
    assert(numel(vals) == nExp, 'Mesh size mismatch: %s (got %d, expected %d)', path, numel(vals), nExp);

    A = reshape(vals, [Ny, Nx, Nz]);   % file ordering assumption
    A = permute(A, [2 1 3]);           % -> (Nx,Ny,Nz)
end
