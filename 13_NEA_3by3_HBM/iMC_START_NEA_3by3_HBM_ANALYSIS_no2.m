% =====================================================================
% HBM vs REF comparison (Y–Z slices @ fixed X)
%  - |HBM - REF|
%  - propagated uncertainty: sqrt(σ_HBM^2 + σ_REF^2)
%  - SAME clim for diff and sigma (per quantity) for fair comparison
% =====================================================================

clear; clc;

%% ====================== USER SETTINGS ======================
HBM_index = 3;
x_list    = [8, 25, 42];

Nx = 51; Ny = 51; Nz = 20;

% --- original clim upper bounds (from your script)
CLIM.power.diff = [0 0.03];
CLIM.power.sig  = [0 0.025];

CLIM.fuel.diff  = [0 20];
CLIM.fuel.sig   = [0 10];

CLIM.bulk.diff  = [0 0.5];
CLIM.bulk.sig   = [0 0.2];

% --- NEW: common clim per quantity (use the larger upper bound)
CLIMC.power = [0 max(CLIM.power.diff(2), CLIM.power.sig(2))];   % [0 0.03]
CLIMC.fuel  = [0 max(CLIM.fuel.diff(2),  CLIM.fuel.sig(2))];    % [0 20]
CLIMC.bulk  = [0 max(CLIM.bulk.diff(2),  CLIM.bulk.sig(2))];    % [0 0.5]

%% ====================== PATH HANDLING ======================
if HBM_index < 10
    pathHead = sprintf('./NEA_3by3_HBM_N_B_50_TRIAL0%d',HBM_index);
else
    pathHead = sprintf('./NEA_3by3_HBM_N_B_50_TRIAL%d',HBM_index);
end

files.power.avg = '/OECD_NEA_3by3_MESH_POWER_MC_AVG.out';
files.power.std = '/OECD_NEA_3by3_HBM_MESH_POWER_MC_STD.out';

files.fuel.avg  = '/OECD_NEA_3by3_MESH_FUEL_TEMP_MC.out';
files.fuel.std  = '/OECD_NEA_3by3_HBM_MESH_FUEL_TEMP_MC_STD.out';

files.bulk.avg  = '/OECD_NEA_3by3_MESH_BULK_TEMP_MC.out';
files.bulk.std  = '/OECD_NEA_3by3_HBM_MESH_BULK_TEMP_MC_STD.out';

%% ====================== LOAD HBM DATA ======================
power_avg = load_mesh(fullfile(pathHead,files.power.avg),Nx,Ny,Nz);
power_std = load_mesh(fullfile(pathHead,files.power.std),Nx,Ny,Nz);

fuel_avg  = load_mesh(fullfile(pathHead,files.fuel.avg),Nx,Ny,Nz);
fuel_std  = load_mesh(fullfile(pathHead,files.fuel.std),Nx,Ny,Nz);

bulk_avg  = load_mesh(fullfile(pathHead,files.bulk.avg),Nx+1,Ny+1,Nz+1);
bulk_std  = load_mesh(fullfile(pathHead,files.bulk.std),Nx+1,Ny+1,Nz+1);

%% ====================== LOAD REFERENCE (AVG + STD) ======================
load('NEA_3by3_REF_DATA.mat', ...
     'power_avg_ref','power_std_ref', ...
     'Tfuel_avg_ref','Tfuel_std_ref', ...
     'Tbulk_avg_ref','Tbulk_std_ref');

%% ====================== GLOBAL PLOT STYLE ======================
set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesFontSize',13);
set(groot,'defaultColorbarFontSize',11);

%% ====================== MAIN LOOP ======================
for ix = 1:numel(x_list)
    x0 = x_list(ix);

    % ---------------- POWER ----------------
    diff_yz = abs_slice(power_avg, power_avg_ref, x0);
    sig_yz  = sqrt( slice_x(power_std,x0).^2 + slice_x(power_std_ref,x0).^2 );
    plot_YZ_block_sameclim(diff_yz, sig_yz, x0, 'Power', CLIMC.power);

    % ------------- FUEL TEMPERATURE -------------
    diff_yz = abs_slice(fuel_avg, Tfuel_avg_ref, x0);
    sig_yz  = sqrt( slice_x(fuel_std,x0).^2 + slice_x(Tfuel_std_ref,x0).^2 );
    plot_YZ_block_sameclim(diff_yz, sig_yz, x0, 'Fuel temperature [K]', CLIMC.fuel);

    % ------------- BULK TEMPERATURE (staggered) -------------
    xb = min(x0, size(bulk_avg,1));
    diff_yz = abs_slice(bulk_avg, Tbulk_avg_ref, xb);
    sig_yz  = sqrt( slice_x(bulk_std,xb).^2 + slice_x(Tbulk_std_ref,xb).^2 );
    plot_YZ_block_sameclim(diff_yz, sig_yz, xb, 'Coolant temperature [K]', CLIMC.bulk);
end

%% ====================== FUNCTIONS ======================

function A = load_mesh(path,Nx,Ny,Nz)
    txt  = fileread(path);
    vals = sscanf(txt,'%f');
    assert(numel(vals)==Nx*Ny*Nz,'Size mismatch: %s',path);
    A = reshape(vals,[Ny,Nx,Nz]);
    A = permute(A,[2 1 3]);   % (Nx,Ny,Nz)
end

function yz = slice_x(A,x0)
    yz = squeeze(A(x0,:,:));
end

function yz = abs_slice(A, Aref, x0)
    yz = abs(squeeze(A(x0,:,:)) - squeeze(Aref(x0,:,:)));
end

function plot_YZ_block_sameclim(diff_yz, sig_yz, x0, label, clim_common)

    mask = isfinite(diff_yz) & isfinite(sig_yz);
    diff_yz(~mask) = NaN;
    sig_yz(~mask)  = NaN;

    f = figure('Name',sprintf('%s |Δ| vs σprop @ x=%d (same clim)',label,x0), ...
               'Position',[360 120 850 700]);

    t = tiledlayout(f,2,1,"Padding","compact","TileSpacing","compact");

    % -------- |HBM - REF|
    nexttile;
    imagesc(diff_yz.');
    axis image xy;
    title(['|HBM - REF| : ' label],'FontSize',15);
    cb = colorbar; cb.Label.String = label;
    clim(clim_common);

    % -------- propagated uncertainty
    nexttile;
    imagesc(sig_yz.');
    axis image xy;
    title(['Propagated uncertainty : ' label],'FontSize',15);
    cb = colorbar; cb.Label.String = label;
    clim(clim_common);
end
