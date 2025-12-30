% =================================================================== %
% iMC에서 계산된 출력밀도 구배정보 + 핵연료/냉각재 온도 구배정보 읽고 출력 %
% =================================================================== %
HBM_index = 3;
if(HBM_index < 10)
    pathHead = strcat('./NEA_3by3_HBM_N_B_50_TRIAL0',num2str(HBM_index));
else
    pathHead = strcat('./NEA_3by3_HBM_N_B_50_TRIAL', num2str(HBM_index));
end

powerFile = '/OECD_NEA_3by3_MESH_POWER_MC_AVG.out'; % MC_AVG 파일
TfuelFile = '/OECD_NEA_3by3_MESH_FUEL_TEMP_MC.out'; % 핵연료 온도 구배
TbulkFile = '/OECD_NEA_3by3_MESH_BULK_TEMP_MC.out'; % 냉각재 온도 구배

HBMpower_std = '/OECD_NEA_3by3_HBM_MESH_POWER_MC_STD.out';     % MC_AVG 파일
HBMTfuel_std = '/OECD_NEA_3by3_HBM_MESH_FUEL_TEMP_MC_STD.out'; % 핵연료 온도 구배
HBMTbulk_std = '/OECD_NEA_3by3_HBM_MESH_BULK_TEMP_MC_STD.out'; % 냉각재 온도 구배

APPpower_std = '/OECD_NEA_3by3_MESH_POWER_MC_STD.out'; % MC_AVG 파일

Nx = 51; Ny = 51; Nz = 20;

% --- RE-DECLARE THE PATH FOR EACH DATA RESULT
powerFile = strcat(pathHead,powerFile);
TfuelFile = strcat(pathHead,TfuelFile);
TbulkFile = strcat(pathHead,TbulkFile);

HBMpower_std = strcat(pathHead,HBMpower_std);
HBMTfuel_std = strcat(pathHead,HBMTfuel_std);
HBMTbulk_std = strcat(pathHead,HBMTbulk_std);

APPpower_std = strcat(pathHead,APPpower_std);

%% ====================== LOAD FILES =========================
power_avg = load_mesh(powerFile, Nx, Ny, Nz);  
fuel_dist = load_mesh(TfuelFile, Nx, Ny, Nz);
bulk_dist = load_mesh(TbulkFile, Nx+1, Ny+1, Nz+1);

power_std_HBM = load_mesh(HBMpower_std, Nx, Ny, Nz);  
fuel_std_HBM = load_mesh(HBMTfuel_std, Nx, Ny, Nz);  
bulk_std_HBM = load_mesh(HBMTbulk_std, Nx+1, Ny+1, Nz+1); 

power_std_APP = load_mesh(APPpower_std, Nx, Ny, Nz);  

% %% ====================== PLOT THE DISTRIBUTION INFORMATION ======================
% 전역 보기 옵션(배경/폰트) — 필요 없으면 지워도 됨
set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesFontSize',12);
set(groot,'defaultColorbarFontSize',10);

% 중앙 슬라이스 인덱스
[iP,jP,kP] = central_ijk(size(power_avg));
[iF,jF,kF] = central_ijk(size(fuel_dist));
[iB,jB,kB] = central_ijk(size(bulk_dist));

% % ===================== (1) Averages (1x3) =====================
% fAvg = figure('Name','3D slices: Averages (central planes)','Position',[80 60 1800 560]);
% tAvg = tiledlayout(fAvg,1,3,"Padding","compact","TileSpacing","compact");
% 
% % Power mean
% ax = nexttile(tAvg); set(ax,'PositionConstraint','outerposition');
% Sliceplot3D_Volume(power_avg, iP, jP, kP, true);   % 0/NaN 투명
% title('Power (mean)');
% cb = colorbar(ax,'Location','eastoutside'); thinColorbar(cb,0.6);
% clim([0.0 2.2]);
% 
% % Fuel temperature mean
% ax = nexttile(tAvg); set(ax,'PositionConstraint','outerposition');
% Sliceplot3D_Volume(fuel_dist, iF, jF, kF, true);  % 0/NaN만 투명
% title('Fuel temperature (mean)');
% cb = colorbar(ax,'Location','eastoutside'); thinColorbar(cb,0.6);
% clim([600.0 1250.0]);
% 
% % Coolant temperature mean
% ax = nexttile(tAvg); set(ax,'PositionConstraint','outerposition');
% Sliceplot3D_Volume(bulk_dist, iB, jB, kB, false);
% title('Coolant temperature (mean)');
% cb = colorbar(ax,'Location','eastoutside'); thinColorbar(cb,0.6);
% clim([560.0 610.0]);
% 
% ===================== (2) Standard deviations (1x3) =====================
% fStd = figure('Name','3D slices: Standard deviations (central planes)','Position',[100 80 1800 560]);
% tStd = tiledlayout(fStd,1,3,"Padding","compact","TileSpacing","compact");
% 
% % Power std
% ax = nexttile(tStd); set(ax,'PositionConstraint','outerposition');
% Sliceplot3D_Volume(power_std_HBM, iP, jP, kP, true);   % 0/NaN 투명
% title('Power STD');
% cb = colorbar(ax,'Location','eastoutside'); thinColorbar(cb,0.6); 
% clim([0.0 0.025]);
% 
% % Fuel temperature std
% ax = nexttile(tStd); set(ax,'PositionConstraint','outerposition');
% Sliceplot3D_Volume(fuel_std_HBM, iF, jF, kF, true);
% title('Fuel temperature STD');
% cb = colorbar(ax,'Location','eastoutside'); thinColorbar(cb,0.6);
% clim([0.0 10.0]);
% 
% % Coolant temperature std
% ax = nexttile(tStd); set(ax,'PositionConstraint','outerposition');
% Sliceplot3D_Volume(bulk_std_HBM, iB, jB, kB, false);
% title('Coolant temperature STD');
% cb = colorbar(ax,'Location','eastoutside'); thinColorbar(cb,0.6);
% clim([0.0 0.20]);

%% ====================== METHOD 1: L2-norm uncertainty (no propagation) ======================
% E_abs = ||sigma||_2,  E_rel = ||sigma||_2 / ||mu||_2   (NaN 무시)

[L2M1.power.abs,  L2M1.power.rel,  L2M1.power.N]  = L2_uncert_m1(power_avg,  power_std_HBM,  'Power');
[L2M1.Tfuel.abs,  L2M1.Tfuel.rel,  L2M1.Tfuel.N]  = L2_uncert_m1(fuel_dist,  fuel_std_HBM,  'Fuel temperature');
[L2M1.Tbulk.abs,  L2M1.Tbulk.rel,  L2M1.Tbulk.N]  = L2_uncert_m1(bulk_dist,  bulk_std_HBM,  'Coolant temperature');
[L2M1.power_app.abs,  L2M1.power_app.rel,  L2M1.power_app.N]  = L2_uncert_m1(power_avg,  power_std_APP,  'Power [APP]');

% --- COMPARISON WITH THE REFERENCE DATA: HOW MUCH DEVIATION DOES IT HAVE
load('NEA_3by3_REF_DATA.mat')
diff_power = power_avg - power_avg_ref;
diff_Tfuel = fuel_dist - Tfuel_avg_ref;
diff_Tbulk = bulk_dist - Tbulk_avg_ref;

[L2M1.power_diff.abs,  L2M1.power_diff.rel,  L2M1.power_diff.N]  = L2_uncert_m1(power_avg_ref,  diff_power,  'Power Diff');
[L2M1.Tfuel_diff.abs,  L2M1.Tfuel_diff.rel,  L2M1.Tfuel_diff.N]  = L2_uncert_m1(Tfuel_avg_ref,  diff_Tfuel,  'Fuel Diff');
[L2M1.Tbulk_diff.abs,  L2M1.Tbulk_diff.rel,  L2M1.Tbulk_diff.N]  = L2_uncert_m1(Tbulk_avg_ref,  diff_Tbulk,  'Bulk Diff');

fprintf('\n=== Method 1 (plain L2: HBM) ===\n');
fprintf('Power       : abs=%.6e, rel=%.6e, N=%d\n', L2M1.power.abs,  L2M1.power.rel,  L2M1.power.N);
fprintf('Power [APP] : abs=%.6e, rel=%.6e, N=%d\n', L2M1.power_app.abs,  L2M1.power_app.rel,  L2M1.power_app.N);
fprintf('Fuel T      : abs=%.6e, rel=%.6e, N=%d\n', L2M1.Tfuel.abs,  L2M1.Tfuel.rel,  L2M1.Tfuel.N);
fprintf('Coolant     : abs=%.6e, rel=%.6e, N=%d\n', L2M1.Tbulk.abs,  L2M1.Tbulk.rel,  L2M1.Tbulk.N);
fprintf('Power [DIF] : abs=%.6e, rel=%.6e, N=%d\n', L2M1.power_diff.abs,  L2M1.power_diff.rel,  L2M1.power_diff.N);
fprintf('Fuel  [DIF] : abs=%.6e, rel=%.6e, N=%d\n', L2M1.Tfuel_diff.abs,  L2M1.Tfuel_diff.rel,  L2M1.Tfuel_diff.N);
fprintf('BULK  [DIF] : abs=%.6e, rel=%.6e, N=%d\n', L2M1.Tbulk_diff.abs,  L2M1.Tbulk_diff.rel,  L2M1.Tbulk_diff.N);

% %% 그림그리기 (VISUALIZE: diff_power, diff_Tfuel, diff_Tbulk)
% % 전역 보기 옵션(배경/폰트)
% set(groot,'defaultFigureColor','w');
% set(groot,'defaultAxesFontSize',12);
% set(groot,'defaultColorbarFontSize',10);
% 
% % 중앙 슬라이스 인덱스 (ref 사이즈 기준)
% [iP,jP,kP] = central_ijk(size(power_avg_ref));
% [iF,jF,kF] = central_ijk(size(Tfuel_avg_ref));
% [iB,jB,kB] = central_ijk(size(Tbulk_avg_ref));
% 
% % ===================== ABSOLUTE DIFFERENCES (1x3) =====================
% fStd = figure('Name','3D slices: ABSOLUTE DIFFERENCES (central planes)','Position',[100 80 1800 560]);
% tStd = tiledlayout(fStd,1,3,"Padding","compact","TileSpacing","compact");
% 
% % Power std
% ax = nexttile(tStd); set(ax,'PositionConstraint','outerposition');
% Sliceplot3D_Volume(abs(diff_power), iP, jP, kP, true);   % 0/NaN 투명
% title('Power [DIF]');
% cb = colorbar(ax,'Location','eastoutside'); thinColorbar(cb,0.6);
% clim([0.0 0.025]);
% 
% % Fuel temperature std
% ax = nexttile(tStd); set(ax,'PositionConstraint','outerposition');
% Sliceplot3D_Volume(abs(diff_Tfuel), iF, jF, kF, true);
% title('Fuel temperature [DIF]');
% cb = colorbar(ax,'Location','eastoutside'); thinColorbar(cb,0.6);
% clim([0.0 10.0]);
% 
% % Coolant temperature std
% ax = nexttile(tStd); set(ax,'PositionConstraint','outerposition');
% Sliceplot3D_Volume(abs(diff_Tbulk), iB, jB, kB, false);
% title('Coolant temperature [DIF]');
% cb = colorbar(ax,'Location','eastoutside'); thinColorbar(cb,0.6);
% clim([0.0 0.20]);

%% ====================== Y–Z STD SLICES @ FIXED X (HBM, VERTICAL) ======================
x_list = [8, 25, 42];
Nx_ref = size(power_std_HBM,1);

for ix = 1:numel(x_list)
    x0 = x_list(ix);
    x0 = max(1, min(Nx_ref, x0));   % safety clamp

    fYZ = figure('Name',sprintf('HBM STD Y–Z slices @ x = %d', x0), ...
                 'Position',[250 40 750 950]);   % 세로로 길쭉

    tYZ = tiledlayout(fYZ,3,1, ...
        "Padding","compact","TileSpacing","compact");

    % ================= Power STD =================
    Pstd_yz = squeeze(power_std_HBM(x0, :, :));
    Pstd_yz(~isfinite(Pstd_yz) | Pstd_yz==0) = NaN;

    ax = nexttile(tYZ);
    imagesc(Pstd_yz.');
    axis image; axis xy;
    xlabel('Y index'); ylabel('Z index');
    title(sprintf('Power STD (Y–Z @ x = %d)', x0),'FontSize',15);
    cb = colorbar; cb.Label.String='[W]';
    clim([0.0 0.025]);
    set(ax,'FontSize',13);

    % ================= Fuel temperature STD =================
    Fstd_yz = squeeze(fuel_std_HBM(x0, :, :));
    Fstd_yz(~isfinite(Pstd_yz)) = NaN;   % power-based mask (일관성)

    ax = nexttile(tYZ);
    imagesc(Fstd_yz.');
    axis image; axis xy;
    xlabel('Y index'); ylabel('Z index');
    title(sprintf('Fuel temperature STD (Y–Z @ x = %d)', x0),'FontSize',15);
    cb = colorbar; cb.Label.String='[K]';
    clim([0.0 10.0]);
    set(ax,'FontSize',13);

    % ================= Coolant temperature STD =================
    xB = min(x0, size(bulk_std_HBM,1));   % staggered grid
    Bstd_yz = squeeze(bulk_std_HBM(xB, :, :));

    ax = nexttile(tYZ);
    imagesc(Bstd_yz.');
    axis image; axis xy;
    xlabel('Y index'); ylabel('Z index');
    title(sprintf('Coolant temperature STD (Y–Z @ x = %d)', x0),'FontSize',15);
    cb = colorbar; cb.Label.String='[K]';
    clim([0.0 0.20]);
    set(ax,'FontSize',13);
end

%% ====================== FUNCTIONS ==========================
function A = load_mesh(path, Nx, Ny, Nz)
    % 17x17 행렬이 z=1..Nz 순서로 20개 이어진 텍스트를 (Nx,Ny,Nz)로 변환
    txt  = fileread(path);
    vals = sscanf(txt, '%f');              % 전부 실수로 읽음 (열벡터)
    nExp = Nx*Ny*Nz;
    assert(numel(vals)==nExp, 'Value count mismatch: %s', path);
    A = reshape(vals, [Ny, Nx, Nz]);       % 줄 단위가 y-방향이라고 가정
    A = permute(A, [2 1 3]);               % (Nx,Ny,Nz)
end

function [i0,j0,k0] = central_ijk(sz)
    i0 = ceil(sz(1)/2); j0 = ceil(sz(2)/2); k0 = ceil(sz(3)/2);
end

function Sliceplot3D_Volume(vol, i, j, k, maskZero)
    % 셀 중심 좌표(0.5..N-0.5) + 축 0..N (하단 잘림 방지)
    if nargin<5, maskZero=false; end
    [nx,ny,nz] = size(vol);
    xc = (0.5:1:nx-0.5); yc = (0.5:1:ny-0.5); zc = (0.5:1:nz-0.5);
    i = max(1,min(nx,i)); j = max(1,min(ny,j)); k = max(1,min(nz,k));

    h = slice(xc, yc, zc, vol, xc(i), yc(j), zc(k));
    shading flat; axis vis3d;
    xlim([0 nx]); ylim([0 ny]); zlim([0 nz]);
    set(gca,'XDir','normal','YDir','normal','ZDir','normal');

    % 0(또는 NaN) 셀만 투명 처리 (데이터는 유지)
    if all(isgraphics(h))
        hs = h(:).';
        for hh = hs
            c = get(hh,'CData');                        % 2D slice의 값들
            if maskZero
                alphaMask = isfinite(c) & (c ~= 0);     % 0은 투명
            else
                alphaMask = isfinite(c);                % NaN만 투명
            end
            set(hh,'AlphaData', double(alphaMask), ...
                   'AlphaDataMapping','none', ...       % 0/1 그대로 사용
                   'FaceAlpha','flat', ...
                   'EdgeColor','none');
        end
    end
end

function thinColorbar(cb, scale)
    if nargin<2, scale=0.6; end
    drawnow; p = cb.Position; p(3) = p(3)*scale; cb.Position = p;
end

function [E_abs, E_rel, Ncells] = L2_uncert_m1(mu3d, sigma3d, label)
    % mu3d, sigma3d : same size 3D (또는 3D에 준함), NaN은 무시
    m = isfinite(mu3d) & isfinite(sigma3d);
    mu = mu3d(m);
    sg = sigma3d(m);
    Ncells = numel(mu);

    A2 = sum(sg.^2, 'omitnan');     % ||sigma||_2^2
    D2 = sum(mu.^2, 'omitnan');     % ||mu||_2^2
    E_abs = sqrt(max(A2,0));
    denom = sqrt(max(D2, realmin));
    E_rel = E_abs / denom;

    if nargin >= 3 && ~isempty(label)
        fprintf('[%s] abs=%.6e, rel=%.6e, N=%d\n', label, E_abs, E_rel, Ncells);
    end
end

