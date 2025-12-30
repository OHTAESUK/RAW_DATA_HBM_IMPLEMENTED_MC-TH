% ==================================== %
% iMC 배치 계산 결과를 토대로 기준해 계산 %
% ==================================== %
num_bat = 40;
Nx = 51;
Ny = 51;
Nz = 20;

%% READ BATCH WISE DATA 
mat_power = zeros([num_bat, Nx, Ny, Nz]);
mat_Tfuel = zeros([num_bat, Nx, Ny, Nz]);
mat_Tbulk = zeros([num_bat, Nx+1, Ny+1, Nz+1]);

for i = 1:num_bat
    if(i < 10)
        head_path = strcat('TRIAL0',num2str(i));
    else
        head_path = strcat('TRIAL',num2str(i));
    end
    avgFile = 'OECD_NEA_3by3_MESH_POWER_MC_AVG.out';   % MC_AVG 파일
    TfuelFile = 'OECD_NEA_3by3_MESH_FUEL_TEMP_MC.out'; % 핵연료 온도 구배
    TbulkFile = 'OECD_NEA_3by3_MESH_BULK_TEMP_MC.out'; % 냉각재 온도 구배
    % READ THE DATA
    power_avg = load_mesh(strcat(head_path,'\',avgFile), Nx, Ny, Nz);
    fuel_dist = load_mesh(strcat(head_path,'\',TfuelFile), Nx, Ny, Nz);
    bulk_dist = load_mesh(strcat(head_path,'\',TbulkFile), Nx+1, Ny+1, Nz+1);
    % APPEND THE DATA
    mat_power(i,:,:,:) = power_avg;
    mat_Tfuel(i,:,:,:) = fuel_dist;
    mat_Tbulk(i,:,:,:) = bulk_dist;
end

%% FROM THE BATCH WISE DATA, ACQUIRE THE TRUE MEAN & STDDEV
power_avg_ref = squeeze(mean(mat_power,1));
power_std_ref = squeeze(std (mat_power,1,1));

Tfuel_avg_ref = squeeze(mean(mat_Tfuel,1));
Tfuel_std_ref = squeeze(std (mat_Tfuel,1,1));

Tbulk_avg_ref = squeeze(mean(mat_Tbulk,1));
Tbulk_std_ref = squeeze(std (mat_Tbulk,1,1));

%% ====================== PLOT THE DISTRIBUTION INFORMATION ======================
% 전역 보기 옵션(배경/폰트)
set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesFontSize',12);
set(groot,'defaultColorbarFontSize',10);

% 중앙 슬라이스 인덱스 (ref 사이즈 기준)
[iP,jP,kP] = central_ijk(size(power_avg_ref));
[iF,jF,kF] = central_ijk(size(Tfuel_avg_ref));
[iB,jB,kB] = central_ijk(size(Tbulk_avg_ref));

% % ===================== (1) Averages (1x3) =====================
% fAvg = figure('Name','3D slices: Averages (central planes)','Position',[80 60 1800 560]);
% tAvg = tiledlayout(fAvg,1,3,"Padding","compact","TileSpacing","compact");
% 
% % Power mean
% ax = nexttile(tAvg); set(ax,'PositionConstraint','outerposition');
% Sliceplot3D_Volume(power_avg_ref, iP, jP, kP, true);   % 0/NaN 투명
% title('Power (mean)');
% cb = colorbar(ax,'Location','eastoutside'); thinColorbar(cb,0.6);
% clim([0.0 2.2]);
% 
% % Fuel temperature mean
% ax = nexttile(tAvg); set(ax,'PositionConstraint','outerposition');
% Sliceplot3D_Volume(Tfuel_avg_ref, iF, jF, kF, true);  % 0/NaN만 투명
% title('Fuel temperature (mean)');
% cb = colorbar(ax,'Location','eastoutside'); thinColorbar(cb,0.6);
% clim([600.0 1250.0]);
% 
% % Coolant temperature mean
% ax = nexttile(tAvg); set(ax,'PositionConstraint','outerposition');
% Sliceplot3D_Volume(Tbulk_avg_ref, iB, jB, kB, false);
% title('Coolant temperature (mean)');
% cb = colorbar(ax,'Location','eastoutside'); thinColorbar(cb,0.6);
% clim([560.0 610.0]);
% 
% % ===================== (2) Standard deviations (1x3) =====================
% fStd = figure('Name','3D slices: Standard deviations (central planes)','Position',[100 80 1800 560]);
% tStd = tiledlayout(fStd,1,3,"Padding","compact","TileSpacing","compact");
% 
% % Power std
% ax = nexttile(tStd); set(ax,'PositionConstraint','outerposition');
% Sliceplot3D_Volume(power_std_ref, iP, jP, kP, true);   % 0/NaN 투명
% title('Power STD');
% cb = colorbar(ax,'Location','eastoutside'); thinColorbar(cb,0.6);
% clim([0.0 0.025]);
% 
% % Fuel temperature std
% ax = nexttile(tStd); set(ax,'PositionConstraint','outerposition');
% Sliceplot3D_Volume(Tfuel_std_ref, iF, jF, kF, true);
% title('Fuel temperature STD');
% cb = colorbar(ax,'Location','eastoutside'); thinColorbar(cb,0.6);
% clim([0.0 10.0]);
% 
% % Coolant temperature std
% ax = nexttile(tStd); set(ax,'PositionConstraint','outerposition');
% Sliceplot3D_Volume(Tbulk_std_ref, iB, jB, kB, false);
% title('Coolant temperature STD');
% cb = colorbar(ax,'Location','eastoutside'); thinColorbar(cb,0.6);
% clim([0.0 0.20]);

%% ====================== PLOT THE 2D (CONDENSED) IMAGE ====================== 
% 2D 평균
P2D = mean(power_avg_ref, 3, 'omitnan');      % [Nx, Ny]
F2D = mean(Tfuel_avg_ref, 3, 'omitnan');      % [Nx, Ny]
B2D = mean(Tbulk_avg_ref, 3, 'omitnan');      % [(Nx+1), (Ny+1)]

P2D(P2D==0) = NaN;
F2D(F2D==0) = NaN;
B2D(B2D==0) = NaN;

% 마스크: 출력(파워) 없는 위치는 투명 처리 (Power/Fuel에만 적용)
maskP = isfinite(P2D) & (P2D > 0);

% 보기 옵션
set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesFontSize',12);
set(groot,'defaultColorbarFontSize',10);

f2d = figure('Name','2D z-mean (masked by power)','Position',[100 80 1800 560]);
t2d = tiledlayout(f2d,1,3,"Padding","compact","TileSpacing","compact");

% --- Power (z-mean) ---
ax = nexttile(t2d);
h  = imagesc(P2D);
set(h,'AlphaData', double(maskP), 'AlphaDataMapping','none');   % 출력 0인 곳 투명
axis image; axis xy; set(ax,'Color','w');
title('Power (z-mean)');
cb = colorbar(ax,'Location','eastoutside'); cb.Label.String='[W]';
clim([0.6 1.4]);

% --- Fuel temperature (z-mean) ---
ax = nexttile(t2d);
h  = imagesc(F2D);
set(h,'AlphaData', double(maskP), 'AlphaDataMapping','none');   % 출력 0인 곳 투명
axis image; axis xy; set(ax,'Color','w');
title('Fuel temperature (z-mean)');
cb = colorbar(ax,'Location','eastoutside'); cb.Label.String='[K]';
clim([700 1050]);

% --- Coolant temperature (z-mean) ---
ax = nexttile(t2d);
h  = imagesc(B2D);
% 냉각재는 전체 영역을 보이게(마스크 적용 안 함)
axis image; axis xy; set(ax,'Color','w');
title('Coolant temperature (z-mean)');
cb = colorbar(ax,'Location','eastoutside'); cb.Label.String='[K]';
clim([570 590]);

%% ====================== Y–Z SLICES @ FIXED X (VERTICAL 3x1) ======================
x_list = [8, 25, 42];   % S1 / S2 / S3
Nx_ref = size(power_avg_ref,1);

for ix = 1:numel(x_list)
    x0 = x_list(ix);
    x0 = max(1, min(Nx_ref, x0));   % safety clamp

    fYZ = figure('Name',sprintf('Y–Z slices @ x=%d',x0), ...
                 'Position',[250 40 750 950]);

    tYZ = tiledlayout(fYZ,3,1, ...
        "Padding","compact","TileSpacing","compact");

    % ================= Power =================
    P_yz = squeeze(power_avg_ref(x0, :, :));   % (Ny, Nz)
    P_yz(P_yz==0) = NaN;

    ax = nexttile(tYZ);
    imagesc(P_yz.');
    axis image; axis xy;
    xlabel('Y index'); ylabel('Z index');
    title(sprintf('Power (Y–Z @ x=%d)', x0),'FontSize',16);
    cb = colorbar; cb.Label.String='[W]';
    clim([0.0 2.2]);
    set(ax,'FontSize',14);

    % ================= Fuel temperature =================
    F_yz = squeeze(Tfuel_avg_ref(x0, :, :));
    F_yz(isnan(P_yz)) = NaN;

    ax = nexttile(tYZ);
    imagesc(F_yz.');
    axis image; axis xy;
    xlabel('Y index'); ylabel('Z index');
    title(sprintf('Fuel temperature (Y–Z @ x=%d)', x0),'FontSize',16);
    cb = colorbar; cb.Label.String='[K]';
    clim([600 1250]);
    set(ax,'FontSize',14);

    % ================= Coolant temperature =================
    xB = min(x0, size(Tbulk_avg_ref,1));   % staggered grid
    B_yz = squeeze(Tbulk_avg_ref(xB, :, :));

    ax = nexttile(tYZ);
    imagesc(B_yz.');
    axis image; axis xy;
    xlabel('Y index'); ylabel('Z index');
    title(sprintf('Coolant temperature (Y–Z @ x=%d)', x0),'FontSize',16);
    cb = colorbar; cb.Label.String='[K]';
    clim([560 610]);
    set(ax,'FontSize',14);
end

%% ====================== Y–Z STD SLICES @ FIXED X (VERTICAL 3x1) ======================
x_list = [8, 25, 42];   % S1 / S2 / S3
Nx_ref = size(power_std_ref,1);

for ix = 1:numel(x_list)
    x0 = x_list(ix);
    x0 = max(1, min(Nx_ref, x0));   % safety clamp

    fYZ = figure('Name',sprintf('Y–Z STD slices @ x=%d',x0), ...
                 'Position',[250 40 750 950]);

    tYZ = tiledlayout(fYZ,3,1, ...
        "Padding","compact","TileSpacing","compact");

    % ================= Power STD =================
    Pstd_yz = squeeze(power_std_ref(x0, :, :));
    Pstd_yz(~isfinite(Pstd_yz)) = NaN;

    ax = nexttile(tYZ);
    imagesc(Pstd_yz.');
    axis image; axis xy;
    xlabel('Y index'); ylabel('Z index');
    title(sprintf('Power STD (Y–Z @ x=%d)', x0),'FontSize',16);
    cb = colorbar; cb.Label.String='[W]';
    clim([0.0 0.025]);   % 🔥 이미 네가 위에서 쓴 값
    set(ax,'FontSize',14);

    % ================= Fuel temperature STD =================
    Fstd_yz = squeeze(Tfuel_std_ref(x0, :, :));
    Fstd_yz(~isfinite(Pstd_yz)) = NaN;   % power-based mask (선택)

    ax = nexttile(tYZ);
    imagesc(Fstd_yz.');
    axis image; axis xy;
    xlabel('Y index'); ylabel('Z index');
    title(sprintf('Fuel temperature STD (Y–Z @ x=%d)', x0),'FontSize',16);
    cb = colorbar; cb.Label.String='[K]';
    clim([0.0 10.0]);
    set(ax,'FontSize',14);

    % ================= Coolant temperature STD =================
    xB = min(x0, size(Tbulk_std_ref,1));   % staggered grid
    Bstd_yz = squeeze(Tbulk_std_ref(xB, :, :));

    ax = nexttile(tYZ);
    imagesc(Bstd_yz.');
    axis image; axis xy;
    xlabel('Y index'); ylabel('Z index');
    title(sprintf('Coolant temperature STD (Y–Z @ x=%d)', x0),'FontSize',16);
    cb = colorbar; cb.Label.String='[K]';
    clim([0.0 0.20]);
    set(ax,'FontSize',14);
end


%% ====================== METHOD 1: L2-norm uncertainty (no propagation) ======================
% E_abs = ||sigma||_2,  E_rel = ||sigma||_2 / ||mu||_2   (NaN 무시)

[L2M1.power.abs,  L2M1.power.rel,  L2M1.power.N]  = L2_uncert_m1(power_avg_ref,  power_std_ref,  'Power');
[L2M1.Tfuel.abs,  L2M1.Tfuel.rel,  L2M1.Tfuel.N]  = L2_uncert_m1(Tfuel_avg_ref,  Tfuel_std_ref,  'Fuel temperature');
[L2M1.Tbulk.abs,  L2M1.Tbulk.rel,  L2M1.Tbulk.N]  = L2_uncert_m1(Tbulk_avg_ref,  Tbulk_std_ref,  'Coolant temperature');

fprintf('\n=== Method 1 (plain L2) ===\n');
fprintf('Power   : abs=%.6e, rel=%.6e, N=%d\n', L2M1.power.abs,  L2M1.power.rel,  L2M1.power.N);
fprintf('Fuel T  : abs=%.6e, rel=%.6e, N=%d\n', L2M1.Tfuel.abs,  L2M1.Tfuel.rel,  L2M1.Tfuel.N);
fprintf('Coolant : abs=%.6e, rel=%.6e, N=%d\n', L2M1.Tbulk.abs,  L2M1.Tbulk.rel,  L2M1.Tbulk.N);

% SAVINGS (MATLAB 데이터 형식으로 저장)
save('NEA_3by3_REF_DATA','power_avg_ref','Tfuel_avg_ref','Tbulk_avg_ref',...
                         'power_std_ref','Tfuel_std_ref','Tbulk_std_ref')

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