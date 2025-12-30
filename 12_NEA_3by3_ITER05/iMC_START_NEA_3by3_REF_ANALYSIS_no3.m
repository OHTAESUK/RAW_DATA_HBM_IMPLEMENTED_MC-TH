% =================================================================== %
% 특정 배치 계산 결과 읽어서 가분산 분포 그림 그려주는 스크립트         %
% =================================================================== %
num_bat = 15;
Nx = 51; Ny = 51; Nz = 20;

i = num_bat;
if(i < 10)
    head_path = strcat('TRIAL0',num2str(i));
else
    head_path = strcat('TRIAL',num2str(i));
end
avgFile = strcat(head_path,'\OECD_NEA_3by3_MESH_POWER_MC_AVG.out');   % MC_AVG 파일
stdFile = strcat(head_path,'\OECD_NEA_3by3_MESH_POWER_MC_STD.out');   % MC_STD 파일
TfuelFile = strcat(head_path,'\OECD_NEA_3by3_MESH_FUEL_TEMP_MC.out'); % 핵연료 온도 구배
TbulkFile = strcat(head_path,'\OECD_NEA_3by3_MESH_BULK_TEMP_MC.out'); % 냉각재 온도 구배

%% ====================== LOAD FILES =========================
power_avg = load_mesh(avgFile, Nx, Ny, Nz);   % 평균 파워 [W]
power_std = load_mesh(stdFile, Nx, Ny, Nz);   % 표준편차 [W]
fuel_dist = load_mesh(TfuelFile, Nx, Ny, Nz);
bulk_dist = load_mesh(TbulkFile, Nx+1, Ny+1, Nz+1);

%% ====================== PLOT THE DISTRIBUTION INFORMATION ======================
% 전역 보기 옵션(배경/폰트) — 필요 없으면 지워도 됨
set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesFontSize',12);
set(groot,'defaultColorbarFontSize',10);

% 중앙 슬라이스 인덱스
[iP,jP,kP] = central_ijk(size(power_avg));
[iF,jF,kF] = central_ijk(size(fuel_dist));
[iB,jB,kB] = central_ijk(size(bulk_dist));

% % ---- 3D 슬라이스 (1x3) : Power / Fuel T / Coolant T ----
% f3d = figure('Name','3D slices (central planes)','Position',[80 60 1800 560]);
% t3d = tiledlayout(f3d,1,3,"Padding","compact","TileSpacing","compact");
% 
% % Power
% ax = nexttile(t3d); set(ax,'PositionConstraint','outerposition');
% Sliceplot3D_Volume(power_avg, iP, jP, kP, true);  % <- true: NaN 투명 처리
% title('Power (central slices)');
% cb = colorbar(ax,'Location','eastoutside'); thinColorbar(cb,0.6);
% clim([0.0 2.2]);
% % cb.Label.String='W';
% % (컬러맵 지정 안 함: 디폴트 사용)
% 
% % Fuel temperature
% ax = nexttile(t3d); set(ax,'PositionConstraint','outerposition');
% Sliceplot3D_Volume(fuel_dist, iF, jF, kF, true); % NaN 투명 필요 없으면 false
% title('Fuel temperature (central slices)');
% cb = colorbar(ax,'Location','eastoutside'); thinColorbar(cb,0.6);
% clim([600.0 1250.0]);
% 
% % Coolant temperature
% ax = nexttile(t3d); set(ax,'PositionConstraint','outerposition');
% Sliceplot3D_Volume(bulk_dist, iB, jB, kB, false);
% title('Coolant temperature (central slices)');
% cb = colorbar(ax,'Location','eastoutside'); thinColorbar(cb,0.6);
% clim([560.0 610.0]);
% 
% % ---- 3D slice: Power STD (power_avg와 동일 스타일) ----
% fStd = figure('Name','3D slice: Power STD (central planes)','Position',[80 60 700 560]);
% ax = axes(fStd); set(ax,'PositionConstraint','outerposition');
% 
% % 중앙 슬라이스는 power_avg에서 구한 (iP,jP,kP) 그대로 사용
% Sliceplot3D_Volume(power_std, iP, jP, kP, true);   % true: 값=0 또는 NaN 투명 처리
% title('Power STD (central slices)');
% cb = colorbar(ax,'Location','eastoutside'); thinColorbar(cb,0.6); clim([0.0 0.025]);
% % 컬러맵은 지정하지 않음(디폴트 사용), 단위 라벨도 생략(원하면 cb.Label.String='W'; 추가)

%% ====================== Power STD: Y–Z SLICES @ x = 8, 25, 42 (ONE FIGURE) ======================
x_list = [8, 25, 42];
Nx_ref = size(power_std,1);

fYZ = figure('Name','Power STD Y–Z slices (x = 8, 25, 42)', ...
             'Position',[250 40 750 950]);   % 세로로 길게

tYZ = tiledlayout(fYZ,3,1, ...
    "Padding","compact","TileSpacing","compact");

for ix = 1:numel(x_list)
    x0 = x_list(ix);
    x0 = max(1, min(Nx_ref, x0));   % safety clamp

    % ---- Power STD Y–Z slice ----
    Pstd_yz = squeeze(power_std(x0, :, :));   % (Ny, Nz)
    Pstd_yz(~isfinite(Pstd_yz) | Pstd_yz==0) = NaN;

    ax = nexttile(tYZ);
    imagesc(Pstd_yz.');
    axis image; axis xy;
    xlabel('Y index'); ylabel('Z index');
    title(sprintf('Power STD (Y–Z @ x = %d)', x0),'FontSize',15);

    cb = colorbar;
    cb.Label.String = '[W]';
    clim([0.0 0.025]);   % 위에서 사용한 STD 범위와 동일
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