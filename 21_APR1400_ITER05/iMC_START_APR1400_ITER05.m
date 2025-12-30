%% ===================== READ DATA =====================
iter_max = 5; 
Nx = 136; Ny = 136; Nz = 38;

filename = './ITER05_TRIAL01/log.out';
fid = fopen(filename); lines = textscan(fid, '%s', 'Delimiter', '\n'); lines = lines{1}; fclose(fid);

power    = zeros(iter_max,   Nx,   Ny,   Nz);      % 선출력 [W/m]가 계산됨
fueltemp = zeros(iter_max+1, Nx,   Ny,   Nz);      % temp: ITER+1 매핑
bulktemp = zeros(iter_max+1, Nx+1, Ny+1, Nz+1);    % bulk: 경계격자라고 가정

lastTempIter = NaN;  % ZVAL 통해 최신 temp iter(+1) 추적
for i = 1:length(lines)
    line = strtrim(lines{i});

    if startsWith(line, 'POWERZ')
        tokens = sscanf(line, 'POWERZ %d %d');   % z, iter
        if numel(tokens) >= 2
            z = tokens(1); iter = tokens(2);
            power(iter,:,:,z) = Data_Extract(lines,i+1,Nx,Ny);
        end

    elseif startsWith(line, 'ZVAL')
        tokens = sscanf(line, 'ZVAL %d %d');     % z, iterBase(0-based)
        if numel(tokens) >= 2
            z = tokens(1); iter = tokens(2)+1;   % temp는 +1
            fueltemp(iter,:,:,z) = Data_Extract(lines,i+2,Nx,Ny);
            lastTempIter = iter;                 % BULKTEMP iter 추적용
        end

    elseif startsWith(line, 'BULKTEMP')
        % 허용 포맷: "BULKTEMP z" 또는 "BULKTEMP z iterBase"
        ints = sscanf(line, 'BULKTEMP %d %d');
        if numel(ints) >= 1
            z = ints(1);
            if numel(ints) >= 2
                iter = ints(2)+1;               % 로그가 0-based면 +1
            else
                % iter가 없으면 직전에 읽은 temp iter를 따름
                if isnan(lastTempIter)
                    warning('BULKTEMP iter 불명확: 직전 temp iter가 없어 스킵합니다 (line %d).', i);
                    continue;
                end
                iter = lastTempIter;
            end
            bulktemp(iter,:,:,z) = Data_Extract(lines,i+1,Nx+1,Ny+1);
        end
    end
end

% 선출력 -> 그냥 출력으로 변경
power = power * 365.76/100.0 / 20.;

% % 마스킹
% fueltemp(fueltemp<10) = NaN;
% power(power==0)       = NaN;

% % 각 이터레이션 출력분포 L2 norm의 불확도
E_n_rel1 = 2.993054e-02;
E_n_rel2 = 4.880705e-03;
E_n_rel3 = 1.706763e-05;

%% ===================== USER CONTROLS =====================
ITER = 5;                                    % 보기용: power=ITER, temp=ITER+1
dx = 1; dy = 1; dz = 1;                      % 격자 간격 (등간격이면 1,1,1)
totalPower_W = 350e6;                        % 총 출력 [W] (균일 baseline용)

%% ===================== GLOBAL PLOT DEFAULTS =====================
set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesFontSize',12);
set(groot,'defaultColorbarFontSize',10);

% %% ===================== (1) 초기 TH & ΔP(2-1) VIEW =====================
% ShowInitialTH(fueltemp, bulktemp);           % temp(iter=1)
% ShowDeltaPower12(power);                     % power(2)-power(1)

%% ===================== (2) 3D SLICES AT ITER =====================
% ShowIter3D(power, fueltemp, bulktemp, ITER, [], [], []);  % [] -> 중앙 자동
ShowIter2D(power, fueltemp, bulktemp, ITER);              % 새 2D (1x3, z-avg)

%% ===================== (3) CONVERGENCE (iter=1 포함) =====================
ConvergencePlot(power, fueltemp, bulktemp, [dx dy dz], totalPower_W, E_n_rel1,E_n_rel2,E_n_rel3);
% ConvergencePlot(power, fueltemp, bulktemp, [dx dy dz], totalPower_W);

%% ===================== FUNCTIONS =====================
function data = Data_Extract(raw, idx, x, y)
    data = zeros(x, y);
    col = 1; j = 1;
    while col <= x && (idx + j - 1) <= length(raw)
        nums = sscanf(strtrim(raw{idx + j - 1}), '%f')';
        if numel(nums) == y
            data(col, :) = nums; col = col + 1;
        end
        j = j + 1;
    end
end

function ShowInitialTH(fueltemp, bulktemp)
    figure('Name','Initial TH (temp iter=1)','Position',[100 100 1400 520]);
    tiledlayout(1,2,"Padding","compact","TileSpacing","compact");

    % fueltemp(iter=1)
    nexttile;
    volF = squeeze(fueltemp(1,:,:,:));
    [iF,jF,kF] = central_ijk(size(volF));
    % <<< CHANGED >>> 0/NaN 투명화: 연료온도는 출력(0) 구간 생략
    Sliceplot3D_Volume(volF, iF, jF, kF, true);
    title('Fuel temperature (iteration 1)'); cb = colorbar; thinColorbar(cb,0.6);

    % bulktemp(iter=1)
    nexttile;
    volB = squeeze(bulktemp(1,:,:,:));
    [iB,jB,kB] = central_ijk(size(volB));
    % <<< CHANGED >>> 벌크(냉각재)는 NaN만 투명(0도 보이게 하려면 false)
    Sliceplot3D_Volume(volB, iB, jB, kB, false);
    title('Coolant temperature (iteration 1)'); cb = colorbar; thinColorbar(cb,0.6);
end

function ShowDeltaPower12(power)
    if size(power,1) < 2, return; end
    vol1 = squeeze(power(1,:,:,:));
    vol2 = squeeze(power(2,:,:,:));
    dP   = vol2 - vol1;

    figure('Name','ΔPower (iter2 - iter1)','Position',[100 100 900 650]);
    [i0,j0,k0] = central_ijk(size(dP));
    % <<< CHANGED >>> ΔP에서도 0/NaN 투명화
    Sliceplot3D_Volume(dP, i0, j0, k0, true);
    title('Power difference (iteration 2 minus 1)'); cb = colorbar; thinColorbar(cb,0.6);
end

function ShowIter3D(power, fueltemp, bulktemp, ITER, ijk_power, ijk_fuel, ijk_bulk)
    f = figure('Name',sprintf('3D fields at coupling iteration %d',ITER),...
               'Position',[100 100 1800 560]);
    t = tiledlayout(f,1,3,"Padding","compact","TileSpacing","compact");

    % ---- Power (ITER) ----
    ax = nexttile(t); set(ax,'PositionConstraint','outerposition');
    volP = squeeze(power(ITER,:,:,:));
    cnt = sum(~isnan(volP(:)) & volP(:) ~= 0);
    volP = volP/mean(volP,'all','omitnan')*numel(volP)/cnt;
    if isempty(ijk_power), [i,j,k] = central_ijk(size(volP)); else, ijk=clamp_ijk(ijk_power,size(volP)); [i,j,k]=deal(ijk(1),ijk(2),ijk(3)); end
    % <<< CHANGED >>> Power: 0/NaN 투명
    Sliceplot3D_Volume(volP, i, j, k, true);
    title(sprintf('Power (iteration %d)',ITER)); cb = colorbar; thinColorbar(cb,0.6);
    % clim([0.0 2.2]);

    % ---- Fuel Temp (ITER+1) ----
    ax = nexttile(t); set(ax,'PositionConstraint','outerposition');
    volF = squeeze(fueltemp(ITER+1,:,:,:));
    if isempty(ijk_fuel), [i,j,k] = central_ijk(size(volF)); else, ijk=clamp_ijk(ijk_fuel,size(volF)); [i,j,k]=deal(ijk(1),ijk(2),ijk(3)); end
    % <<< CHANGED >>> Fuel: 0/NaN 투명
    Sliceplot3D_Volume(volF, i, j, k, true);
    title(sprintf('Fuel temperature (iteration %d)',ITER)); cb = colorbar; thinColorbar(cb,0.6);
    % clim([600.0 1250.0]);

    % ---- Coolant Temp (ITER+1) ----
    ax = nexttile(t); set(ax,'PositionConstraint','outerposition');
    volB = squeeze(bulktemp(ITER+1,:,:,:));
    if isempty(ijk_bulk), [i,j,k] = central_ijk(size(volB)); else, ijk=clamp_ijk(ijk_bulk,size(volB)); [i,j,k]=deal(ijk(1),ijk(2),ijk(3)); end
    % <<< CHANGED >>> Coolant: NaN만 투명(0은 표시)
    Sliceplot3D_Volume(volB, i, j, k, false);
    title(sprintf('Coolant temperature (iteration %d)',ITER)); cb = colorbar; thinColorbar(cb,0.6);
    % clim([560.0 610.0]);
end

function ShowIter2D(power, fueltemp, bulktemp, ITER)
    f = figure('Name',sprintf('2D (z-avg) at coupling iteration %d',ITER),...
               'Position',[100 100 1800 560]);
    t = tiledlayout(f,1,3,"Padding","compact","TileSpacing","compact");

    volP = squeeze(power(ITER,:,:,:));        % (Nx,Ny,Nz)
    volF = squeeze(fueltemp(ITER+1,:,:,:));   % (Nx,Ny,Nz)
    volB = squeeze(bulktemp(ITER+1,:,:,:));   % (Nx+1,Ny+1,Nz+1)

    mapP = mean(volP, 3, 'omitnan');      % [Nx, Ny]
    mapF = mean(volF, 3, 'omitnan');      % [Nx, Ny]
    mapB = mean(volB, 3, 'omitnan');      % [(Nx+1), (Ny+1)]

    mapP = mapP / mean(mapP,'all','omitnan');

    % <<< CHANGED >>> 0을 NaN으로 돌려 투명화 + 마스크 생성
    mapP(mapP==0) = NaN;
    mapF(mapF==0) = NaN;
    mapB(mapB==0) = NaN;

    maskP = isfinite(mapP) & (mapP > 0);
    maskF = isfinite(mapF) & (mapF > 0);
    % 냉각재는 전체 보이게 하고 싶으면 maskB 생략. (NaN만 투명)
    maskB = isfinite(mapB);  % 필요 시 (mapB>0)로 변경

    % ----- Power <z> -----
    ax = nexttile(t); set(ax,'PositionConstraint','outerposition');
    h = imagesc(mapP); axis image; axis xy;
    % <<< CHANGED >>> 2D에서도 AlphaData로 투명 처리
    set(h,'AlphaData', double(maskF), 'AlphaDataMapping','none');
    set(ax,'Color','w');
    title('Power  ⟨·⟩_z'); cb = colorbar(ax,'Location','eastoutside'); thinColorbar(cb,0.6);
    % clim([0.6 1.4]);

    % ----- Fuel T <z> -----
    ax = nexttile(t); set(ax,'PositionConstraint','outerposition');
    h = imagesc(mapF); axis image; axis xy;
    % <<< CHANGED >>>
    set(h,'AlphaData', double(maskF), 'AlphaDataMapping','none');
    set(ax,'Color','w');
    title('Fuel T ⟨·⟩_z'); cb = colorbar(ax,'Location','eastoutside'); thinColorbar(cb,0.6);
    % clim([700 1050]);

    % ----- Coolant T <z> -----
    ax = nexttile(t); set(ax,'PositionConstraint','outerposition');
    h = imagesc(mapB); axis image; axis xy;
    % <<< CHANGED >>> 냉각재: NaN만 투명(0도 보이게 하려면 maskB를 isfinite 로)
    set(h,'AlphaData', double(maskB), 'AlphaDataMapping','none');
    set(ax,'Color','w');
    title('Coolant T ⟨·⟩_z'); cb = colorbar(ax,'Location','eastoutside'); thinColorbar(cb,0.6);
    % clim([570 590]);
end

function ConvergencePlot(power, fueltemp, bulktemp, dxyz, totalPower_W, relNoiseLine1, relNoiseLine2, relNoiseLine3)
% ConvergencePlot
%   - Power:   ΔP(1)=P1-P0(균일 baseline), ΔP(n)=Pn-P(n-1), n>=2
%   - Temp:    ΔX(n)=X(n+1)-X(n), n>=1  (fuel/coolant)
%   - 왼쪽: 절대 L2, 오른쪽: 상대 L2 (y축 로그)
%   - relNoiseLine (옵션): E_n_rel을 넘기면 오른쪽 플롯에 2×E_n_rel 기준선을 표시

    if nargin<4 || isempty(dxyz), dxyz=[1 1 1]; end
    if nargin<5, totalPower_W = []; end
    if nargin<6
        relNoiseLine1 = []; 
        relNoiseLine2 = []; 
        relNoiseLine3 = [];  
    end

    w = prod(dxyz);                       % 등체적 가정 (비등간격이면 가중 적용부 확장 필요)

    itP = size(power,1);                  % power iteration count (e.g., 10)
    itT = size(fueltemp,1);               % temp iteration count   (e.g., 11)

    % ---------- POWER: ΔP ----------
    pow_abs = nan(itP,1);
    pow_rel = nan(itP,1);

    if itP>=1
        volP1 = squeeze(power(1,:,:,:));
        if ~isempty(totalPower_W) && any(isfinite(volP1(:)))
            % P0: 균일 baseline (총출력 totalPower_W 균등분배; NaN 제외)
            P0 = make_uniform_baseline_like(volP1, totalPower_W, []);
            [pow_abs(1), pow_rel(1)] = l2_pair(volP1, P0, w);  % n=1: P1 - P0
        else
            % totalPower_W 미지정 시 첫 점은 NaN 유지
            % warning('totalPower_W 미지정: ΔP(1)=P1-P0는 생략됩니다.');
        end
    end
    for n = 2:itP
        [pow_abs(n), pow_rel(n)] = l2_pair( squeeze(power(n,:,:,:)), squeeze(power(n-1,:,:,:)), w );
    end

    % ---------- FUEL / COOLANT: ΔX ----------
    fuel_abs = nan(max(itT-1,0),1);
    fuel_rel = nan(max(itT-1,0),1);
    bulk_abs = nan(max(itT-1,0),1);
    bulk_rel = nan(max(itT-1,0),1);

    for n = 1:itT-1
        [fuel_abs(n), fuel_rel(n)] = l2_pair( squeeze(fueltemp(n+1,:,:,:)), squeeze(fueltemp(n,:,:,:)), w );
        [bulk_abs(n), bulk_rel(n)] = l2_pair( squeeze(bulktemp(n+1,:,:,:)), squeeze(bulktemp(n,:,:,:)), w );
    end

    % ---------- PLOT ----------
    figure('Name','Convergence (L2 norms)','Position',[100 100 1200 450]);
    tiledlayout(1,2,"Padding","compact","TileSpacing","compact");

    % (Left) Absolute L2
    nexttile; hold on;
    semilogy(1:itP,   pow_abs,  '-o','DisplayName','Power  ||P_n - P_{n-1}||_2');
    semilogy(1:numel(fuel_abs), fuel_abs, '-s','DisplayName','Fuel   ||F_{n+1} - F_n||_2');
    semilogy(1:numel(bulk_abs), bulk_abs, '-^','DisplayName','Coolant||C_{n+1} - C_n||_2');
    set(gca,'YScale','log');
    grid on; legend('Location','best');
    xlabel('Iteration index n'); ylabel('Absolute L2 norm');

    % (Right) Relative L2
    nexttile; hold on;
    semilogy(1:itP,   pow_rel,  '-o','DisplayName','Power (relative)');
    semilogy(1:numel(fuel_rel), fuel_rel, '-s','DisplayName','Fuel (relative)');
    semilogy(1:numel(bulk_rel), bulk_rel, '-^','DisplayName','Coolant (relative)');
    set(gca,'YScale','log'); grid on;

    ax = gca;                           % 혹은 해당 axes 핸들
    CO = ax.ColorOrder;                 % MATLAB 기본 색상표 (Nx3)

    % 원하는 순서: 하늘색(6) -> 주황(2) -> 노랑(3)
    ord = [6 2 3];
    ord = ord(ord <= size(CO,1));       % 방어: 색상표 길이 확인
    cols = CO(ord, :);

    % 예: power / fuel / coolant 기준선 (2×σ 라인)
    if exist('relNoiseLine1','var') && ~isempty(relNoiseLine1) && isfinite(relNoiseLine1) && relNoiseLine1>0
        yline(ax, 2*relNoiseLine1, '--', 'Color', cols(1,:), 'LineWidth', 1.5, ...
              'DisplayName','2*\sigma_{n=5}^{HBM}(power)');
    end
    if exist('relNoiseLine2','var') && ~isempty(relNoiseLine2) && isfinite(relNoiseLine2) && relNoiseLine2>0
        yline(ax, 2*relNoiseLine2, '--', 'Color', cols(2,:), 'LineWidth', 1.5, ...
              'DisplayName','2*\sigma_{n=5}^{HBM}(fuel)');
    end
    if exist('relNoiseLine3','var') && ~isempty(relNoiseLine3) && isfinite(relNoiseLine3) && relNoiseLine3>0
        yline(ax, 2*relNoiseLine3, '--', 'Color', cols(3,:), 'LineWidth', 1.5, ...
              'DisplayName','2*\sigma_{n=5}^{HBM}(coolant)');
    end
    legend(ax,'Location','best');
    legend('Location','best');
    xlabel('Iteration index n'); ylabel('Relative L2 norm');

    % ---------- CONSOLE SUMMARY ----------
    il = @(v) find(isfinite(v),1,'last');
    ip = il(pow_abs);  if ~isempty(ip), fprintf('Power   last: abs=%.3e, rel=%.3e (n=%d)\n',  pow_abs(ip),  pow_rel(ip),  ip); end
    ifl = il(fuel_abs); if ~isempty(ifl), fprintf('Fuel    last: abs=%.3e, rel=%.3e (n=%d)\n', fuel_abs(ifl), fuel_rel(ifl), ifl); end
    ibl = il(bulk_abs); if ~isempty(ibl), fprintf('Coolant last: abs=%.3e, rel=%.3e (n=%d)\n',  bulk_abs(ibl), bulk_rel(ibl), ibl); end
end

function P0 = make_uniform_baseline_like(volLike, totalPower_W, V)
    % volLike: power(1,:,:,:)와 동일 크기 (NaN=비활성 셀)
    % totalPower_W: 총 출력 [W]
    % V: (옵션) 셀 체적 3D 배열; []면 등체적
    P0 = nan(size(volLike));
    mask = isfinite(volLike);
    if ~any(mask(:)), warning('baseline: 유효 셀이 없습니다.'); return; end
    if nargin>=3 && ~isempty(V)
        W = V; W(~mask)=0; S = sum(W(:));
        if S==0, error('체적 가중 합이 0입니다.'); end
        P0(mask) = totalPower_W * W(mask) / S;
    else
        n = nnz(mask);
        P0(mask) = totalPower_W / n;
    end
end

function [absL2, relL2] = l2_pair(B, A, weight)
    % L2 of (B-A), relative = ||B-A|| / ||B||  (NaN 무시)
    D = B - A; m = isfinite(D) & isfinite(B);
    if ~any(m(:)), absL2=NaN; relL2=NaN; return; end
    d2 = sum(D(m).^2)*weight; b2 = sum(B(m).^2)*weight + eps;
    absL2 = sqrt(d2); relL2 = sqrt(d2)/sqrt(b2);
end

function [i0,j0,k0] = central_ijk(sz)
    i0 = ceil(sz(1)/2); j0 = ceil(sz(2)/2); k0 = ceil(sz(3)/2);
end

function ijk = clamp_ijk(ijk, sz)
    ijk = max([1 1 1], min(sz, ijk));
end

% <<< CHANGED >>> Sliceplot3D_Volume: 0/NaN 투명화 지원 (maskZero=true면 0도 투명)
function Sliceplot3D_Volume(vol, i, j, k, maskZero)
    if nargin<5, maskZero = false; end
    [nx,ny,nz] = size(vol);
    % 좌표 벡터 방식으로 slice 호출
    h = slice(1:nx, 1:ny, 1:nz, vol, i, j, k);
    shading flat; axis vis3d;
    xlim([0 nx]); ylim([0 ny]); zlim([0 nz]);
    set(gca,'XDir','normal','YDir','normal','ZDir','normal');

    % ***** 핵심: 0/NaN 투명화 *****
    if all(isgraphics(h))
        hs = h(:).';
        for hh = hs
            c = get(hh,'CData');              % 2D slice 상의 값
            if maskZero
                alphaMask = isfinite(c) & (c ~= 0);  % 0도 투명
            else
                alphaMask = isfinite(c);             % NaN만 투명
            end
            set(hh,'AlphaData', double(alphaMask), ...
                   'AlphaDataMapping','none', ...
                   'FaceAlpha','flat', ...
                   'EdgeColor','none');
        end
    end
end

function thinColorbar(cb, scale)
    if nargin<2, scale=0.6; end
    drawnow; pos = cb.Position; pos(3) = pos(3)*scale; cb.Position = pos;
end
