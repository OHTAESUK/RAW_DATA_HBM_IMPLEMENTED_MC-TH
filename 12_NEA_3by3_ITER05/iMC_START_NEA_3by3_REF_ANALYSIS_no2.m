% ======================================== %
% iMC 배치 계산 결과별 가분산 데이터를 정리함
% ======================================== %
num_bat = 40;
Nx = 51;
Ny = 51;
Nz = 20;

%% READ BATCH WISE DATA & CALCULATE THE RELATIVE L2-NORM LIKE ERROR FOR EACH BATCH
vec_E_rel = zeros([num_bat,1]);

for i = 1:num_bat
    if(i < 10)
        head_path = strcat('TRIAL0',num2str(i));
    else
        head_path = strcat('TRIAL',num2str(i));
    end
    avgFile = 'OECD_NEA_3by3_MESH_POWER_MC_AVG.out';   % MC_AVG 파일
    stdFile = 'OECD_NEA_3by3_MESH_POWER_MC_STD.out';   % MC_STD 파일
    % READ THE DATA
    power_avg = load_mesh(strcat(head_path,'\',avgFile), Nx, Ny, Nz);
    power_std = load_mesh(strcat(head_path,'\',stdFile), Nx, Ny, Nz);
    % CALCULATE THE RELATIVE L2-NORM LIKE ERROR
    [err_abs,err_rel,tmp] = L2_uncert_m1(power_avg,power_std);
    fprintf('Batch %02d: rel = %.6f\n', i, err_rel);
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

function [E_abs, E_rel, Ncells] = L2_uncert_m1(mu3d, sigma3d)
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
end

