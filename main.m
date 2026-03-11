% 2-CRT 序列 AAoI 主程式 (解析解 vs 蒙地卡羅)
clc; clear; close all;

% 加入 utils 資料夾路徑
addpath('utils');

%% 設定參數
K_list = [5, 6, 8];
num_runs = 2000;    % MC 隨機 τ 抽樣次數
N_period = 500;     % 每次 run 模擬的週期數

for K = K_list
    [Is, L, pK] = build_2CRT(K);
    fprintf('=== Building 2-CRT for K = %d (Length = %d) ===\n', K, L);
    
    %% 定義要比較的 set（對應 Table I 或自訂）
    if K == 6
        % Set 1: s0, s2, s3, s4, s5, s7 
        C_idx_1 = [1, 3, 4, 5, 6, 8];
        % Set 2: s1, s2, s3, s4, s5, s6
        C_idx_2 = [2, 3, 4, 5, 6, 7];
        sets = {C_idx_1, C_idx_2};
    elseif K == 8
        % Use Set 1 indices matching user script Is[2] to Is[9] (s2 to s9)
        sets = {[3, 4, 5, 6, 7, 8, 9, 10]};
    else
        % K = 5
        sets = {[1, 2, 3, 4, 5]};
    end
    
    for set_id = 1:length(sets)
        C_idx = sets{set_id};
        U = length(C_idx); % = K
        
        aoi_analytic_periodic = zeros(1, U);
        aoi_analytic_gatw = zeros(1, U);
        aoi_mc_periodic = zeros(1, U);
        aoi_mc_gatw = zeros(1, U);
        
        for user_idx = 1:U
            %% 解析公式計算
            aoi_analytic_periodic(user_idx) = compute_aoi_periodic(Is, C_idx, user_idx, K, L);
            aoi_analytic_gatw(user_idx) = compute_aoi_gatw(Is, C_idx, user_idx, K, L);
            
            %% MC 模擬計算
            aoi_mc_periodic(user_idx) = simulate_aoi_mc(Is, C_idx, user_idx, K, L, 'periodic', num_runs, N_period);
            aoi_mc_gatw(user_idx) = simulate_aoi_mc(Is, C_idx, user_idx, K, L, 'gatw', num_runs, N_period);
            
            fprintf('Set %d User %d done.\n', set_id, user_idx);
        end
        
        % 取得序列名稱 s_g (0-indexed)
        seq_names = cell(1, U);
        for j = 1:U
            seq_names{j} = sprintf('s%d', C_idx(j) - 1);
        end
        
        %% 整理結果並繪圖
        results_periodic = struct();
        results_periodic.analytical = aoi_analytic_periodic;
        results_periodic.mc = aoi_mc_periodic;
        results_periodic.seq_names = seq_names;
        results_periodic.set_id = set_id;
        
        results_gatw = struct();
        results_gatw.analytical = aoi_analytic_gatw;
        results_gatw.mc = aoi_mc_gatw;
        results_gatw.seq_names = seq_names;
        results_gatw.set_id = set_id;
        
        plot_results(results_periodic, K, 'Periodic (T=L)');
        plot_results(results_gatw, K, 'Generate-at-will (T=1)');
    end
end
