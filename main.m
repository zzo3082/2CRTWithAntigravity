% 2-CRT 序列 AAoI 主程式
clc; clear; close all;

% 加入 utils 資料夾路徑
addpath('utils');

K_list = [5, 6];

for K = K_list
    [Is, L, pK] = build_2CRT(K);
    fprintf('=== Building 2-CRT for K = %d (Length = %d) ===\n', K, L);
    
    results_gatw = struct();
    results_periodic = struct();
    
    if K == 6
        % Set 1: s0, s2, s3, s4, s5, s7 
        C_idx_1 = [1, 3, 4, 5, 6, 8];
        % Set 2: s1, s2, s3, s4, s5, s6
        C_idx_2 = [2, 3, 4, 5, 6, 7];
        
        aoi_gatw_1 = zeros(1, 6);
        aoi_periodic_1 = zeros(1, 6);
        for user_idx = 1:length(C_idx_1)
            aoi_gatw_1(user_idx) = compute_aoi_gatw(Is, C_idx_1, user_idx, K, L);
            aoi_periodic_1(user_idx) = compute_aoi_periodic(Is, C_idx_1, user_idx, K, L);
            fprintf('Set 1 User %d done.\n', user_idx);
        end
        results_gatw.set1 = aoi_gatw_1;
        results_periodic.set1 = aoi_periodic_1;
        
        aoi_gatw_2 = zeros(1, 6);
        aoi_periodic_2 = zeros(1, 6);
        for user_idx = 1:length(C_idx_2)
            aoi_gatw_2(user_idx) = compute_aoi_gatw(Is, C_idx_2, user_idx, K, L);
            aoi_periodic_2(user_idx) = compute_aoi_periodic(Is, C_idx_2, user_idx, K, L);
            fprintf('Set 2 User %d done.\n', user_idx);
        end
        results_gatw.set2 = aoi_gatw_2;
        results_periodic.set2 = aoi_periodic_2;
        
    else
        % K = 5
        C_idx_1 = [1, 2, 3, 4, 5];
        aoi_gatw_1 = zeros(1, 5);
        aoi_periodic_1 = zeros(1, 5);
        for user_idx = 1:length(C_idx_1)
            aoi_gatw_1(user_idx) = compute_aoi_gatw(Is, C_idx_1, user_idx, K, L);
            aoi_periodic_1(user_idx) = compute_aoi_periodic(Is, C_idx_1, user_idx, K, L);
            fprintf('Set 1 User %d done.\n', user_idx);
        end
        results_gatw.set1 = aoi_gatw_1;
        results_periodic.set1 = aoi_periodic_1;
    end
    
    plot_results(results_gatw, K, 'Generate-at-will (T=1)');
    plot_results(results_periodic, K, 'Periodic (T=L)');
end
