function AoI_mc = simulate_aoi_mc(Is, C_idx, user_idx, K, L, traffic, num_runs, N_period)
    % simulate_aoi_mc: 用蒙地卡羅模擬法計算單一使用者的 AAoI
    % 支援 generate-at-will 與 periodic 兩種 traffic 模型
    
    gamma = 2; % MPR capability = 2
    AoI_total = 0;
    
    T_sim = N_period * L;
    warmup_end = ceil(0.1 * T_sim); % 捨棄前 10% 的 slot 做為 warmup
    
    % 取得所有使用者的序列集合
    U = length(C_idx);
    seqs = cell(1, U);
    for j = 1:U
        seqs{j} = Is{C_idx(j)};
    end
    my_user = user_idx;
    
    for run = 1:num_runs
        % Step 1: 隨機抽取各使用者的相對時間偏移量 τ
        % 目標使用者偏移設為 0
        tau_rel = zeros(1, U);
        for j = 1:U
            if j ~= my_user
                tau_rel(j) = randi([0, L-1]);
            end
        end
        
        % Step 2: 建立傳輸排程矩陣
        % tx(j, t) = 1 iff 使用者 j 在 slot t 傳輸
        t_vec = 0:(T_sim-1);
        tx = zeros(U, T_sim);
        for j = 1:U
            shifted_seq = mod(seqs{j} + tau_rel(j), L);
            % 判斷每個 slot t 是否在這個使用者平移後的特徵集合中
            tx(j, :) = ismember(mod(t_vec, L), shifted_seq);
        end
        
        % Step 3: 計算每個 slot 的碰撞狀況
        total_tx = sum(tx, 1); % 計算每一個 slot 的總傳輸人數
        % 成功傳輸條件：自己有傳，且總人數 <= gamma (MPR limit)
        success_i = (tx(my_user, :) == 1) & (total_tx <= gamma);
        
        % Step 4: 逐 slot 計算 AoI 時間序列 Ai(t)
        Ai_series = zeros(1, T_sim);
        current_AoI = 0; % MC 初始值從 0 開始，依靠 warmup 洗掉
        
        if strcmp(traffic, 'gatw')
            for t = 1:T_sim
                if success_i(t)
                    % Generate-at-will: 在發送的當下生成新封包 (S_i^q = 0)
                    % 一旦發送成功，Age 瞬間降為 0
                    current_AoI = 0;
                else
                    current_AoI = current_AoI + 1;
                end
                Ai_series(t) = current_AoI;
            end
            
        elseif strcmp(traffic, 'periodic')
            for t = 1:T_sim
                % slot count 從 0 開始，對應 t-1
                real_t = t - 1;
                % 每個週期 L 的開頭會產生新的 packet（generation time）
                g_t = floor(real_t / L) * L;
                
                if success_i(t)
                    % 如果成功傳輸，AoI 會降到 (目前時間 - 生成時間)
                    current_AoI = real_t - g_t;
                else
                    current_AoI = current_AoI + 1;
                end
                Ai_series(t) = current_AoI;
            end
        end
        
        % Step 5: 捨棄前 warm-up 期，取平均加總
        AoI_run = mean(Ai_series(warmup_end:end));
        AoI_total = AoI_total + AoI_run;
    end
    
    % 取所有 run 的平均
    AoI_mc = AoI_total / num_runs;
end
