# 2-CRT 序列 AAoI 模擬規格文件

> 本文件供 Claude Code 實作使用，對應論文：
> *User-Irrepressible Sequences for Multiple-Packet Reception with MPR Capability 2: Constructions and Age of Information*

---

## 1. 專案目標

計算 K=5（prime）與 K=6（composite）兩種情境下，2-CRT 序列在以下兩種流量模型的 AAoI：

- **Generate-at-will traffic**（論文 Section IV-A，公式 (23)）
- **Periodic traffic**（論文 Section IV-B，公式 (27)）

使用**兩種方法**計算並互相對照：

| 方法 | 說明 | 對應模組 |
|---|---|---|
| **解析公式法** | 直接套用論文推導的閉合公式 | `compute_aoi_periodic.m`, `compute_aoi_gatw.m` |
| **蒙地卡羅模擬法** | 隨機抽取 τ offset，模擬 AoI 時間序列，取 10^6 次平均 | `simulate_aoi_mc.m` |

兩種方法的結果應在容許誤差內吻合，輸出對應論文 **Table I** 格式的數值對照表與折線圖。

---

## 2. 目錄結構

```
project/
├── main.m                  % 主程式：呼叫所有模組、輸出對照結果
├── build_2CRT.m            % 模組1：建構 2-CRT(K) 特徵集合
├── compute_hamming.m       % 模組2：計算 Hamming cross-correlation
├── compute_e.m             % 模組3：計算 e(A, B)（公式 (20)(21)(22)）
├── compute_aoi_periodic.m  % 模組4：Periodic traffic AAoI 解析公式
├── compute_aoi_gatw.m      % 模組5：Generate-at-will AAoI 解析公式
├── simulate_aoi_mc.m       % 模組6：蒙地卡羅模擬（兩種 traffic 通用）
├── plot_results.m          % 模組7：繪圖與表格輸出（含解析 vs MC 對照）
└── utils/
    ├── next_prime.m        % 工具：找比 K 大的最小質數 pK
    ├── mod_inv.m           % 工具：在 Fp 中計算模反元素
    └── crt_map.m           % 工具：CRT 映射 θ(t) = (t mod pK, t mod K)
```

---

## 3. 數學符號對照表

| MATLAB 變數名 | 論文符號 | 說明 |
|---|---|---|
| `K` | K | 使用者數量 |
| `pK` | p_K | 比 K 大的最小質數（K=5 時 pK=5，K=6 時 pK=7） |
| `L` | L | 序列長度 = pK * K |
| `Is{g}` | I_{s_g} | 第 g 個序列的特徵集合（1-indexed，g=1..K+1） |
| `delta` | δ | 特徵集合的連續差分向量 |
| `C_delta` | C_δ(d) | 自相關函數（公式 (28)） |
| `h_n` | h_n | 公式 (30) 的係數 |
| `h_hat_n` | ĥ_n | 公式 (31) 的累積分布值 |
| `e_val` | e(A,B) | 公式 (20) 的重疊計數 |

---

## 4. 模組規格

### 4.1 `build_2CRT.m`

**功能**：建構 K 個使用者的 2-CRT 序列特徵集合。

**輸入**：
- `K`（integer）：使用者數量

**輸出**：
- `Is`（cell array, size 1×(K+1)）：`Is{g}` 為第 g 個序列特徵集合，元素為 Z_L 的子集（0-indexed 整數集合）
- `L`（integer）：序列長度 = pK * K
- `pK`（integer）：比 K 大的最小質數

**邏輯**：

```
1. 用 next_prime.m 求 pK（K 為質數時 pK = K，K 為合數時 pK = 比 K 大的最小質數）
   ※ 注意：K=5 是質數，直接用 Definition 2（序列長度 L = p^2 = 25）
            K=6 是合數，用 Definition 3（序列長度 L = pK * K = 7*6 = 42）

2. [K 為質數時 - Definition 2]
   for g = 0 to p-1:
     Is{g+1} = { mod(j * (1 + g*p), p^2) : j = 0,1,...,p-1 }
   Is{p+1} = { mod(j * p, p^2) : j = 0,1,...,p-1 }

3. [K 為合數時 - Definition 3]
   使用 CRT 映射 θ: Z_{pK*K} → Z_{pK} × Z_K
   for g = 0 to pK-1:
     Ibs_g = { mod(j * [g, 1], [pK, K]) : j = 0,...,K-1 }  (pair-wise mod)
     Is{g+1} = { t ∈ Z_L : crt_map(t) ∈ Ibs_g }
   Ibs_{pK} = { mod(j * [1, 0], [pK, K]) : j = 0,...,K-1 }
   Is{pK+1} = { t ∈ Z_L : crt_map(t) ∈ Ibs_{pK} }
```

**驗證**（build 後自動 assert）：
- 每個 `Is{g}` 的元素個數 = K
- 所有特徵集合的聯集不超過 Z_L
- 對應 Example 1（K=5）與 Example 2（K=6）的特徵集合值需與論文一致

---

### 4.2 `compute_hamming.m`

**功能**：計算兩序列之間的 Hamming cross-correlation H_{si, sj}。

**輸入**：
- `Isi`（array）：第 i 個序列特徵集合
- `Isj`（array）：第 j 個序列特徵集合
- `L`（integer）：序列長度

**輸出**：
- `H`（integer）：max over τ of |Isi ∩ (Isj + τ)|

**邏輯**：
```
H = 0
for τ = 0 to L-1:
  shifted = mod(Isj + τ, L)
  overlap = length(intersect(Isi, shifted))
  H = max(H, overlap)
```

---

### 4.3 `compute_e.m`

**功能**：計算 e(A, B)（公式 (20)(21)(22)），即 |A ∩ (B + τ)| = 2 的 τ 數量。

**輸入**：
- `A`（array）：子集，來自某個序列特徵集合（可以是 Is_g 的子集）
- `B`（array）：另一序列 Is_h 的特徵集合
- `L`（integer）：序列長度
- `K`（integer）：使用者數量

**輸出**：
- `e_val`（integer）：符合 |A ∩ (B+τ)| = 2 的 τ 個數

**邏輯**：
```
使用公式 (21)(22)：
1. 計算 d*(A) 與 d*(B) 的交集，找出 {±x}
2. 若交集為空，e_val = 0
3. 否則，x = |i*g| 的形式，找出對應的 i
4. 計算 f(A; x) 與 f(B; x)（公式 (19)）
5. e_val = f(A; x) * f(B; x)（公式 (20)）
   若 A 是子集，用公式 (22)
```

---

### 4.4 `compute_aoi_periodic.m`

**功能**：計算 Periodic traffic 下使用者 i 的 AAoI（公式 (27)）。

**輸入**：
- `Is`（cell array）：所有序列特徵集合
- `C_idx`（array）：被選用的 K 個序列的 index（從 Is 中選）
- `user_idx`（integer）：要計算 AAoI 的使用者在 C_idx 中的位置
- `K`（integer）：使用者數量
- `L`（integer）：序列長度

**輸出**：
- `AoI`（double）：該使用者的 AAoI

**邏輯**：

```
公式 (27)：AoI = E[S_i^q] + L/2 - 1/2

[Type I 序列（Hsg,sh = 1 for all h ≠ g）]
用公式 (28)(29)(30)：
  1. 計算 delta = 特徵集合的連續差分
  2. 計算 C_delta(d) = sum_{j=0}^{K-1} delta_j * delta_{j-d}（公式 (28)）
  3. 計算 h_n（公式 (30)）：hn = sum_{i=1}^{n} (-1)^{i+1} * C(n,i) * (1 - i/pK)^{K-1}
  4. 計算 E[S_i^q]（公式 (29)）

[Type II 序列（存在 Hsg,sh = 2 的情況）]
用公式 (31)(32)(33)：
  1. 對每個 τ ∈ Z_L，計算 I_sg^τ = {x1,...,xK}（排序後）
  2. 計算 h_hat_n（公式 (31)）：
     對所有非空子集 A ⊆ {x1,...,xn}：
       h_hat_n += (-1)^{|A|+1} * prod_{sh ∈ C\sg} (1 - (|A|*K - e(A,sh)) / L)
  3. 計算 E[S_i^q | τ]（公式 (32)）
  4. 對所有 τ 平均得 E[S_i^q]（公式 (33)）
```

**注意**：Type I / Type II 的判斷由 `compute_hamming.m` 預先計算，max cross-correlation = 1 為 Type I，否則為 Type II。

---

### 4.5 `compute_aoi_gatw.m`

**功能**：計算 Generate-at-will traffic 下使用者 i 的 AAoI（公式 (23)）。

**輸入**：
- `Is`（cell array）：所有序列特徵集合
- `C_idx`（array）：被選用的 K 個序列的 index
- `user_idx`（integer）：目標使用者在 C_idx 中的位置
- `K`（integer）
- `L`（integer）

**輸出**：
- `AoI`（double）

**邏輯**：
```
公式 (23)：AoI = E[(Y_i^q)^2] / (2 * E[Y_i^q]) - 1/2

1. 計算 Ti（吞吐量）：
   Ti = (成功傳輸的 slot 數) / L
   需計算在所有可能的 τ offset 下，使用者 i 平均成功傳輸幾次

2. E[Y_i^q] = L / Ti（公式 (24)）

3. 計算 E[(Y_i^q)^2]（公式 (26)）：
   對每個 τ ∈ Z_L 與每種其他使用者 collision 組合：
     - 計算 reverse Hamming cross-correlation αj(A)
     - 找出成功傳輸的 slot 集合 I_si^τ ∩ ((Is_i)^2 \ B)
     - 計算 η(I) = 各成功傳輸間隔的平方和
   加總除以 L^{K+1}，再除以 Ti

4. AoI = E[(Y_i^q)^2] / (2 * E[Y_i^q]) - 1/2
```

---

### 4.6 `simulate_aoi_mc.m`

**功能**：用蒙地卡羅模擬法計算單一使用者的 AAoI，支援 generate-at-will 與 periodic 兩種 traffic。

**輸入**：
- `Is`（cell array）：所有序列特徵集合
- `C_idx`（array）：被選用的 K 個序列的 index
- `user_idx`（integer）：目標使用者在 C_idx 中的位置
- `K`（integer）：使用者數量
- `L`（integer）：序列長度
- `traffic`（string）：`'gatw'` 或 `'periodic'`
- `num_runs`（integer）：隨機 τ 抽樣次數，預設 `2000`
- `N_period`（integer）：每次 run 模擬的週期數，預設 `500`

**輸出**：
- `AoI_mc`（double）：該使用者的模擬 AAoI

**核心邏輯**：

```
AoI_total = 0;

for run = 1 to num_runs:

  %% Step 1: 隨機抽取各使用者的相對時間偏移量 τ
  % 目標使用者 i 偏移設為 0（不失一般性）
  % 其他 K-1 個使用者各自均勻隨機抽 tau_j ∈ Z_L
  tau_rel = [0, randi([0, L-1], 1, K-1)];

  %% Step 2: 建立長度 T_sim = N_period * L 的傳輸排程矩陣
  % tx(j, t) = 1 iff 使用者 j 在 slot t 傳輸
  % 即 mod(t, L) ∈ mod(Is{C_idx(j)} + tau_rel(j), L)
  % 利用 MATLAB 向量化：
  %   t_vec = 0 : T_sim-1
  %   for each user j: tx(j,:) = ismember(mod(t_vec, L), shifted_Is{j})

  %% Step 3: 計算每個 slot 的碰撞狀況
  % total_tx(t) = sum(tx(:, t))  % 同時傳輸的使用者總數
  % success_i(t) = tx(i, t) AND total_tx(t) <= gamma  % gamma=2

  %% Step 4: 逐 slot 計算 AoI 時間序列 Ai(t)
  % [Generate-at-will]
  %   生成時間 g_t = t（每次傳輸時才生成新 packet）
  %   Ai(t) 初始化為大數（如 T_sim）
  %   若 success_i(t)=1：Ai(t) = 0（S_i^q = 0，立即生成立即送）
  %   否則：Ai(t) = Ai(t-1) + 1

  % [Periodic]
  %   生成時間 g_t = floor(t / L) * L（每個週期開頭生成）
  %   若 success_i(t)=1：Ai(t) = t - g_t
  %   否則：Ai(t) = Ai(t-1) + 1
  %   注意：若整個週期都沒成功，AoI 會持續累積到下個週期

  %% Step 5: 捨棄前 warm-up 期（建議捨棄前 10% 的 slot）
  AoI_run = mean(Ai_series(warmup_end:end));
  AoI_total = AoI_total + AoI_run;

end

AoI_mc = AoI_total / num_runs;
```

**碰撞判斷（MPR capability γ=2）**：

```
在 slot t：
  同時傳輸的使用者總數 = total_tx(t)
  成功條件：tx(i,t)=1 且 total_tx(t) <= gamma（= 2）
  即：使用者 i 傳輸，且同時傳輸的人數（含自己）≤ 2
```

**效能建議**：

- 預先用 `ismember` 向量化建構 `tx` 矩陣，避免雙層 for loop
- `N_period=500`、`num_runs=2000` 在 K=6 時約 30 秒內可完成
- 可加 `parfor` 加速外層 run loop（需 Parallel Computing Toolbox）

---

### 4.7 `plot_results.m`

**功能**：輸出解析公式與蒙地卡羅模擬的對照表格與圖形。

**輸入**：
- `results`（struct）：包含以下欄位
  - `results.analytical`：各使用者解析 AAoI（K×1 array）
  - `results.mc`：各使用者蒙地卡羅 AAoI（K×1 array）
  - `results.seq_names`：序列名稱（如 `{'s0','s2','s3','s4','s5','s7'}`）
- `K`（integer）
- `traffic`（string）：`'periodic'` 或 `'gatw'`

**輸出**：

1. **數值對照表**（印到 command window）：

```
=== K=6, Periodic (T=L), Set 1 ===
         s0      s2      s3      s4      s5      s7    Avg
Analyt: 25.711  27.072  26.568  26.250  26.908  26.054 26.427
MC:     25.69   27.05   26.55   26.24   26.89   26.03  26.41
Diff:    0.021   0.022   0.018   0.010   0.018   0.024  0.017
```

2. **折線圖**：x 軸為序列編號，y 軸為 AAoI，分別畫解析值（實線）與 MC 值（虛線＋圓點），同一張圖方便對照。圖例標示 `Analytical` 與 `Monte Carlo (N=2000)`。

---

### 4.7 `utils/next_prime.m`

```matlab
% 回傳比 n 大的最小質數（若 n 本身為質數則回傳 n）
function p = next_prime(n)
```

### 4.8 `utils/mod_inv.m`

```matlab
% 在 Fp（模 p 的有限域）中計算 a 的乘法反元素
% 使用 Extended Euclidean Algorithm
function inv = mod_inv(a, p)
```

### 4.9 `utils/crt_map.m`

```matlab
% CRT 映射：θ(t) = (t mod pK, t mod K)
% 輸入 t（scalar 或 array），回傳 [t mod pK, t mod K]
function out = crt_map(t, pK, K)
```

---

## 5. 主程式 `main.m` 流程

```matlab
%% 設定參數
K_list   = [5, 6];
num_runs = 2000;    % MC 隨機 τ 抽樣次數
N_period = 500;     % 每次 run 模擬的週期數
gamma    = 2;       % MPR capability

for K = K_list

  %% Step 1: 建構序列
  [Is, L, pK] = build_2CRT(K);

  %% Step 2: 定義要比較的 set（對應 Table I）
  % K=6: Set1={s0,s2,s3,s4,s5,s7}, Set2={s1,s2,s3,s4,s5,s6}
  % K=5: Set1=前5個序列, Set2=後5個序列（共6個）
  sets = define_sets(K);   % 回傳 cell array of C_idx

  for set_id = 1:length(sets)
    C_idx = sets{set_id};

    for user_idx = 1:K

      %% Step 3a: 解析公式
      aoi_analytic_periodic(user_idx) = compute_aoi_periodic(Is, C_idx, user_idx, K, L, pK);
      aoi_analytic_gatw(user_idx)     = compute_aoi_gatw(Is, C_idx, user_idx, K, L);

      %% Step 3b: 蒙地卡羅模擬
      aoi_mc_periodic(user_idx) = simulate_aoi_mc(Is, C_idx, user_idx, K, L, 'periodic', num_runs, N_period);
      aoi_mc_gatw(user_idx)     = simulate_aoi_mc(Is, C_idx, user_idx, K, L, 'gatw',     num_runs, N_period);

    end

    %% Step 4: 整理結果並輸出
    results_periodic.analytical = aoi_analytic_periodic;
    results_periodic.mc         = aoi_mc_periodic;
    results_periodic.seq_names  = get_seq_names(C_idx);

    results_gatw.analytical = aoi_analytic_gatw;
    results_gatw.mc         = aoi_mc_gatw;
    results_gatw.seq_names  = get_seq_names(C_idx);

    plot_results(results_periodic, K, 'periodic');
    plot_results(results_gatw,     K, 'gatw');

  end
end
```

---

## 6. 預期輸出（對應論文 Table I，K=6）

```
=== K=6, Generate-at-will (T=1) ===
         s0      s2      s3      s4      s5      s7    Avg
Set 1:  5.163   6.505   6.019   5.699   6.345   5.493  5.871
Set 2:  17.205  6.529   6.067   5.726   6.383   6.141  8.008

=== K=6, Periodic (T=L) ===
         s0      s2      s3      s4      s5      s7    Avg
Set 1:  25.711  27.072  26.568  26.250  26.908  26.054 26.427
Set 2:  37.423  27.137  26.635  26.279  26.970  26.731 28.529
```

---

## 7. 測試計畫

### Test 1：序列建構正確性

| 測試項目 | 預期結果 |
|---|---|
| K=5，`Is{1}` = Is0 | `{0,1,2,3,4}` |
| K=5，`Is{2}` = Is1 | `{0,6,12,18,24}` |
| K=5，`Is{6}` = Is5 | `{0,5,10,15,20}` |
| K=4，`Is{1}` = Is0 | `{0,5,10,15}` |
| K=4，`Is{2}` = Is1 | `{0,1,2,3}` |
| 每個集合元素個數 | = K |
| 序列長度 L（K=5）| 25 |
| 序列長度 L（K=6）| 42 |

### Test 2：Hamming Cross-Correlation

| 測試項目 | 預期結果 |
|---|---|
| K=5，H(s0, s5) | 1（Type I pair） |
| K=5，H(s0, s1) | 2（Type II pair） |
| K=6，H(s0, s1) | 1（s0 為 Type I） |
| K=6，H(s1, s2) | ≤ 2 |

### Test 3：解析公式數值驗證（K=6，對照論文 Table I）

| 測試項目 | 容許誤差 |
|---|---|
| Set 1 generate-at-will 平均 = 5.871 | ±0.01 |
| Set 2 generate-at-will 平均 = 8.008 | ±0.01 |
| Set 1 periodic 平均 = 26.427 | ±0.01 |
| Set 2 periodic 平均 = 28.529 | ±0.01 |

### Test 4：蒙地卡羅模擬收斂性

| 測試項目 | 預期結果 |
|---|---|
| MC 結果與解析值的差（每個使用者） | ≤ 0.05（num_runs=2000） |
| MC 結果與解析值的差（平均） | ≤ 0.02 |
| 增加 num_runs 後誤差應縮小 | 驗證收斂性 |

**收斂性驗證方法**：對 K=6 Set 1，分別跑 `num_runs = [500, 1000, 2000, 5000]`，畫出 MC 誤差 vs num_runs 的折線圖，確認趨近於 0。

### Test 5：e(A, B) 計算

- 手動驗算 K=5 中兩個 Type II 序列對的 e 值，對照公式 (21)
- 確認 e(Is_g, Is_h) = i*(K-i) 其中 i = (h-g)^{-1} in Fp

### Test 6：2-UI 性質驗證（sanity check）

對 K=6 的任意 K-subset，以蒙地卡羅方式驗證：在 10^4 組隨機 τ 中，每位使用者都至少有 1 個成功傳輸 slot（即碰撞數 ≤ 1 的 slot 存在）。

---

## 8. 注意事項與常見陷阱

1. **Index 起始**：論文公式使用 0-indexed，MATLAB 預設 1-indexed。建議特徵集合元素保持 0-indexed 儲存，僅 cell array 的 index 用 1-indexed。

2. **K=5 是質數**：K=5 時直接用 Definition 2，序列長度 L = p^2 = 25（不是 pK * K）。K=6 時用 Definition 3，L = 7 * 6 = 42。

3. **Type I / Type II 判斷**：s0 和 spK（最後一個序列）通常是 Type I。建議在 build 後對每對序列預先算好 Hamming 值並儲存為矩陣。

4. **公式 (31) 的子集枚舉**：`{x1,...,xn}` 的非空子集數為 2^n - 1，K=6 時最多 2^6-1=63 個，計算量可接受。

5. **公式 (26) 的計算複雜度**：對 K=6，需枚舉 L 個 τ 值 × 多種 collision 組合，建議預先計算並快取 αj(A) 的值。

6. **mod 運算**：MATLAB 的 `mod()` 回傳非負值，符合 Z_L 的需求。但注意 `mod(a - b, L)` 計算差分時需確保結果在 [0, L-1]。

7. **MC 的 AoI 初始值**：模擬開始時 AoI 初始化為 0 即可，搭配 warm-up 捨棄可避免初始值影響結果。建議捨棄前 `ceil(0.1 * N_period * L)` 個 slot。

8. **Generate-at-will 的 AoI 定義**：論文設 S_i^q = 0（packet 在傳輸當下生成），所以成功傳輸時 Ai(t) 直接歸零，而非計算 t - generation_time。

9. **Periodic traffic 的邊界條件**：若一個週期內使用者 i 完全沒有成功傳輸（可能因為每次都碰撞），AoI 會從上個週期繼續累積，不能自動重置。這是 2-CRT 在 Set 2（含 s1）時 AoI 較高的原因（s1 是 Type II 且 cross-correlation 較高）。

10. **MC 與解析值不一致時的排查順序**：
    - 先確認序列建構正確（Test 1）
    - 再確認碰撞判斷邏輯（γ=2 表示同時傳輸 ≤ 2 人才成功）
    - 最後確認 AoI 更新公式對應正確的 traffic 模型

---

## 9. 參考公式索引

| 公式編號 | 內容 | 使用位置 |
|---|---|---|
| (2) | AoI 更新規則 Ai(t) | `simulate_aoi_mc.m` |
| (3) | AAoI 時間平均定義 | `simulate_aoi_mc.m` |
| (19) | f(A; x) 定義 | `compute_e.m` |
| (20) | e(A, B) 定義 | `compute_e.m` |
| (21) | e(Is_g, Is_h) 解析式 | `compute_e.m` |
| (22) | e(A, Is_h) 子集版 | `compute_e.m` |
| (23) | Generate-at-will AAoI | `compute_aoi_gatw.m` |
| (24) | E[Y] = L/T | `compute_aoi_gatw.m` |
| (26) | E[Y^2] 計算 | `compute_aoi_gatw.m` |
| (27) | Periodic AAoI | `compute_aoi_periodic.m` |
| (28) | C_δ(d) 定義 | `compute_aoi_periodic.m` |
| (29) | Type I E[S] | `compute_aoi_periodic.m` |
| (30) | h_n 係數 | `compute_aoi_periodic.m` |
| (31) | ĥ_n（Type II CDF） | `compute_aoi_periodic.m` |
| (32) | E[S\|τ] | `compute_aoi_periodic.m` |
| (33) | E[S] 平均 | `compute_aoi_periodic.m` |
