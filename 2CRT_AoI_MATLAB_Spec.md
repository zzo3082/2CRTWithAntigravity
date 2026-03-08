# 2-CRT 序列 AAoI 模擬規格文件

> 本文件供 Claude Code 實作使用，對應論文：
> *User-Irrepressible Sequences for Multiple-Packet Reception with MPR Capability 2: Constructions and Age of Information*

---

## 1. 專案目標

根據論文中的解析公式，計算 K=5（prime）與 K=6（composite）兩種情境下，2-CRT 序列在以下兩種流量模型的 AAoI：

- **Generate-at-will traffic**（論文 Section IV-A，公式 (23)）
- **Periodic traffic**（論文 Section IV-B，公式 (27)）

輸出對應論文 **Table I** 格式的數值表格，並畫出每位使用者的 AAoI 比較圖。

---

## 2. 目錄結構

```
project/
├── main.m                  % 主程式：呼叫所有模組、輸出結果
├── build_2CRT.m            % 模組1：建構 2-CRT(K) 特徵集合
├── compute_hamming.m       % 模組2：計算 Hamming cross-correlation
├── compute_e.m             % 模組3：計算 e(A, B)（公式 (20)(21)(22)）
├── compute_aoi_periodic.m  % 模組4：Periodic traffic AAoI（公式 (27)(28)(29)(30)(31)(32)(33)）
├── compute_aoi_gatw.m      % 模組5：Generate-at-will AAoI（公式 (23)(24)(26)）
├── plot_results.m          % 模組6：繪圖與表格輸出
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

### 4.6 `plot_results.m`

**功能**：輸出數值表格（對應 Table I）與 AAoI 比較圖。

**輸入**：
- `results`（struct）：包含所有使用者的 AAoI 結果
- `K`（integer）
- `mode`（string）：`'periodic'` 或 `'gatw'`

**輸出**：
- 印出 Table I 格式的數值表（含各序列 AAoI 與平均值）
- 畫出各序列 AAoI 的 bar chart 或折線圖

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

```
%% 設定參數
K_list = [5, 6];

for each K in K_list:
  %% Step 1: 建構序列
  [Is, L, pK] = build_2CRT(K);

  %% Step 2: 選取 Set 1 與 Set 2（對應 Table I）
  % K=6 時：
  %   Set 1 = {s0, s2, s3, s4, s5, s7} → C_idx = [1, 3, 4, 5, 6, 8]（1-indexed）
  %   Set 2 = {s1, s2, s3, s4, s5, s6} → C_idx = [2, 3, 4, 5, 6, 7]
  % K=5 時：任選 5 個序列（共 6 個可選）

  %% Step 3: 計算每個使用者的 AAoI（periodic + generate-at-will）
  for each set:
    for each user in set:
      aoi_periodic(user) = compute_aoi_periodic(Is, C_idx, user, K, L);
      aoi_gatw(user)     = compute_aoi_gatw(Is, C_idx, user, K, L);
    end
  end

  %% Step 4: 輸出結果
  plot_results(results, K, 'periodic');
  plot_results(results, K, 'gatw');
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

### Test 3：AAoI 數值驗證（K=6）

| 測試項目 | 容許誤差 |
|---|---|
| Set 1 generate-at-will 平均 = 5.871 | ±0.01 |
| Set 2 generate-at-will 平均 = 8.008 | ±0.01 |
| Set 1 periodic 平均 = 26.427 | ±0.01 |
| Set 2 periodic 平均 = 28.529 | ±0.01 |

### Test 4：e(A, B) 計算

- 手動驗算 K=5 中兩個 Type II 序列對的 e 值，對照公式 (21)
- 確認 e(Is_g, Is_h) = i*(K-i) 其中 i = (h-g)^{-1} in Fp

---

## 8. 注意事項與常見陷阱

1. **Index 起始**：論文公式使用 0-indexed，MATLAB 預設 1-indexed。建議特徵集合元素保持 0-indexed 儲存，僅 cell array 的 index 用 1-indexed。

2. **K=5 是質數**：K=5 時直接用 Definition 2，序列長度 L = p^2 = 25（不是 pK * K）。K=6 時用 Definition 3，L = 7 * 6 = 42。

3. **Type I / Type II 判斷**：s0 和 spK（最後一個序列）通常是 Type I。建議在 build 後對每對序列預先算好 Hamming 值並儲存為矩陣。

4. **公式 (31) 的子集枚舉**：`{x1,...,xn}` 的非空子集數為 2^n - 1，K=6 時最多 2^6-1=63 個，計算量可接受。

5. **公式 (26) 的計算複雜度**：對 K=6，需枚舉 L 個 τ 值 × 多種 collision 組合，建議預先計算並快取 αj(A) 的值。

6. **mod 運算**：MATLAB 的 `mod()` 回傳非負值，符合 Z_L 的需求。但注意 `mod(a - b, L)` 計算差分時需確保結果在 [0, L-1]。

---

## 9. 參考公式索引

| 公式編號 | 內容 | 使用位置 |
|---|---|---|
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
