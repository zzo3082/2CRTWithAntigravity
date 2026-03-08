# 2-CRT Sequences for Multiple-Packet Reception (AAoI Simulation)

本專案提供了一套基於 MATLAB 的模擬框架，用於實作並驗證多重封包接收（Multiple-Packet Reception, MPR $\gamma=2$）能力下的 **2-CRT 序列 (2-Chinese Remainder Theorem sequences)** 的平均資訊更新年齡（Average Age of Information, 簡稱 AAoI）。

相關理論與推導對應於論文：
*User-Irrepressible Sequences for Multiple-Packet Reception with MPR Capability 2: Constructions and Age of Information*

## 🌟 核心特色

本模擬專案不僅遵循了論文的理論基礎，還實作了**快速機率分佈卷積法 (Probability Distribution Convolution)**：
- 避免了直接嵌套 6 層以上的指數級碰撞窮舉。
- 完美適配論文中複雜條件下的 **Generate-at-will traffic** 及 **Periodic traffic** 模型。
- 高效且精準地計算由於網路延遲 $\tau$ 所引發的使用者之間互相干擾與成功傳送間距 (Gaps) 之期望值 $E[W]$ 與 $E[\eta]$。

## 📂 專案結構

```text
├── main.m                  % 主程式：呼叫所有模組並輸出最終報表與圖表
├── build_2CRT.m            % 建構模組：建立 2-CRT(K) 序列的特徵集合
├── compute_hamming.m       % 計算兩個序列之間的 Hamming cross-correlation
├── compute_e.m             % 計算 Subset 碰撞下的 e 函數計數
├── compute_aoi_periodic.m  % 計算 Periodic Traffic 下的 AAoI 值
├── compute_aoi_gatw.m      % 計算 Generate-at-will Traffic 下的 AAoI 值
├── plot_results.m          % 繪製控制台數值表格與對照長條圖
├── 2CRT_AoI_MATLAB_Spec.md % 實作與推導原始規格文件
└── utils/                  % 工具庫
    ├── next_prime.m        % 尋找下一個質數 
    ├── mod_inv.m           % 模反元素運算
    ├── crt_map.m           % CRT 參數座標對應 (t mod pK, t mod K)
    └── compute_h_hat.m     % (Type II) 期望次數輔助運算
```

## 🚀 執行方式

1. 確保您的電腦上已安裝 **MATLAB**。
2. 使用 MATLAB 開啟本專案目錄。
3. 在指令視窗 (Command Window) 或直接開啟 `main.m`，並點擊執行 (Run)。
4. 程式將會在命令列自動印出針對 $K=5$ 與 $K=6$ 的對照表格（與論文 Table I 吻合），並彈出長條圖 (Bar charts)。

```matlab
% 直接在 MATLAB 輸入以下指令即可執行
>> main
```

## 📊 驗證與預期結果 (K = 6 範例)

本模擬產出的數據與論文中的表格展示數據精準對應：

```text
=== K=6, Generate-at-will (T=1) ===
         s0      s2      s3      s4      s5      s7    Avg
Set 1:  5.163   6.505   6.019   5.699   6.345   5.493  5.871
Set 2:  17.205  6.529   6.067   5.726   6.383   6.141  8.008

=== K=6, Periodic (T=L) ===
         s0      s2      s3      s4      s5      s7    Avg
Set 1:  25.711  27.072  26.568  26.250  26.908  26.054 26.427
Set 2:  37.423  27.137  26.635  26.279  26.970  26.731 28.529
```
