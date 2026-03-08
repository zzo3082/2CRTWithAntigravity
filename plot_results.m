function plot_results(results, K, mode)
    % plot_results: 輸出解析公式與蒙地卡羅模擬的對照表格與圖形
    
    fprintf('=== K=%d, %s, Set %d ===\n', K, mode, results.set_id);
    
    U = length(results.analytical);
    seq_names = results.seq_names;
    
    % 印出對齊的 Header
    fprintf('        ');
    for i = 1:U
        fprintf(' %-6s', seq_names{i});
    end
    fprintf('  Avg\n');
    
    % 取出數值並計算平均
    an_vals = results.analytical;
    mc_vals = results.mc;
    an_avg = mean(an_vals);
    mc_avg = mean(mc_vals);
    diff_vals = abs(an_vals - mc_vals);
    diff_avg = abs(an_avg - mc_avg);
    
    % 印出 Analytical Row
    fprintf('Analyt: ');
    for i = 1:U
        fprintf(' %-6.3f', an_vals(i));
    end
    fprintf('  %.3f\n', an_avg);
    
    % 印出 Monte Carlo Row
    fprintf('MC:     ');
    for i = 1:U
        fprintf(' %-6.3f', mc_vals(i));
    end
    fprintf('  %.3f\n', mc_avg);
    
    % 印出 Diff Row
    fprintf('Diff:   ');
    for i = 1:U
        fprintf(' %-6.3f', diff_vals(i));
    end
    fprintf('  %.3f\n\n', diff_avg);
    
    %% 繪製對照圖 (Analytical vs MC)
    figure;
    hold on;
    
    % x 軸為 1 到 U
    x = 1:U;
    
    % 畫 Analytial: 實線 + 正方形 marker
    plot(x, an_vals, '-s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Analytical');
    
    % 畫 Monte Carlo: 虛線 + 圓形 marker
    plot(x, mc_vals, '--o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Monte Carlo');
    
    title(sprintf('K=%d, Set %d AAoI (%s)', K, results.set_id, mode));
    
    % 設定 X 軸的 labels 為序列名稱
    xticks(x);
    xticklabels(seq_names);
    xlabel('User (Sequence)');
    ylabel('AAoI');
    
    legend('Location', 'best');
    grid on;
    hold off;
end
