function plot_results(results, K, mode)
    fprintf('=== K=%d, %s ===\n', K, mode);
    
    if K == 6
        % Labels for Set 1 and Set 2 based on spec
        labels1 = {'s0', 's2', 's3', 's4', 's5', 's7', 'Avg'};
        labels2 = {'s1', 's2', 's3', 's4', 's5', 's6', 'Avg'};
        
        % Print header
        fprintf('         s0      s2      s3      s4      s5      s7    Avg\n');
        
        vals1 = results.set1;
        avg1 = mean(vals1);
        fprintf('Set 1:  %.3f   %.3f   %.3f   %.3f   %.3f   %.3f  %.3f\n', vals1(1), vals1(2), vals1(3), vals1(4), vals1(5), vals1(6), avg1);
        
        fprintf('         s1      s2      s3      s4      s5      s6    Avg\n');
        vals2 = results.set2;
        avg2 = mean(vals2);
        fprintf('Set 2:  %.3f   %.3f   %.3f   %.3f   %.3f   %.3f  %.3f\n', vals2(1), vals2(2), vals2(3), vals2(4), vals2(5), vals2(6), avg2);
        
        % Plotting
        figure;
        bar([vals1', vals2']);
        title(sprintf('K=%d AAoI Comparison (%s)', K, mode));
        xlabel('User index in Set');
        ylabel('AAoI');
        legend('Set 1', 'Set 2');
        grid on;
    else
        % For K=5
        fprintf('         seq1    seq2    seq3    seq4    seq5    Avg\n');
        vals1 = results.set1;
        avg1 = mean(vals1);
        fprintf('Set 1:  %.3f   %.3f   %.3f   %.3f   %.3f  %.3f\n', vals1(1), vals1(2), vals1(3), vals1(4), vals1(5), avg1);
        
        figure;
        bar(vals1');
        title(sprintf('K=%d AAoI Comparison (%s)', K, mode));
        xlabel('User index in Set');
        ylabel('AAoI');
        legend('Set 1');
        grid on;
    end
    fprintf('\n');
end
