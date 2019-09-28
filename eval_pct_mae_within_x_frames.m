function [pct_begin, pct_subtask_2, pct_subtask_3, pct_end] = eval_pct_mae_within_x_frames(mae_arr, x)

pct_begin = mean([mae_arr.begin] <= x) * 100;
pct_subtask_2 = mean([mae_arr.subtask_2] <= x) * 100;
pct_subtask_3 = mean([mae_arr.subtask_3] <= x) * 100;
pct_end = mean([mae_arr.end] <= x) * 100;

fprintf('(<= %d frames) begin: %.0f%%; subtask_2: %.0f%%; subtask_3: %.0f%%; end: %.0f%%\n', ...
    x, pct_begin, pct_subtask_2, pct_subtask_3, pct_end);

