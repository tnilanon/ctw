function mae = compute_alignment_mae_from_float(alignment, indice_a, indice_b, column_idx_a, column_idx_b)

path_a = alignment.P(:, column_idx_a);
path_b = alignment.P(:, column_idx_b);

mae_begin = mean([ ...
    abs(interp1(path_a, path_b, indice_a.index_begin, 'linear', 'extrap') - indice_b.index_begin), ...
    abs(interp1(path_b, path_a, indice_b.index_begin, 'linear', 'extrap') - indice_a.index_begin) ...
    ]);
mae_subtask_2 = mean([ ...
    abs(interp1(path_a, path_b, indice_a.index_subtask_2, 'linear', 'extrap') - indice_b.index_subtask_2), ...
    abs(interp1(path_b, path_a, indice_b.index_subtask_2, 'linear', 'extrap') - indice_a.index_subtask_2) ...
    ]);
mae_subtask_3 = mean([ ...
    abs(interp1(path_a, path_b, indice_a.index_subtask_3, 'linear', 'extrap') - indice_b.index_subtask_3), ...
    abs(interp1(path_b, path_a, indice_b.index_subtask_3, 'linear', 'extrap') - indice_a.index_subtask_3) ...
    ]);
mae_end = mean([ ...
    abs(interp1(path_a, path_b, indice_a.index_end, 'linear', 'extrap') - indice_b.index_end), ...
    abs(interp1(path_b, path_a, indice_b.index_end, 'linear', 'extrap') - indice_a.index_end) ...
    ]);

mae = struct( ...
    'begin', mae_begin, ...
    'subtask_2', mae_subtask_2, ...
    'subtask_3', mae_subtask_3, ...
    'end', mae_end ...
    );

