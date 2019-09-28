%% set verbosity
prSet(0); % 0 | 1 | 2 | ...

%% grid search parameter
par_cca_d_list = [0.8, 0.85, 0.9, 0.95];
par_cca_lambda_list = 0.1:0.1:0.9;

%% algorithm parameter
% parDtw = [];
% parCca = st('d', .8, 'lams', .6); % CCA: reduce dimension to keep at least 0.8 energy, set the regularization weight to .6
% parCtw = st('debg', 'n');

% parCca will be instantiate in the loop
parGN = st('nItMa', 2, 'inp', 'linear'); % Gauss-Newton: 2 iterations to update the weight in GTW,
parGtw = st('nItMa', 20);

%% data
kinect_task_1_analysis_patient_list = { ...
    'P001'; 'P003'; 'P005'; 'P006'; 'P007'; 'P009'; 'P010'; 'P011'; ...
    'P012'; 'P014'; 'P015'; 'P016'; 'P017'; 'P018'; 'P019'; 'P023'; ...
    'P024'; 'P026'; 'P027'; 'P028'; 'P031'; 'P034'; 'P035'; 'P036'; ...
    'P037'; 'P038'; 'P039'; 'P040'; 'P041'; 'P042'; 'P043'; 'P044'; ...
    'P045'; 'P047'; 'P048'; 'P049'; 'P050'; 'P053'; 'P054'; 'P055'};

analysis_indice_task_1_T = readtable('analysis_indice_task_1.csv');
analysis_indice_task_1_T.index_begin = analysis_indice_task_1_T.index_begin + 1;
analysis_indice_task_1_T.index_subtask_2 = analysis_indice_task_1_T.index_subtask_2 + 1;
analysis_indice_task_1_T.index_subtask_3 = analysis_indice_task_1_T.index_subtask_3 + 1;
analysis_indice_task_1_T.index_end = analysis_indice_task_1_T.index_end + 1;

skeleton_analysis_feature_list = { ...
    'spine_base_x'; 'spine_base_y'; 'spine_base_z'; ...
    'spine_mid_x'; 'spine_mid_y'; 'spine_mid_z'; ...
    'spine_shoulder_x'; 'spine_shoulder_y'; 'spine_shoulder_z'; ...
    'neck_x'; 'neck_y'; 'neck_z'; ...
    'head_x'; 'head_y'; 'head_z'; ...
    'shoulder_left_x'; 'shoulder_left_y'; 'shoulder_left_z'; ...
    'elbow_left_x'; 'elbow_left_y'; 'elbow_left_z'; ...
    'shoulder_right_x'; 'shoulder_right_y'; 'shoulder_right_z'; ...
    'elbow_right_x'; 'elbow_right_y'; 'elbow_right_z'; ...
    'hip_left_x'; 'hip_left_y'; 'hip_left_z'; ...
    'knee_left_x'; 'knee_left_y'; 'knee_left_z'; ...
    'hip_right_x'; 'hip_right_y'; 'hip_right_z'; ...
    'knee_right_x'; 'knee_right_y'; 'knee_right_z' ...
    };

n_patients = numel(kinect_task_1_analysis_patient_list);
n_pairs = n_patients * (n_patients - 1) * 2;

pair_id_arr = cell(n_pairs, 1);
gtw_mae_arr = struct( ...
    'begin', cell(n_pairs, 1), ...
    'subtask_2', cell(n_pairs, 1), ...
    'subtask_3', cell(n_pairs, 1), ...
    'end', cell(n_pairs, 1));

n_grid_searches = numel(par_cca_d_list) * numel(par_cca_lambda_list);

grid_search_result = struct( ...
    'd', cell(n_grid_searches, 1), ...
    'lambda', cell(n_grid_searches, 1), ...
    'mae_begin', cell(n_grid_searches, 1), ...
    'mae_subtask_2', cell(n_grid_searches, 1), ...
    'mae_subtask_3', cell(n_grid_searches, 1), ...
    'mae_end', cell(n_grid_searches, 1), ...
    'pct_15f_begin', cell(n_grid_searches, 1), ...
    'pct_15f_subtask_2', cell(n_grid_searches, 1), ...
    'pct_15f_subtask_3', cell(n_grid_searches, 1), ...
    'pct_15f_end', cell(n_grid_searches, 1), ...
    'pct_30f_begin', cell(n_grid_searches, 1), ...
    'pct_30f_subtask_2', cell(n_grid_searches, 1), ...
    'pct_30f_subtask_3', cell(n_grid_searches, 1), ...
    'pct_30f_end', cell(n_grid_searches, 1));

%% grid search
idx_result = 1;
for par_cca_d = par_cca_d_list
    for par_cca_lambda = par_cca_lambda_list
        parCca = st('d', par_cca_d, 'lams', par_cca_lambda);
        fprintf('d: %.2f; lambda: %.1f;\n', par_cca_d, par_cca_lambda);

        % for-looping over pairs of recordings
        fprintf('total %4d: ', n_pairs);
        m = 1;
        for i = 1:n_patients
            patient_id_a = kinect_task_1_analysis_patient_list{i};
            for j = i+1:n_patients
                patient_id_b = kinect_task_1_analysis_patient_list{j};
                for k = 1:2
                    visit_a = sprintf('visit_%1d', k);
                    for l = 1:2
                        visit_b = sprintf('visit_%1d', l);

                        pair_id_arr(m) = {sprintf('%s %s, %s %s', patient_id_a, visit_a, patient_id_b, visit_b)};

                        recording_a = struct('patient_id', patient_id_a, 'visit', visit_a);
                        recording_b = struct('patient_id', patient_id_b, 'visit', visit_b);

                        indice_a = get_analysis_indice(analysis_indice_task_1_T, recording_a.patient_id, recording_a.visit);
                        indice_b = get_analysis_indice(analysis_indice_task_1_T, recording_b.patient_id, recording_b.visit);

                        arr_a = get_data_array(recording_a, skeleton_analysis_feature_list);
                        arr_b = get_data_array(recording_b, skeleton_analysis_feature_list);

                        Xs = {arr_a, arr_b};

                        % monotonic basis
                        ns = cellDim(Xs, 2);
                        len = round(max(ns) * 1.1);
                        bas = baTems(len, ns, 'pol', [5 .5], 'tan', [5 1 1]); % 5 polynomial and 5 tangent functions

                        % utw (initialization for Procrustes and gtw)
                        aliUtw = utw(Xs, bas, []);

        %                 % dtw
        %                 aliDtw = dtw(Xs, [], parDtw);
        % 
        %                 % ctw
        %                 aliCtw = ctw(Xs, aliDtw, [], parCtw, parCca, parDtw);

                        % gtw
                        aliGtw = gtw(Xs, bas, aliUtw, [], parGtw, parCca, parGN);

                        % compute mae
        %                 dtw_mae_arr(m) = compute_alignment_mae(aliDtw, indice_a, indice_b, 1, 2);
        %                 ctw_mae_arr(m) = compute_alignment_mae(aliCtw, indice_a, indice_b, 1, 2);
                        gtw_mae_arr(m) = compute_alignment_mae_from_float(aliGtw, indice_a, indice_b, 1, 2);

                        if mod(m, 100) == 0
                            fprintf('%4d, ', m);
                        end
                        m = m + 1;
                    end
                end
            end
        end
        fprintf('done\n');
        
        grid_search_result(idx_result).mae_begin = mean([gtw_mae_arr.begin]);
        grid_search_result(idx_result).mae_subtask_2 = mean([gtw_mae_arr.subtask_2]);
        grid_search_result(idx_result).mae_subtask_3 = mean([gtw_mae_arr.subtask_3]);
        grid_search_result(idx_result).mae_end = mean([gtw_mae_arr.end]);
        
        [grid_search_result(idx_result).pct_15f_begin, ...
            grid_search_result(idx_result).pct_15f_subtask_2, ...
            grid_search_result(idx_result).pct_15f_subtask_3, ...
            grid_search_result(idx_result).pct_15f_end ...
            ] = feval(@(x) x{:}, num2cell(eval_pct_mae_within_x_frames(gtw_mae_arr, 15)));
        
        [grid_search_result(idx_result).pct_30f_begin, ...
            grid_search_result(idx_result).pct_30f_subtask_2, ...
            grid_search_result(idx_result).pct_30f_subtask_3, ...
            grid_search_result(idx_result).pct_30f_end ...
            ] = feval(@(x) x{:}, num2cell(eval_pct_mae_within_x_frames(gtw_mae_arr, 30)));
        
        idx_result = idx_result + 1;
    end
end

