%% set verbosity
prSet(0); % 0 | 1 | 2 | ...

%% algorithm parameter
parDtw = [];
parCca = st('d', .8, 'lams', .6); % CCA: reduce dimension to keep at least 0.8 energy, set the regularization weight to .6
parCtw = st('debg', 'n');

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
dtw_mae_arr = struct( ...
    'begin', cell(n_pairs, 1), ...
    'subtask_2', cell(n_pairs, 1), ...
    'subtask_3', cell(n_pairs, 1), ...
    'end', cell(n_pairs, 1));
ctw_mae_arr = struct( ...
    'begin', cell(n_pairs, 1), ...
    'subtask_2', cell(n_pairs, 1), ...
    'subtask_3', cell(n_pairs, 1), ...
    'end', cell(n_pairs, 1));

%% for-looping over pairs of recordings
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

                % utw (initialization)
                aliUtw = utw(Xs, [], []);

                % dtw
                aliDtw = dtw(Xs, [], parDtw);

                % ctw
                aliCtw = ctw(Xs, aliDtw, [], parCtw, parCca, parDtw);

                % compute mae
                dtw_mae_arr(m) = compute_alignment_mae(aliDtw, indice_a, indice_b, 1, 2);
                ctw_mae_arr(m) = compute_alignment_mae(aliCtw, indice_a, indice_b, 1, 2);

                if mod(m, 10) == 0
                    fprintf('%4d/%4d\n', m, n_pairs);
                end
                m = m + 1;
            end
        end
    end
end
fprintf('done\n');

%%
patient_id_a = 'P001';
patient_id_b = 'P003';
visit_a = 'visit_1';
visit_b = 'visit_1';
recording_a = struct('patient_id', patient_id_a, 'visit', visit_a);
recording_b = struct('patient_id', patient_id_b, 'visit', visit_b);
indice_a = get_analysis_indice(analysis_indice_task_1_T, recording_a.patient_id, recording_a.visit);
indice_b = get_analysis_indice(analysis_indice_task_1_T, recording_b.patient_id, recording_b.visit);
arr_a = get_data_array(recording_a, skeleton_analysis_feature_list);
arr_b = get_data_array(recording_b, skeleton_analysis_feature_list);
Xs = {arr_a, arr_b};
D = conDst(Xs{1}, Xs{2});
[v, S] = dtwFord(D);
P = dtwBack(S);
t.P = P;
compute_alignment_mae(t, indice_a, indice_b)

%% save results
% save('2019-09-26-1800.mat', 'pair_id_arr', 'dtw_mae_arr', 'ctw_mae_arr');

%% get average MAE
fprintf('\n');
fprintf('[dtw MAE] begin: %.2f; subtask_2: %.2f; subtask_3: %.2f; end: %.2f\n', ...
    mean([dtw_mae_arr.begin]), mean([dtw_mae_arr.subtask_2]), mean([dtw_mae_arr.subtask_3]), mean([dtw_mae_arr.end]));
fprintf('[dtw MAE] '); eval_pct_mae_within_x_frames(dtw_mae_arr, 15);
fprintf('[dtw MAE] '); eval_pct_mae_within_x_frames(dtw_mae_arr, 30);

fprintf('\n');
fprintf('[ctw MAE] begin: %.2f; subtask_2: %.2f; subtask_3: %.2f; end: %.2f\n', ...
    mean([ctw_mae_arr.begin]), mean([ctw_mae_arr.subtask_2]), mean([ctw_mae_arr.subtask_3]), mean([ctw_mae_arr.end]));
fprintf('[ctw MAE] '); eval_pct_mae_within_x_frames(ctw_mae_arr, 15);
fprintf('[ctw MAE] '); eval_pct_mae_within_x_frames(ctw_mae_arr, 30);

%% show sequences
% shAlis2d({aliDtw, aliCtw}, 'legs', {'dtw', 'ctw'});
