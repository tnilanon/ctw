%% set verbosity
prSet(0); % 0 | 1 | 2 | ...

%% algorithm parameter
parDtw = [];
% parCca = st('d', .8, 'lams', .6); % CCA: reduce dimension to keep at least 0.8 energy, set the regularization weight to .6
% [gtw MAE] begin: 24.16; subtask_2: 31.45; subtask_3: 33.25; end: 30.85
% [gtw MAE] (<= 15 frames) begin: 44%; subtask_2: 33%; subtask_3: 31%; end: 31%
% [gtw MAE] (<= 30 frames) begin: 71%; subtask_2: 59%; subtask_3: 56%; end: 57%
parCca = st('d', .6, 'lams', .98);
parCtw = st('debg', 'n');

parGN = st('nItMa', 2, 'inp', 'linear'); % Gauss-Newton: 2 iterations to update the weight in GTW,
parGtw = st('nItMa', 20);

%% data
kinect_task_1_analysis_recording_list = {
    'P001', 'visit_1';
    'P001', 'visit_2';
    'P006', 'visit_1';
    'P006', 'visit_2';
    'P007', 'visit_1';
    'P009', 'visit_2';
    'P010', 'visit_1';
    'P010', 'visit_2';
    'P011', 'visit_1';
    'P011', 'visit_2';
    'P012', 'visit_1';
    'P012', 'visit_2';
    'P014', 'visit_1';
    'P014', 'visit_2';
    'P015', 'visit_1';
    'P015', 'visit_2';
    'P016', 'visit_1';
    'P017', 'visit_1';
    'P017', 'visit_2';
    'P018', 'visit_1';
    'P018', 'visit_2';
    'P019', 'visit_1';
    'P019', 'visit_2';
    'P023', 'visit_1';
    'P023', 'visit_2';
    'P024', 'visit_1';
    'P024', 'visit_2';
    'P026', 'visit_1';
    'P027', 'visit_1';
    'P027', 'visit_2';
    'P028', 'visit_1';
    'P028', 'visit_2';
    'P031', 'visit_1';
    'P034', 'visit_1';
    'P034', 'visit_2';
    'P035', 'visit_1';
    'P035', 'visit_2';
    'P036', 'visit_1';
    'P036', 'visit_2';
    'P037', 'visit_1';
    'P038', 'visit_1';
    'P038', 'visit_2';
    'P039', 'visit_1';
    'P040', 'visit_1';
    'P040', 'visit_2';
    'P042', 'visit_1';
    'P042', 'visit_2';
    'P043', 'visit_1';
    'P043', 'visit_2';
    'P044', 'visit_1';
    'P044', 'visit_2';
    'P045', 'visit_2';
    'P047', 'visit_1';
    'P048', 'visit_2';
    'P049', 'visit_1';
    'P049', 'visit_2';
    'P050', 'visit_1';
    'P050', 'visit_2';
    'P053', 'visit_1';
    'P053', 'visit_2';
    'P054', 'visit_1';
    'P054', 'visit_2';
    'P055', 'visit_1';
    'P055', 'visit_2';
};

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

n_recordings = size(kinect_task_1_analysis_recording_list, 1);
n_pairs = n_recordings * (n_recordings - 1) / 2;

pair_id_arr = cell(n_pairs, 1);
dtw_error_arr = zeros(n_pairs, 18);
ctw_error_arr = zeros(n_pairs, 18);
gtw_error_arr = zeros(n_pairs, 18);

%% for-looping over pairs of recordings
fprintf('total %4d: ', n_pairs);
m = 1;
for i = 1:n_recordings
    patient_id_a = kinect_task_1_analysis_recording_list{i, 1};
    visit_a = kinect_task_1_analysis_recording_list{i, 2};
    for j = i+1:n_recordings
        patient_id_b = kinect_task_1_analysis_recording_list{j, 1};
        visit_b = kinect_task_1_analysis_recording_list{j, 2};

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

        % dtw
        aliDtw = dtw(Xs, [], parDtw);

        % ctw
        aliCtw = ctw(Xs, aliDtw, [], parCtw, parCca, parDtw);

        % gtw
        aliGtw = gtw(Xs, bas, aliUtw, [], parGtw, parCca, parGN);

        % compute mae
        dtw_error_arr(m, :) = compute_alignment_error_v2(aliDtw, indice_a, indice_b, 1, 2);
        ctw_error_arr(m, :) = compute_alignment_error_v2(aliCtw, indice_a, indice_b, 1, 2);
        gtw_error_arr(m, :) = compute_alignment_error_interp_v2(aliGtw, indice_a, indice_b, 1, 2);

        if mod(m, 10) == 0
            fprintf('%4d, ', m);
        end
        m = m + 1;
%         break
    end
%     break
end
fprintf('done\n');

%%
% patient_id_a = 'P001';
% patient_id_b = 'P003';
% visit_a = 'visit_1';
% visit_b = 'visit_1';
% recording_a = struct('patient_id', patient_id_a, 'visit', visit_a);
% recording_b = struct('patient_id', patient_id_b, 'visit', visit_b);
% indice_a = get_analysis_indice(analysis_indice_task_1_T, recording_a.patient_id, recording_a.visit);
% indice_b = get_analysis_indice(analysis_indice_task_1_T, recording_b.patient_id, recording_b.visit);
% arr_a = get_data_array(recording_a, skeleton_analysis_feature_list);
% arr_b = get_data_array(recording_b, skeleton_analysis_feature_list);
% Xs = {arr_a, arr_b};
% D = conDst(Xs{1}, Xs{2});
% [v, S] = dtwFord(D);
% P = dtwBack(S);
% t.P = P;
% compute_alignment_mae(t, indice_a, indice_b)

%% save results
% save('2019-09-26-1800.mat', 'pair_id_arr', 'dtw_mae_arr', 'ctw_mae_arr');

%% get average MAE
temp = mean(dtw_error_arr, 1);
fprintf('\n[dtw MAE]\n');
fprintf('subtask_1: %.1f %.1f %.1f (%.3f %.3f %.3f);\nsubtask_2: %.1f %.1f %.1f (%.3f %.3f %.3f);\nsubtask_3: %.1f %.1f %.1f (%.3f %.3f %.3f);\n', ...
    temp);

temp = mean(ctw_error_arr, 1);
fprintf('\n[ctw MAE]\n');
fprintf('subtask_1: %.1f %.1f %.1f (%.3f %.3f %.3f);\nsubtask_2: %.1f %.1f %.1f (%.3f %.3f %.3f);\nsubtask_3: %.1f %.1f %.1f (%.3f %.3f %.3f);\n', ...
    temp);

temp = mean(gtw_error_arr, 1);
fprintf('\n[gtw MAE]\n');
fprintf('subtask_1: %.1f %.1f %.1f (%.3f %.3f %.3f);\nsubtask_2: %.1f %.1f %.1f (%.3f %.3f %.3f);\nsubtask_3: %.1f %.1f %.1f (%.3f %.3f %.3f);\n', ...
    temp);

%% show sequences
% shAlis2d({aliDtw, aliCtw}, 'legs', {'dtw', 'ctw'});
