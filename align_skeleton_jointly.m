%% set verbosity
prSet(0); % 0 | 1 | 2 | ...

%% algorithm parameter
parDtw = [];
parPimw = st('lA', 1, 'lB', 1); % IMW: regularization weight
% parCca = st('d', 3, 'lams', .1); % CCA: reduce dimension to keep at least 0.95 energy
% parCca = st('d', .8, 'lams', .6); % CCA: reduce dimension to keep at least 0.8 energy, set the regularization weight to .6
% parCca = st('d', 3, 'lams', .1); % best
% parCca = st('d', .8, 'lams', .6); % worse
% parCca = st('d', .9, 'lams', .6); % worse
parPctw = [];
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

Xs = cell(1, n_patients * 2);
X_id = cell(1, n_patients * 2);
k = 1;
for i = 1:n_patients
    patient_id = kinect_task_1_analysis_patient_list{i};
    for j = 1:2
        visit = sprintf('visit_%1d', j);

        recording = struct('patient_id', patient_id, 'visit', visit);

        Xs{k} = get_data_array(recording, skeleton_analysis_feature_list);
        X_id{k} = sprintf('%s %s', patient_id, visit);

        k = k + 1;
    end
end

pair_id_arr = cell(n_pairs, 1);
pdtw_mae_arr = struct( ...
    'begin', cell(n_pairs, 1), ...
    'subtask_2', cell(n_pairs, 1), ...
    'subtask_3', cell(n_pairs, 1), ...
    'end', cell(n_pairs, 1));
pddtw_mae_arr = struct( ...
    'begin', cell(n_pairs, 1), ...
    'subtask_2', cell(n_pairs, 1), ...
    'subtask_3', cell(n_pairs, 1), ...
    'end', cell(n_pairs, 1));
pimw_mae_arr = struct( ...
    'begin', cell(n_pairs, 1), ...
    'subtask_2', cell(n_pairs, 1), ...
    'subtask_3', cell(n_pairs, 1), ...
    'end', cell(n_pairs, 1));
pctw_mae_arr = struct( ...
    'begin', cell(n_pairs, 1), ...
    'subtask_2', cell(n_pairs, 1), ...
    'subtask_3', cell(n_pairs, 1), ...
    'end', cell(n_pairs, 1));
gtw_mae_arr = struct( ...
    'begin', cell(n_pairs, 1), ...
    'subtask_2', cell(n_pairs, 1), ...
    'subtask_3', cell(n_pairs, 1), ...
    'end', cell(n_pairs, 1));

%% monotonic basis
ns = cellDim(Xs, 2);
len = round(max(ns) * 1.1);
bas = baTems(len, ns, 'pol', [5 .5], 'tan', [5 1 1]); % 5 polynomial and 5 tangent functions

%% optimize jointly
% utw (initialization)
fprintf('started utw @ %s\n', datetime);
aliUtw = utw(Xs, bas, []);
fprintf('done with utw @ %s\n', datetime);

% pdtw
fprintf('started pdtw @ %s\n', datetime);
aliPdtw = pdtw(Xs, aliUtw, [], parDtw);
fprintf('done with pdtw @ %s\n', datetime);

% pddtw
fprintf('started pddtw @ %s\n', datetime);
aliPddtw = pddtw(Xs, aliUtw, [], parDtw);
fprintf('done with pddtw @ %s\n', datetime);

% pimw
fprintf('started pimw @ %s\n', datetime);
aliPimw = pimw(Xs, aliUtw, [], parPimw, parDtw);
fprintf('done with pimw @ %s\n', datetime);
%%
% pctw
fprintf('started pctw @ %s\n', datetime);
aliPctw = pctw(Xs, aliUtw, [], parPctw, parCca, parDtw);
fprintf('done with pctw @ %s\n', datetime);

% gtw
fprintf('started gtw @ %s\n', datetime);
aliGtw = gtw(Xs, bas, aliUtw, [], parGtw, parCca, parGN);
fprintf('done with gtw @ %s\n', datetime);

%% save intermediate results

% started utw @ 27-Sep-2019 06:01:15
% done with utw @ 27-Sep-2019 06:01:15
% started pdtw @ 27-Sep-2019 06:01:15
% done with pdtw @ 27-Sep-2019 06:02:25
% started pddtw @ 27-Sep-2019 06:02:25
% done with pddtw @ 27-Sep-2019 06:02:34
% started pimw @ 27-Sep-2019 06:02:34
% done with pimw @ 27-Sep-2019 06:09:18
% started pctw @ 27-Sep-2019 07:34:44
% done with pctw @ 27-Sep-2019 07:37:22
% started gtw @ 27-Sep-2019 07:37:22
% done with gtw @ 27-Sep-2019 07:38:43

% save('align_skeleton_jointly.warping_path.2019-09-27-0613.mat', ...
%     'aliUtw', 'aliPdtw', 'aliPddtw', 'aliPimw', 'aliPctw', 'aliGtw');

%% debugging
% patient_id_a = 'P001';
% visit_a = 'visit_1';
% 
% patient_id_b = 'P005';
% visit_b = 'visit_2';
% 
% recording_a = struct('patient_id', patient_id_a, 'visit', visit_a);
% recording_b = struct('patient_id', patient_id_b, 'visit', visit_b);
% 
% indice_a = get_analysis_indice(analysis_indice_task_1_T, recording_a.patient_id, recording_a.visit);
% indice_b = get_analysis_indice(analysis_indice_task_1_T, recording_b.patient_id, recording_b.visit);
% 
% column_idx_a = 1;
% column_idx_b = 2;
% 
% alignment = aliGtw;
% 
% path_a = alignment.P(:, column_idx_a);
% path_b = alignment.P(:, column_idx_b);

%% for-looping over pairs of recordings for evaluation
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

                % compute mae
                pdtw_mae_arr(m) = compute_alignment_mae(aliPdtw, indice_a, indice_b, 2*(i-1)+k, 2*(j-1)+l);
                pddtw_mae_arr(m) = compute_alignment_mae(aliPddtw, indice_a, indice_b, 2*(i-1)+k, 2*(j-1)+l);
                pimw_mae_arr(m) = compute_alignment_mae(aliPimw, indice_a, indice_b, 2*(i-1)+k, 2*(j-1)+l);
                pctw_mae_arr(m) = compute_alignment_mae(aliPctw, indice_a, indice_b, 2*(i-1)+k, 2*(j-1)+l);
                gtw_mae_arr(m) = compute_alignment_mae_from_float(aliGtw, indice_a, indice_b, 2*(i-1)+k, 2*(j-1)+l);

                if mod(m, 100) == 0
                    fprintf('%4d/%4d\n', m, n_pairs);
                end
                m = m + 1;
            end
        end
    end
end
fprintf('done\n');

%% save results
% save('align_skeleton_jointly.results.2019-09-26-1800.mat', 'pair_id_arr', ...
%     'pdtw_mae_arr', 'pddtw_mae_arr', 'pimw_mae_arr', 'pctw_mae_arr', 'gtw_mae_arr');

%% get average MAE
fprintf('\n');
fprintf('[pdtw MAE] begin: %.2f; subtask_2: %.2f; subtask_3: %.2f; end: %.2f\n', ...
    mean([pdtw_mae_arr.begin]), mean([pdtw_mae_arr.subtask_2]), mean([pdtw_mae_arr.subtask_3]), mean([pdtw_mae_arr.end]));
fprintf('[pdtw MAE] '); eval_pct_mae_within_x_frames(pdtw_mae_arr, 15);
fprintf('[pdtw MAE] '); eval_pct_mae_within_x_frames(pdtw_mae_arr, 30);

fprintf('\n');
fprintf('[pddtw MAE] begin: %.2f; subtask_2: %.2f; subtask_3: %.2f; end: %.2f\n', ...
    mean([pddtw_mae_arr.begin]), mean([pddtw_mae_arr.subtask_2]), mean([pddtw_mae_arr.subtask_3]), mean([pddtw_mae_arr.end]));
fprintf('[pddtw MAE] '); eval_pct_mae_within_x_frames(pddtw_mae_arr, 15);
fprintf('[pddtw MAE] '); eval_pct_mae_within_x_frames(pddtw_mae_arr, 30);

fprintf('\n');
fprintf('[pimw MAE] begin: %.2f; subtask_2: %.2f; subtask_3: %.2f; end: %.2f\n', ...
    mean([pimw_mae_arr.begin]), mean([pimw_mae_arr.subtask_2]), mean([pimw_mae_arr.subtask_3]), mean([pimw_mae_arr.end]));
fprintf('[pimw MAE] '); eval_pct_mae_within_x_frames(pimw_mae_arr, 15);
fprintf('[pimw MAE] '); eval_pct_mae_within_x_frames(pimw_mae_arr, 30);

fprintf('\n');
fprintf('[pctw MAE] begin: %.2f; subtask_2: %.2f; subtask_3: %.2f; end: %.2f\n', ...
    mean([pctw_mae_arr.begin]), mean([pctw_mae_arr.subtask_2]), mean([pctw_mae_arr.subtask_3]), mean([pctw_mae_arr.end]));
fprintf('[pctw MAE] '); eval_pct_mae_within_x_frames(pctw_mae_arr, 15);
fprintf('[pctw MAE] '); eval_pct_mae_within_x_frames(pctw_mae_arr, 30);

fprintf('\n');
fprintf('[gtw MAE] begin: %.2f; subtask_2: %.2f; subtask_3: %.2f; end: %.2f\n', ...
    mean([gtw_mae_arr.begin]), mean([gtw_mae_arr.subtask_2]), mean([gtw_mae_arr.subtask_3]), mean([gtw_mae_arr.end]));
fprintf('[gtw MAE] '); eval_pct_mae_within_x_frames(gtw_mae_arr, 15);
fprintf('[gtw MAE] '); eval_pct_mae_within_x_frames(gtw_mae_arr, 30);
