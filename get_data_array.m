function arr = get_data_array(recording_st, feature_list)

T = readtable(sprintf( ...
    '~/codes/kinect-skeleton-viewer-webapp/kinect_x_aligned_skeleton_files/%s_%s_TASK_1.csv', ...
    recording_st.patient_id, upper(recording_st.visit)));
arr = table2array(T(:, feature_list))';

