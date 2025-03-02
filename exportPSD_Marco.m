 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Brainsense Survey Analysis over time for all channels %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input
% a folder with .mat files containing brainsense survey data
% (LMTD = LFP montage time domain)
% .mat files must have been extracted before from medtronic json files using the
% percieve tool by Julian Neumann

% Operations:
% Load LFP data files from a folder, calculate PSD, extract beta peaks, extract
% frequency corresponding to beta peaks

% Outputs:
% line graph with peak beta amplitude for all channels over time
% PSDs for all channels with peak amplitudes highlighted

%% marco.treven@meduniwien.ac.at; v1.0

% Prompt a window to select a folder that contains all the '.mat' files
folder_path = uigetdir;

% Find all .mat files in the folder and its subfolders
mat_files = dir(fullfile(folder_path, '**/*.mat')); % '**/' is used to include subfolders

% Extract file names
mat_files = arrayfun(@(x) fullfile(x.folder, x.name), mat_files, 'UniformOutput', false);

% Filter file names containing "LMTD"
mat_files_LMTD = mat_files(cellfun(@(x) contains(x, 'LMTD'), mat_files));

% Remove file extensions (.mat) from the filtered file names
mat_files_noext = cellfun(@(x) x(1:end-4), mat_files_LMTD, 'UniformOutput', false);

% Extract datetime information from file names
datetime_str = cellfun(@(x) x(end-15:end-2), mat_files_noext, 'UniformOutput', false);
dates = datetime(datetime_str, 'InputFormat', 'yyyyMMddHHmmss');

% Number of files
num_files = length(mat_files_LMTD);

% Initialize channel_labels variable
channel_labels = {'LFP_Stn_0_1_LEFT_RING', 'LFP_Stn_0_2_LEFT_RING', 'LFP_Stn_0_3_LEFT_RING', ...
                  'LFP_Stn_1_2_LEFT_RING', 'LFP_Stn_1_3_LEFT_RING', 'LFP_Stn_2_3_LEFT_RING', ...
                  'LFP_Stn_0_1_RIGHT_RING', 'LFP_Stn_0_2_RIGHT_RING', 'LFP_Stn_0_3_RIGHT_RING', ...
                  'LFP_Stn_1_2_RIGHT_RING', 'LFP_Stn_1_3_RIGHT_RING', 'LFP_Stn_2_3_RIGHT_RING', ...
                  'LFP_Stn_1_A_1_B_LEFT_SEGMENT', 'LFP_Stn_1_B_1_C_LEFT_SEGMENT', 'LFP_Stn_1_A_1_C_LEFT_SEGMENT', ...
                  'LFP_Stn_2_A_2_B_LEFT_SEGMENT', 'LFP_Stn_2_B_2_C_LEFT_SEGMENT', 'LFP_Stn_2_A_2_C_LEFT_SEGMENT', ...
                  'LFP_Stn_1_A_1_B_RIGHT_SEGMENT', 'LFP_Stn_1_B_1_C_RIGHT_SEGMENT', 'LFP_Stn_1_A_1_C_RIGHT_SEGMENT', ...
                  'LFP_Stn_2_A_2_B_RIGHT_SEGMENT', 'LFP_Stn_2_B_2_C_RIGHT_SEGMENT', 'LFP_Stn_2_A_2_C_RIGHT_SEGMENT', ...
                  'LFP_Stn_1_A_2_A_LEFT_SEGMENT', 'LFP_Stn_1_B_2_B_LEFT_SEGMENT', 'LFP_Stn_1_C_2_C_LEFT_SEGMENT', ...
                  'LFP_Stn_1_A_2_A_RIGHT_SEGMENT', 'LFP_Stn_1_B_2_B_RIGHT_SEGMENT', 'LFP_Stn_1_C_2_C_RIGHT_SEGMENT'};

% Total number of channels
total_channels = length(channel_labels);

% Preallocate the channel_side_all matrix
channel_side_all = nan(total_channels, num_files);
% Preallocate the psd_results_all and smoothed_psd_results_all cell arrays
psd_results_all = cell(total_channels, num_files);
smoothed_psd_results_all = cell(total_channels, num_files);
detrended_psd_results_all = cell(total_channels, num_files);

%% Iterating through each file to fill channel_side_all and calculate PSD
for file_idx = 1:num_files
    % Load LFP data
    lfpdata = load(mat_files_LMTD{file_idx});

    % Get the current number of channels from the loaded data
    current_num_channels = size(lfpdata.data.label, 1);
    
    % Find the indices of the current channels in the total channels list
    current_indices = find(ismember(channel_labels, lfpdata.data.label));

    % Fill channel_side_all with 1=left and 2=right
    for j = 1:current_num_channels
        if contains(lfpdata.data.label{j}, '_LEFT_')
            channel_side_all(current_indices(j), file_idx) = 1; % 'L'
        elseif contains(lfpdata.data.label{j}, '_RIGHT_')
            channel_side_all(current_indices(j), file_idx) = 2; % 'R'
        end
    end

    % Get the LFP data for each channel (6x5288 double)
    lfp_data = lfpdata.data.trial{1};

    % Get the sampling frequency
    fs = lfpdata.data.fsample;

    % Compute the PSD for each channel and store the results
    for channel_idx = 1:size(lfp_data, 1)
        % Get the current LFP data
        current_lfp = lfp_data(channel_idx, :);
        
        % Apply a high-pass filter at 5 Hz
        [b_high, a_high] = butter(4, 5 / (fs / 2), 'high');
        current_lfp_highpass = filtfilt(b_high, a_high, current_lfp);
        
        % Apply a low-pass filter at 90 Hz
        [b_low, a_low] = butter(4, 90 / (fs / 2), 'low');
        current_lfp_filtered = filtfilt(b_low, a_low, current_lfp_highpass);
        
        % Calculate PSD for the filtered LFP data
        [pxx, f] = pwelch(current_lfp_filtered, [], [], [], fs);  

        % Save the PSD results for the current channel in the corresponding
        % row and column of the psd_results_all cell array
        psd_results_all{current_indices(channel_idx), file_idx} = pxx;

        % Smooth the PSD data with a window size of x
        pxx_smoothed = smooth(pxx, 15);

        % Save the smoothed PSD results for the current channel in the
        % corresponding row and column of the smoothed_psd_results_all cell array
        smoothed_psd_results_all{current_indices(channel_idx), file_idx} = pxx_smoothed;

        % Detrend the log-transformed PSD data
        log_pxx = log10(pxx_smoothed);
        poly_coeffs = polyfit(f, log_pxx, 1);
        fitted_values = polyval(poly_coeffs, f);

        % Subtract the fitted values from the original log-transformed PSD data to obtain the detrended PSD data
        log_pxx_detrended = log_pxx - fitted_values;

        % Transform the detrended data back to its original scale by taking the exponent of the detrended data
        pxx_detrended = 10 .^ log_pxx_detrended;
        
        % Save the detrended PSD results for the current channel in the
        % corresponding row and column of the detrended_psd_results_all cell array
        detrended_psd_results_all{current_indices(channel_idx), file_idx} = pxx_detrended;

    end
end


%% Extract peak values, peak frequencies, and peak prominence from the smoothed PSD data in the cell array
%% Determine those betamax values that are above noise

% Preallocate matrices for peak values, peak frequencies, and peak prominences
peak_value_all = nan(total_channels, num_files);
peak_frequency_all = nan(total_channels, num_files);
peak_prominence_all = nan(total_channels, num_files);
significant_peak_all = nan(total_channels, num_files);

% Define frequency range of interest
f_min = 10;
f_max = 35;

% Find the indices of the frequency range of interest
idx_f_min = find(f >= f_min, 1, 'first');
idx_f_max = find(f <= f_max, 1, 'last');

% Iterate through each file and channel
for file_idx = 1:num_files
    for channel_idx = 1:total_channels
        % Check if the cell contains PSD data
        if ~isempty(detrended_psd_results_all{channel_idx, file_idx})
            % Extract PSD data in the frequency range of interest
            pxx_detrended = detrended_psd_results_all{channel_idx, file_idx}(idx_f_min:idx_f_max);
            freq_range = f(idx_f_min:idx_f_max);
            

            % Find the peak value, peak frequency, and peak prominence
            [peak_value, peak_loc, ~, peak_prominence] = findpeaks(pxx_detrended, freq_range, 'NPeaks', 1, 'SortStr', 'descend');

            % Store the results in the corresponding matrices
            if ~isempty(peak_value)
                peak_value_all(channel_idx, file_idx) = peak_value;
                peak_frequency_all(channel_idx, file_idx) = peak_loc;
                peak_prominence_all(channel_idx, file_idx) = peak_prominence;

                % Calculate the local mean and standard deviation based on the frequency range of interest
                local_mean = mean(pxx_detrended);
                local_std = std(pxx_detrended);

                % Calculate the z-score
                z_score = (peak_value - local_mean) / local_std;

                % Check if the z-score is higher than the threshold
                z_threshold = 2;
                if z_score > z_threshold
                   % The peak is significantly higher than the surrounding average signal
                   significant_peak_all(channel_idx, file_idx) = 1;
                   else
                   significant_peak_all(channel_idx, file_idx) = 0;
                end
            end
        end
    end
end

%%%%%%%%%%%%%%
%% PLOTTING %%
%%%%%%%%%%%%%%

%% Plot all the peak values and frequencies in a single figure, with significant peaks in green
figure; % Create a new figure

% Iterate through each file and channel
for file_idx = 1:num_files
    for channel_idx = 1:total_channels
        % Check if there is data in the corresponding fields
        if ~isnan(peak_value_all(channel_idx, file_idx)) && ~isnan(peak_frequency_all(channel_idx, file_idx))
            % Set the color based on the value in significant_peak_all
            if significant_peak_all(channel_idx, file_idx) == 1
                point_color = 'g'; % Green
            else
                point_color = 'r'; % Red
            end

            % Plot the point with the appropriate color
            scatter(peak_frequency_all(channel_idx, file_idx), peak_value_all(channel_idx, file_idx), point_color, 'filled');
            hold on;
        end
    end
end

xlabel('Peak Frequency (Hz)');
ylabel('Peak Value');
title('Peak Value vs. Peak Frequency');
set(gca, 'YScale', 'log'); % Set the y-axis to logarithmic scale
grid on; % Enable the grid

%% 3D scatter plot with the peak prominence on the z-axis
figure; % Create a new figure

% Iterate through each file and channel
for file_idx = 1:num_files
    for channel_idx = 1:total_channels
        % Check if there is data in the corresponding fields
        if ~isnan(peak_value_all(channel_idx, file_idx)) && ~isnan(peak_frequency_all(channel_idx, file_idx)) && ~isnan(peak_prominence_all(channel_idx, file_idx))
            % Set the color based on the value in significant_peak_all
            if significant_peak_all(channel_idx, file_idx) == 1
                point_color = 'g'; % Green
            else
                point_color = 'r'; % Red
            end

            % Plot the point with the appropriate color
            scatter3(peak_frequency_all(channel_idx, file_idx), peak_value_all(channel_idx, file_idx), peak_prominence_all(channel_idx, file_idx), point_color, 'filled');
            hold on;
        end
    end
end

xlabel('Peak Frequency (Hz)');
ylabel('Peak Value');
zlabel('Peak Prominence');
title('Peak Value vs. Peak Frequency vs. Peak Prominence');
set(gca, 'YScale', 'log');
set(gca, 'ZScale', 'log');
grid on; % Enable the grid



%% Plot all the channels per file in the same schema
% Set the number of subplot rows and columns
subplot_rows = 5;
subplot_columns = 6;

% Define the size of the figure
figure_width = 1800; % pixels
figure_height = 1500; % pixels

% Loop through all files
for file_idx = 1:num_files
    % Create a new figure for the current file
    fig = figure;
    
    % Set the figure size and position
    fig.Position = [100, 100, figure_width, figure_height];
    
    % Set the figure title to the datetime of the file
    sgtitle(datestr(dates(file_idx), 'yyyy-mm-dd HH:MM:SS'));
    
    % Loop through all channels
    for channel_idx = 1:total_channels
        % Check if the current cell contains data
        if ~isempty(detrended_psd_results_all{channel_idx, file_idx})
            % Get the smoothed PSD data and corresponding frequencies
            pxx_smoothed = detrended_psd_results_all{channel_idx, file_idx};
            f_data = f;
            
            % Create a subplot for the current channel
            subplot(subplot_rows, subplot_columns, channel_idx);
            
            % Plot the PSD curve on a log y axis
            %semilogy(f_data, pxx_smoothed);
            plot(f_data, pxx_smoothed);
            hold on;
            
            % Get the betamax and f_max values for the current channel and file
            betamax = peak_value_all(channel_idx, file_idx);
            f_max = peak_frequency_all(channel_idx, file_idx);
            significant_peak = significant_peak_all(channel_idx, file_idx);
            
            % Choose the circle color based on the significant_peak value
            if significant_peak == 1
                circle_color = 'g';
            else
                circle_color = 'r';
            end
            
            % Plot the circle at the (f_max, betamax) point
            scatter(f_max, betamax, 50, circle_color, 'o', 'filled');
            
            % Set the subplot title to the channel label
            title(channel_labels{channel_idx});
            
            % Set the x-axis and y-axis labels
            xlabel('Frequency (Hz)');
            ylabel('µV^2/Hz');
            
            % Display the betamax and f_max values in the subplot
            text(0.65, 0.9, sprintf('Betamax: %.2f', betamax), 'Units', 'normalized');
            text(0.65, 0.8, sprintf('F_max: %.2f Hz', f_max), 'Units', 'normalized');
            
            % Set the x-axis limits
            xlim([0 100]);

            % log y axis
            set(gca, 'YScale', 'log');
            
            % Enable the grid
            grid on;
        end
    end
end


%% Plot all channel peak values over time
%% Plot channel peak values separately for side (L/R) and rings separately from segments

% Set up a colormap for different channel colors
colors = lines(total_channels);

% Define channel indices for each group
channels_1_to_12 = 1:12;
channels_13_to_30 = 13:30;

% Create four separate figures
for side = 1:2 % 1 for left, 2 for right
    for group = 1:2 % 1 for channels 1 to 12, 2 for channels 13 to 30
        figure; % Create a new figure
        
        % Create an array to store the line plot handles for the legend
        line_handles = gobjects(1, 0);
        
        if group == 1
            channel_indices = channels_1_to_12;
            title_str = 'Channels 1 to 12 ';
        else
            channel_indices = channels_13_to_30;
            title_str = 'Channels 13 to 30 ';
        end
        
        if side == 1
            title_str = [title_str 'Left'];
        else
            title_str = [title_str 'Right'];
        end

        % Iterate through each channel
        for channel_idx = channel_indices
            % Check if the channel is in the desired side
            if any(channel_side_all(channel_idx, :) == side)
                % Get the peak values for the current channel over time
                peak_values = peak_value_all(channel_idx, :);

                % Remove NaN values and get the corresponding dates
                non_nan_indices = ~isnan(peak_values);
                peak_values_no_nan = peak_values(non_nan_indices);
                dates_no_nan = dates(non_nan_indices);

                % Plot the peak values for the current channel over time with a unique color
                line_handle = plot(dates_no_nan, peak_values_no_nan, 'Color', colors(channel_idx, :), 'DisplayName', channel_labels{channel_idx});
                hold on;
                line_handles = [line_handles, line_handle];

                % Get significant and non-significant peak indices
                significant_indices = significant_peak_all(channel_idx, :) == 1;
                non_significant_indices = significant_peak_all(channel_idx, :) == 0;

                % Plot significant peaks with a green dot and non-significant peaks with a red dot
                scatter(dates(significant_indices & non_nan_indices), peak_values(significant_indices & non_nan_indices), 'g', 'filled');
                scatter(dates(non_significant_indices & non_nan_indices), peak_values(non_significant_indices & non_nan_indices), 'r', 'filled');
            end
        end

        xlabel('Time');
        ylabel('Peak Amplitude');
        title(title_str);
        legend(line_handles, 'Location', 'bestoutside');
        set(gca, 'YScale', 'log');
        grid on;
    end
end



%% Save matrices, file info, and figures

% Get the last folder in the 'folder_path' variable
[parent_folder, last_folder] = fileparts(folder_path);

% Create an output directory in the current working folder with a "_PSD" suffix
output_dir = fullfile(pwd, [last_folder, '_PSD']);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end


% Save all the matrices as '.mat' files
save(fullfile(output_dir, 'peak_value_all.mat'), 'peak_value_all');
save(fullfile(output_dir, 'peak_frequency_all.mat'), 'peak_frequency_all');
save(fullfile(output_dir, 'peak_prominence_all.mat'), 'peak_prominence_all');
save(fullfile(output_dir, 'significant_peak_all.mat'), 'significant_peak_all');
save(fullfile(output_dir, 'channel_labels.mat'), 'channel_labels');
save(fullfile(output_dir, 'channel_side_all.mat'), 'channel_side_all');
save(fullfile(output_dir, 'dates.mat'), 'dates');
save(fullfile(output_dir, 'psd_results_all.mat'), 'psd_results_all');
save(fullfile(output_dir, 'significant_peak_all.mat'), 'significant_peak_all');

% Save all the figures as '.png' files
figHandles = findall(groot, 'Type', 'figure');
for fig_idx = 1:numel(figHandles)
    fig_handle = figHandles(fig_idx);
    fig_name = sprintf('figure_%d.png', fig_idx);
    saveas(fig_handle, fullfile(output_dir, fig_name));
end

%% SAVE FILE INFO TO TXT FILE

% Open the text file for writing
fileID = fopen(fullfile(output_dir, 'file_summary.txt'), 'w');

% Iterate through each file
for file_idx = 1:num_files
    % Write the file name
    fprintf(fileID, 'File: %s\n', mat_files_LMTD{file_idx});
    fprintf(fileID, 'Datetime: %s\n', datestr(dates(file_idx)));
    
    % Write the header for the channel information
    fprintf(fileID, 'Channel\t\tSide\tPeak Value\tPeak Frequency\tPeak Prominence\tSignificant\n');
    channel_indices = find(~isnan(channel_side_all(:, file_idx)));
    present_channels = channel_labels(channel_indices);
    fprintf(fileID, 'Channel Labels Present: %s\n', strjoin(present_channels, ', '));
    fprintf(fileID, 'Channel Numbers Present: %s\n', mat2str(channel_indices'));

    % Iterate through each channel
    for channel_idx = 1:total_channels
        % Check if there is data in the corresponding fields
        if ~isnan(channel_side_all(channel_idx, file_idx))
            % Get the side
            if channel_side_all(channel_idx, file_idx) == 1
                side = 'Left';
            else
                side = 'Right';
            end
            
            % Write the channel information
            fprintf(fileID, '%s\t%s\t%.6f\t%.6f\t%.6f\t%d\n', channel_labels{channel_idx}, side, ...
                    peak_value_all(channel_idx, file_idx), peak_frequency_all(channel_idx, file_idx), ...
                    peak_prominence_all(channel_idx, file_idx), significant_peak_all(channel_idx, file_idx));
        end
    end
    
    % Add a separator between files
    fprintf(fileID, '-------------------------------------------------------------\n');
end

% Close the text file
fclose(fileID);


