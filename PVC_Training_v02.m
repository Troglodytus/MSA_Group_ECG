clc; clear all; close all;

%------ SPECIFY DATA ------------------------------------------------------
PATH= 'C:\Users\Daniel\Documents\FHKärnten\Medical Engineering\1. Semester\Applied Medical Signal Analysis\MSA Projekt\Test';

% List of file IDs to process
fileIDs = [100, 101, 105, 103, 104, 106, 112, 114, 115, 116, 117, 200, 208,214,221];

% Initialize variables to store aggregated data
allNormalizedWindows = [];
allTargetMatrix = [];

for fileID = fileIDs
    FileID = num2str(fileID);

    HEADERFILE= [FileID,'.hea'];
    ATRFILE= [FileID,'.atr'];
    DATAFILE=[FileID,'.dat'];
    SAMPLES2READ=200000;

    %------ LOAD HEADER DATA --------------------------------------------------
    fprintf(1,'loading %s ...\n', HEADERFILE);
    signalh= fullfile(PATH, HEADERFILE);
    fid1=fopen(signalh,'r');
    z= fgetl(fid1);
    A= sscanf(z, '%*s %d %d %d',[1,3]);
    nosig= A(1);  % number of signals
    sfreq=A(2);   % sample rate of data
    clear A;
    for k=1:nosig
        z= fgetl(fid1);
        A= sscanf(z, '%*s %d %d %d %d %d',[1,5]);
        dformat(k)= A(1);           % format; here only 212 is allowed
        gain(k)= A(2);              % number of integers per mV
        bitres(k)= A(3);            % bitresolution
        zerovalue(k)= A(4);         % integer value of ECG zero point
        firstvalue(k)= A(5);        % first integer value of signal (to test for errors)
    end;
    fclose(fid1);
    clear A;

    %------ LOAD BINARY DATA --------------------------------------------------
    if dformat~= [212,212], error('this script does not apply binary formats different to 212.'); end;
    signald= fullfile(PATH, DATAFILE);            % data in format 212
    fid2=fopen(signald,'r');
    A= fread(fid2, [3, SAMPLES2READ], 'uint8')';  % matrix with 3 rows, each 8 bits long, = 2*12bit
    fclose(fid2);
    M2H= bitshift(A(:,2), -4);
    M1H= bitand(A(:,2), 15);
    PRL=bitshift(bitand(A(:,2),8),9);     % sign-bit
    PRR=bitshift(bitand(A(:,2),128),5);   % sign-bit
    M( : , 1)= bitshift(M1H,8)+ A(:,1)-PRL;
    M( : , 2)= bitshift(M2H,8)+ A(:,3)-PRR;
    if M(1,:) ~= firstvalue, error('inconsistency in the first bit values'); end;
    switch nosig
    case 2
        M( : , 1)= (M( : , 1)- zerovalue(1))/gain(1);
        M( : , 2)= (M( : , 2)- zerovalue(2))/gain(2);
        TIME=(0:(SAMPLES2READ-1))/sfreq;
    case 1
        M( : , 1)= (M( : , 1)- zerovalue(1));
        M( : , 2)= (M( : , 2)- zerovalue(1));
        M=M';
        M(1)=[];
        sM=size(M);
        sM=sM(2)+1;
        M(sM)=0;
        M=M';
        M=M/gain(1);
        TIME=(0:2*(SAMPLES2READ)-1)/sfreq;
    otherwise  % this case did not appear up to now!
        % here M has to be sorted!!!
        disp('Sorting algorithm for more than 2 signals not programmed yet!');
    end;
    clear A M1H M2H PRR PRL;
    fprintf(1,'Loading data finished \n');

    %------ LOAD ATTRIBUTES DATA ----------------------------------------------
    atrd= fullfile(PATH, ATRFILE);      % attribute file with annotation data
    fid3=fopen(atrd,'r');
    A= fread(fid3, [2, inf], 'uint8')';
    fclose(fid3);
    ATRTIME=[];
    ANNOT=[];
    sa=size(A);
    saa=sa(1);
    i=1;
    while i<=saa
        annoth=bitshift(A(i,2),-2);
        if annoth==59
            ANNOT=[ANNOT;bitshift(A(i+3,2),-2)];
            ATRTIME=[ATRTIME;A(i+2,1)+bitshift(A(i+2,2),8)+...
                    bitshift(A(i+1,1),16)+bitshift(A(i+1,2),24)];
            i=i+3;
        elseif annoth==60
            % nothing to do!
        elseif annoth==61
            % nothing to do!
        elseif annoth==62
            % nothing to do!
        elseif annoth==63
            hilfe=bitshift(bitand(A(i,2),3),8)+A(i,1);
            hilfe=hilfe+mod(hilfe,2);
            i=i+hilfe/2;
        else
            ATRTIME=[ATRTIME;bitshift(bitand(A(i,2),3),8)+A(i,1)];
            ANNOT=[ANNOT;bitshift(A(i,2),-2)];
       end;
       i=i+1;
    end;
    ANNOT(length(ANNOT))=[];       % last line = EOF (=0)
    ATRTIME(length(ATRTIME))=[];   % last line = EOF
    clear A;
    ATRTIME= (cumsum(ATRTIME))/sfreq; % 
    ind= find(ATRTIME <= TIME(end));
    ATRTIMED= ATRTIME(ind);
    ANNOT=round(ANNOT);
    ANNOTD= ANNOT(ind);
    
    %------ DISPLAY DATA ------------------------------------------------------
    figure(1); clf, box on, hold on
    plot(TIME, M(:,1),'r');
    if nosig==2
        plot(TIME, M(:,2),'b');
    end;
    for k=1:length(ATRTIMED)
        text(ATRTIMED(k),0,num2str(ANNOTD(k)));
    end;
    xlim([TIME(1), TIME(end)]);
    xlabel('Time [s]'); ylabel('Voltage [mV]');
    string=['ECG signal ',DATAFILE];
    title(string);
    fprintf(1,'displaying original data \n');
    
    % ------HIGHPASS FILTERING-------------------------------------------------
    
    load test.mat
    Filter_FIR=test;
    %load ECG_FIR_Highpass_02.mat
    %Filter_FIR=ECG_FIR_Highpass_02;
    
    Numerator = Filter_FIR.Numerator;
    sig_filter = filtfilt(Numerator,1,M(:,1));
    
    % Plotting Figure 2 (Original Data, Filtered Data, Atr. Labels) 
    figure(2); clf, box on, hold on
    plot(TIME, sig_filter,'b');
    plot(TIME, M(:,1),'r');
    for k=1:length(ATRTIMED)
        text(ATRTIMED(k),0,num2str(ANNOTD(k)));
    end
    xlim([TIME(1), TIME(end)]);
    xlabel('Time [s]'); ylabel('Voltage [mV]');
    string=['ECG signal ',DATAFILE, ' - high-pass filtered'];
    title(string);
    fprintf(1,"displaying highpass filtered image \n");
    
    % ------ FEATURE EXTRACTION AND WINDOW LABELING -----------------------------
    % Calculating Moving Mean
    window_size=30;
    sig_mov_mean=movmean(sig_filter,window_size); 
    
    % Finding Peaks (of QRS Wave) - 
    [ecgpeaks, peaktimes] = findpeaks(sig_filter, TIME, 'MinPeakProminence', 0.8, 'MinPeakDistance', 0.2);  % Prominence = Threshold
    
    ecgpeaks_norm = ecgpeaks/max(ecgpeaks);
    
    % Creating a dataframe
    dataframe = ones(length(TIME), 4) * 0;
    
    dataframe(:, 1) = TIME';
    dataframe(:, 2) = sig_filter(:,1).^2;  % Square the data
    dataframe(:, 3) = sig_mov_mean(:,1);
    
    % Normalize the squared data between 0 and 1
    dataframe(:, 2) = (dataframe(:, 2) - min(dataframe(:, 2))) / (max(dataframe(:, 2)) - min(dataframe(:, 2)));
    dataframe(:, 2) = sqrt(dataframe(:, 2));
    
    
    figure(3); clf, box on, hold on
    plot(dataframe(:, 1), dataframe(:, 2), 'b');
    % plot(dataframe(:, 1), dataframe(:, 3), 'g');
    for k = 1:length(ATRTIMED)
        text(ATRTIMED(k), 0, num2str(ANNOTD(k)));
    end
    scatter(peaktimes, ecgpeaks_norm, 'r', 'filled');
    xlim([dataframe(1, 1), dataframe(end, 1)]);
    xlabel('Time [s]'); ylabel('Voltage [mV]');
    string=['Filtered ECG signal - ', DATAFILE];
    title(string);
    fprintf(1, 'displaying filtered data\n');
    
    
    
    % Find indices of peaktimes in TIME
    [~, peakIndices] = ismember(peaktimes, TIME);
    
    % Ensure peakIndices are valid
    peakIndices = peakIndices(peakIndices > 1 & peakIndices < length(sig_filter));
    
    % Find the nearest index in TIME for each ATRTIME
    % [~, idx] = min(abs(bsxfun(@minus, ATRTIME(:), TIME(:)')));
    [~, idx] = ismember(ATRTIMED, TIME);
    
    for i = 1:size(idx, 1)
        % Assign ANNOT to the dataframe
        dataframe(idx(i), 4) = ANNOTD(i, 1);
    end
    
    % Find the indices of label 1 in the dataframe
    label_1_indices = find(dataframe(:, 4) == 1);
    
    % Create individual windows
    windows = cell(length(label_1_indices)-1, 1);
    window_labels = zeros(length(label_1_indices)-1, 1);
    
    for i = 1:length(label_1_indices)-1
        % Extract windowed data and time
        window_start = label_1_indices(i);
        window_end = label_1_indices(i+1) - 1;
        
        window_data = dataframe(window_start:window_end, 2);
        
        % Find the most common non-zero label in the window
        labels_in_window = dataframe(window_start:window_end, 4);
        non_zero_labels = labels_in_window(labels_in_window > 1);
        
        if ~isempty(non_zero_labels)
            common_label = mode(non_zero_labels);
        else
            common_label = 0;
        end
        
        % Store the window and label
        windows{i} = window_data';
        window_labels(i) = common_label;
    end
    
    % Create a new array for classes 'normal' and 'arrhythmia'
    class_labels = cell(length(window_labels), 1);
    
    for i = 1:length(window_labels)
        if ismember(window_labels(i), [0, 1, 12, 20, 21, 38])
            class_labels{i} = 'normal';
        else
            class_labels{i} = 'arrhythmia';
        end
    end
    
    % Convert class labels to numerical values if needed
    numeric_labels = ones(size(class_labels));
    numeric_labels(strcmp(class_labels, 'arrhythmia')) = 2;
    
    % Normalize windows in length
    normalized_windows = zeros(length(windows), 256);
    
    for i = 1:length(windows)
            % Interpolate to the 256 length
      normalized_windows(i, :) = interp1(linspace(0, 1, length(windows{i})), windows{i}, linspace(0, 1, 256));
    end
    
        % Append current file's data to aggregated data
    allNormalizedWindows = [allNormalizedWindows; normalized_windows];
    allTargetMatrix = [allTargetMatrix, numeric_labels'];
    
end  
        
        % normalized_windows = normalize(normalized_windows);


%plot windows and labels
figure(4); clf, box on, hold on

temp = 0;
total = 0;
for i = 1:length(windows)
    window_start = length(windows{i}) - (length(windows{i}) - 1) + total;
    window_end = length(windows{i}) - 1 + total;
    
    % Check the class label and plot colored line accordingly
    if strcmp(class_labels{i}, 'normal')  % Use strcmp for string comparison
        plot(dataframe(window_start:window_end, 1), dataframe(window_start:window_end, 2), 'Color', [0, 0.5, 0], 'LineWidth', 1);
    elseif strcmp(class_labels{i}, 'arrhythmia')  % Use strcmp for string comparison
        plot(dataframe(window_start:window_end, 1), dataframe(window_start:window_end, 2), 'r-', 'LineWidth', 1);
    end
    
    % text(dataframe(window_start + round((window_end - window_start) / 2), 1), 0.5, class_labels{i}, 'HorizontalAlignment', 'center');
    
    temp = length(windows{i});
    total = total + temp;
end

% Additional customization 
xlabel('Time [s]'); ylabel('Voltage norm.');
legend('normal', 'arrhythmia');  
string = ['ECG Actual Labels of ', num2str(fileIDs)];
title(string);
fprintf(1, 'displaying classified windows\n');
hold off

% Plotting the data
figure(5);
hold on;

for i = 1:size(allNormalizedWindows, 1)
    if allTargetMatrix(i) == 1
        plot(allNormalizedWindows(i, :), 'Color', [0, 0.5, 0]);
    else
        plot(allNormalizedWindows(i, :), 'r');

    end
end

hold off;

title('RR Intervals - Actual Labels');
xlabel('RR Interval norm. - Index');
ylabel('Voltage norm.');
legend('normal', 'arrhythmia');  


% ------ NEURAL NETWORK TRAINING -------------------------------------------
% Create the target matrix for normalized windows
numOutputClasses = 2;  % Two output classes: 'normal' and 'arrhythmia'
targetMatrix = numeric_labels';

% Define the neural network architecture
net = patternnet([10, 15, 35, 20, 10, 5, numOutputClasses],  'trainlm');
% Modify output layer activation function
net.layers{end}.transferFcn = 'softmax';

% Split the data into training and testing sets
trainingPercentage = 80;
nWindows = size(allNormalizedWindows, 1);
seed = 40;
rng(seed);


cp = cvpartition(nWindows, 'HoldOut', (100 - trainingPercentage) / 100);

% Use logical indexing to select training data
trainingData = allNormalizedWindows(cp.training, :);
trainingTargets = allTargetMatrix(:, cp.training)';

% Use logical indexing to select testing data
testingData = allNormalizedWindows(cp.test, :);
testingTargets = allTargetMatrix(:, cp.test)';

% Train the neural network
net.trainParam.epochs = 20;  % Change this to the desired number of epochs
net.trainParam.showWindow = true;  % Display training progress
net.trainParam.max_fail = 20;  % Maximum validation failures
net.trainParam.lr = 0.06;  % Learning rate
net.trainParam.goal = 1e-6;  % Weight decay

net = train(net, trainingData', trainingTargets');

% Save the trained neural network model to a file
save('ECG_model_DEMO_v01.mat', 'net');

% Test the neural network
predictions = round(net(testingData'),0);

% Convert predictions to indices
[~, predictedLabels] = max(predictions, [], 1);

% Transpose testingTargets if needed
[~, trueLabels] = max(testingTargets', [], 1);

% Create a confusion matrix
confMat = confusionmat(testingTargets', predictedLabels);

% Display the confusion matrix
disp('Confusion Matrix:');
disp(confMat);
% Create and display the confusion chart
figure;
confusionchart(testingTargets', predictedLabels, 'Normalization', 'row-normalized');
title('Confusion Chart - ECG Classification');

% Evaluate performance
accuracy = sum(predictions == testingTargets') / numel(testingTargets');
disp(['Accuracy: ', num2str(accuracy)]);