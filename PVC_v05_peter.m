% This programm reads ECG data which are saved in format 212.
% (e.g., 100.dat from MIT-BIH-DB, cu01.dat from CU-DB,...)
% The data are displayed in a figure together with the annotations.
% The annotations are saved in the vector ANNOT, the corresponding
% times (in seconds) are saved in the vector ATRTIME.
% The annotations are saved as numbers, the meaning of the numbers can
% be found in the codetable "ecgcodes.h" available at www.physionet.org.
%
% ANNOT only contains the most important information, which is displayed
% with the program rdann (available on www.physionet.org) in the 3rd row.
% The 4th to 6th row are not saved in ANNOT.
%
%
%      created on Feb. 27, 2003 by
%      Robert Tratnig (Vorarlberg University of Applied Sciences)
%      (email: rtratnig@gmx.at),
%
%      algorithm is based on a program written by
%      Klaus Rheinberger (University of Innsbruck)
%      (email: klaus.rheinberger@uibk.ac.at)
%
% -------------------------------------------------------------------------
clc; clear all;

%------ SPECIFY DATA ------------------------------------------------------
PATH= 'C:\Users\Daniel\Documents\FHKärnten\Medical Engineering\1. Semester\Applied Medical Signal Analysis\MSA Projekt\Test'; % path, where data are saved
% PATH= 'C:\Users\Daniel\Documents\FHKärnten\Medical Engineering\1. Semester\Applied Medical Signal Analysis\MSA Projekt\Test'; % path, where data are saved

FileID = 101
FileID = num2str(FileID);

HEADERFILE= [FileID,'.hea'];      % header-file in text format
ATRFILE= [FileID,'.atr'];         % attributes-file in binary format
DATAFILE=[FileID,'.dat'];         % data-file
SAMPLES2READ=10000;         % number of samples to be read
                            % in case of more than one signal:
                            % 2*SAMPLES2READ samples are read

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
xlabel('Time / s'); ylabel('Voltage / mV');
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
end;
xlim([TIME(1), TIME(end)]);
xlabel('Time / s'); ylabel('Voltage / mV');
string=['ECG signal ',DATAFILE, ' - zero-phase filtered'];
title(string);
fprintf(1,"displaying highpass filtered image \n");


% ------FEATURE EXTRACTION-------------------------------------------------
% Calculating Moving Mean
window_size=30;
sig_mov_mean=movmean(sig_filter,window_size); 

% Finding Peaks (of QRS Wave) - 
[ecgpeaks, peaktimes] = findpeaks(sig_filter, TIME, 'MinPeakProminence', 0.7, 'MinPeakDistance', 0.15);  % Prominence = Threshold

HeartrateBPM_mean = numel(ecgpeaks)*60/peaktimes(end)
HeartrateBPS_mean = numel(ecgpeaks)/peaktimes(end)

%Creating a dataframe
dataframe = ones(length(TIME), 4) * 50;

dataframe(:, 1) = TIME';
dataframe(:, 2) = sig_filter(:,1);
dataframe(:, 3) = sig_mov_mean(:,1);

% Find indices of peaktimes in TIME
[~, peakIndices] = ismember(peaktimes, TIME);

% Ensure peakIndices are valid
peakIndices = peakIndices(peakIndices > 1 & peakIndices < length(sig_filter));

% Set a maximum number of values to the left and right (adjust as needed)
maxValues = 20;

% Iterate through each peak index
for i = 1:length(peakIndices)
    % Get the index of the current peak
    peakIndex = peakIndices(i);

    % Extract values to the left until reaching 0
    leftValues = sig_filter(peakIndex - 1 : -1 : max(peakIndex - maxValues, 1));
    leftValues = leftValues(leftValues > 0);

    % Extract values to the right until reaching 0
    rightValues = sig_filter(peakIndex + 1 : min(peakIndex + maxValues, end));
    rightValues = rightValues(rightValues > 0);

    dataframe(peakIndex-(size(leftValues)+sfreq/20):peakIndex+(size(rightValues)+sfreq/5), 4) = 1;
end

% figure(4); clf, box on, hold on
% plot(dataframe(:,1), dataframe(:,2),'r');
% 
% plot(TIME, sig_mov_mean,'b', 'DisplayName', 'Moving Mean');
% plot(peaktimes, ecgpeaks, 'mo', 'MarkerSize', 5, 'DisplayName', 'QRS Peaks');
% 
% % Find consecutive regions with the same values in dataframe(:,3)
% consecutiveRegions = bwconncomp(diff(dataframe(:,3)) == 0);
% regionLabels = labelmatrix(consecutiveRegions);
% 
% % Plot colored regions
% 
% uniqueRegions = unique(regionLabels);
% for i = 1:length(uniqueRegions)
%     regionIdx = regionLabels == uniqueRegions(i);
%     scatter(dataframe(regionIdx, 1), dataframe(regionIdx, 2), 10, dataframe(regionIdx, 3), 'filled');
% end
% 
% for k = 1:length(ATRTIMED)
%     text(ATRTIMED(k), 0, num2str(ANNOTD(k)));
% end
% 
% xlim([dataframe(1,1), dataframe(end,1)]);
% xlabel('Time / s'); ylabel('Voltage / mV');
% string = ['ECG signal ', DATAFILE, ' - labeled regions'];
% title(string);
% fprintf(1,'Displaying label regions \n');


% Display or use the selected values as needed

figure(3); clf, box on, hold on

plot(TIME, sig_filter,'r', 'DisplayName', 'Filtered ECG signal');

plot(TIME, sig_mov_mean,'b', 'DisplayName', 'Moving Mean');

plot(peaktimes, ecgpeaks, 'mo', 'MarkerSize', 5, 'DisplayName', 'QRS Peaks');

mask_QRS=dataframe(:,4)==1;

plot(TIME,mask_QRS,'m','LineWidth',1,'DisplayName', 'QRS-Label' );

for k=1:length(ATRTIMED)
    text(ATRTIMED(k),0,num2str(ANNOTD(k)));
end;
dummy_plot = plot(NaN, NaN, 'k.', 'DisplayName', 'Attribute labels');

legend('show');

xlabel('Time / s'); ylabel('Voltage / mV');
string=['Filtered ECG signal ',DATAFILE, ' Feature Extraction'];
title(string);
fprintf(1,"displaying feature extraction filtered image \n");



% ------ FEATURE EXTRACTION AND WINDOW LABELING -----------------------------
% Square the data in dataframe(:,2)
dataframe(:, 2) = dataframe(:, 2) .^ 2;

% Normalize squared voltage data column
dataframe(:, 2) = normalize(dataframe(:, 2), 'norm');

% Set up windows and assign new labels
windowSize = 200; % Adjust as needed

% Initialize variables to store windowed data
windowedData = [];
windowedTime = [];
windowedAnnotations = [];
windowedLabels = [];

% Iterate through each peak index
for i = 1:length(peakIndices)
    % Get the index of the current peak
    peakIndex = peakIndices(i);

    % Define the window indices
    windowStart = max(1, peakIndex - windowSize / 2);
    windowEnd = min(length(dataframe), peakIndex + windowSize / 2);

    % Extract windowed data, time, and annotations
    windowedData = [windowedData; dataframe(windowStart:windowEnd, 2)];
    windowedTime = [windowedTime; dataframe(windowStart:windowEnd, 1)];
    windowedAnnotations = [windowedAnnotations; dataframe(windowStart:windowEnd, 3)];

    % Find the most common annotation in the window
    commonAnnotation = mode(dataframe(windowStart:windowEnd, 3));

    % Assign the common annotation as the label for the window
    windowedLabels = [windowedLabels; commonAnnotation];
end

% ------ NEURAL NETWORK TRAINING -------------------------------------------
% Normalize windowed data
windowedData = normalize(windowedData, 'norm');

% Create the target matrix
numOutputClasses = max(windowedLabels);
targetMatrix = full(ind2vec(windowedLabels'));

% Define the neural network architecture
net = patternnet([20, 25, 50, 55, numOutputClasses]);

% Split the data into training and testing sets
trainingPercentage = 80;
nWindows = size(windowedData, 1);
seed = 40;
rng(seed);

cp = cvpartition(nWindows, 'HoldOut', (100 - trainingPercentage) / 100);

% Use logical indexing to select training data
trainingData = windowedData(training(cp), :);
trainingTargets = targetMatrix(:, training(cp));

% Use logical indexing to select testing data
testingData = windowedData(test(cp), :);
testingTargets = targetMatrix(:, test(cp));

% Train the neural network
net.trainParam.epochs = 25;
net.trainParam.showWindow = true; % Display training progress
net.trainParam.max_fail = 6; % Maximum validation failures
net.trainParam.lr = 0.01; % learning rate
net.trainParam.goal = 1e-5; % weight decay

net = train(net, trainingData', trainingTargets');

% Test the neural network
predictions = net(testingData');

% Convert predictions to indices
[~, predictedLabels] = max(predictions, [], 1);

% Convert testingTargets to indices
[~, trueLabels] = max(testingTargets, [], 1);

% Evaluate performance
accuracy = sum(predictedLabels == trueLabels) / numel(trueLabels);
disp(['Accuracy: ', num2str(accuracy)]);




% % ------CLASSICIFICATION W/ SVM--------------------------------------------
% 
% rand=randperm(size(sig_filter,1)); %returns a row vector containing a random permutation of the integers from 1 to n without repeating elements.
% num_svm=round(size(sig_filter,1))*0.75;
% 
% % Masking the Labels
% mask_QRS=dataframe(:,4)==1;
% 
% points =size(sig_mov_mean,1);
% t=1:1/sfreq:points*1/sfreq+1-1/sfreq;
% 
% % Training Data
% 
% xtr = [sig_mov_mean(rand(1:num_svm))];
% ytr = mask_QRS(rand(1:num_svm));
% 
% % Test Data
% xt = [sig_mov_mean(rand(num_svm:end))];
% yt = mask_QRS(rand(num_svm:end));
% 
% % SVM
% model = fitcsvm(xtr, ytr,'KernelFunction','RBF');
% 
% % Predict
% result = predict(model, xt);
% accuracy = sum(result==yt)/length(yt)*100;
% sp = sprintf("Test accuracy = %.2f", accuracy);
% disp(sp);
% 
% % Plot
% result = predict(model, [sig_mov_mean]);
% 
% 
% figure; set(gcf,'color','w');
% plot(TIME(1:10000),sig_filter(1:10000).*250,'LineWidth',1.5); xlabel('time in s'); ylabel('mV'); 
% hold on; plot(TIME(1:10000),mask_QRS(1:10000).*50,'r','LineWidth',2);
% hold on; plot(TIME(1:10000),result(1:10000).*60,'k:','LineWidth',2 );title(sprintf("Test accuracy = %.2f", accuracy));
% legend('filtered signal','QRS Label','QRS predicted')
% 
% 
% %% SET BREAKPOINT HERE
% breakpoint =1
% 
% %  ------------------------------
% 
% % xtr = [xA xB; yA yB];
% % ytr = mask';
% % 
% % figure; set(gcf,'color','w'); 
% % 
% % plot(xtr(find(mask==1),1),xtr(find(mask==1),2),'+g','MarkerSize',10, 'LineWidth',3);
% % hold on; plot(xtr(find(mask==2),1),xtr(find(mask==2),2),'+b','MarkerSize',10, 'LineWidth',3);
% % xlabel('x'); ylabel('y')
% % % xlim([0 7]); ylim([0 5])
% % title('Training Data')
% % legend('Class 1', 'Class 2')
% % 
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%% SVM %%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % model = fitcsvm(xtr, ytr,'KernelFunction','linear', 'KernelScale','auto');
% % %model = fitcsvm(xtr, ytr,'KernelFunction','polynomial', 'KernelScale','auto');
% % %model = fitcsvm(xtr, ytr,'KernelFunction','rbf');
% % 
% % 
% % sv = model.SupportVectors; % Support vectors
% % %beta = model.Beta % Linear predictor coefficients
% % %b = model.Bias % Bias term
% % 
% % figure; set(gcf,'color','w'); 
% % plot(xtr(find(mask==1),1),xtr(find(mask==1),2),'+g','MarkerSize',10, 'LineWidth',3);
% % hold on; plot(xtr(find(mask==2),1),xtr(find(mask==2),2),'+b','MarkerSize',10, 'LineWidth',3);
% % hold on;  plot(sv(:,1),sv(:,2),'ro','MarkerSize',10, 'LineWidth',3)
% % xlabel('x'); ylabel('y')
% % % xlim([0 7]); ylim([0 5])
% % title('Input Data with Support Vectors')
% % legend('Class 1', 'Class 2', 'Support Vectors')
% % 
% % 
% % %%%%%%  show Classification Pattern
% % 
% % h = 0.2;
% % [X1,X2] = meshgrid(0:h:200,0:h:200);
% % u = [X1(:),X2(:)];
% % result = predict(model,u);
% % scoreGrid = reshape(result,size(X1,1),size(X2,2));
% % 
% % figure; set(gcf,'color','w'); 
% % hold on; plot(u(find(result==1),1),u(find(result==1),2),'+r','MarkerSize',1, 'LineWidth',2);
% % hold on; plot(u(find(result==2),1),u(find(result==2),2),'+k','MarkerSize',1, 'LineWidth',2);
% % hold on; plot(xtr(find(mask==1),1),xtr(find(mask==1),2),'+g','MarkerSize',9, 'LineWidth',2);
% % hold on; plot(xtr(find(mask==2),1),xtr(find(mask==2),2),'+b','MarkerSize',9, 'LineWidth',2);
% % hold on; contour(X1,X2,scoreGrid)
% % 
% % title('Classification Pattern')
% % 
% % 
% 
% 
% % -----------------------------
% 
% 
% % for i = 1:size(idx,1)
% %     temp = idx(i,1);
% %     dataframe(temp,3) = ANNOTD(i,1);
% % end
% 
% windowSize = 70;
% 
% for i = 1:size(idx, 1)
%     centerIndex = idx(i, 1);
%     % Define the window indices
%     windowStart = max(1, centerIndex - windowSize/2);
%     windowEnd = min(size(dataframe, 1), centerIndex + windowSize/2);
% 
%     % Assign ANNOTD to the window
%     dataframe(windowStart:windowEnd, 3) = ANNOTD(i, 1);
% end
% 
% 
% % Find the corresponding indices in the new array for each ATRTIME
% [~, idx] = ismember(ATRTIMED, TIME);
% 
% windowSize = 70;
% 
% for i = 1:size(idx, 1)
%     centerIndex = idx(i, 1);
%     % Define the window indices
%     windowStart = max(1, centerIndex - windowSize/2);
%     windowEnd = min(size(dataframe, 1), centerIndex + windowSize/2);
% 
%     % Assign ANNOTD to the window
%     dataframe(windowStart:windowEnd, 3) = ANNOTD(i, 1);
% end
% 
% 
% % for i=1:size(idx,1)
% %     temp = idx(i,1);
% %     for j=1:windowSize
% %         if idx(1,1) > -windowSize/2 | windowSize/2 < idx(end,1)
% %             dataframe(temp+j,3) = ANNOTD(temp,1);
% %         else
% % 
% %         end
% %     end
% % end
% 
% 
% figure(4); clf, box on, hold on
% plot(dataframe(:,1), dataframe(:,2),'r');
% 
% % Find consecutive regions with the same values in dataframe(:,3)
% consecutiveRegions = bwconncomp(diff(dataframe(:,3)) == 0);
% regionLabels = labelmatrix(consecutiveRegions);
% 
% % Plot colored regions
% uniqueRegions = unique(regionLabels);
% for i = 1:length(uniqueRegions)
%     regionIdx = regionLabels == uniqueRegions(i);
%     scatter(dataframe(regionIdx, 1), dataframe(regionIdx, 2), 10, dataframe(regionIdx, 3), 'filled');
% end
% 
% for k = 1:length(ATRTIMED)
%     text(ATRTIMED(k), 0, num2str(ANNOTD(k)));
% end
% xlim([dataframe(1,1), dataframe(end,1)]);
% xlabel('Time / s'); ylabel('Voltage / mV');
% string = ['ECG signal ', DATAFILE, ' - labeled regions'];
% title(string);
% fprintf(1,'Displaying label regions \n');
% 
% 
% % ------ CLASSIFICATION APPROACH ------------------------------------------------
% % Change sliding window size for classification
% windowSize = 30;
% 
% % Normalize voltage data column
% % dataframe(:, 2) = normalize(dataframe(:, 2), 'norm');
% 
% % Create the target matrix
% labels = dataframe(:, 3);
% 
% % Define the neural network architecture
% numOutputClasses = 50;
% net = patternnet([20, 25, 50, 55, numOutputClasses]);
% 
% % Split the data into training and testing sets
% trainingPercentage = 80;
% n = size(dataframe, 1);
% seed = 40;
% rng(seed);
% 
% cp = cvpartition(n, 'HoldOut', (100 - trainingPercentage)/100);
% 
% % Initialize variables for training and testing data
% trainingData = [];
% trainingLabels = [];
% testingData = [];
% testingLabels = [];
% 
% % Iterate over the data with a sliding window
% for i = 1:(n - windowSize + 1)
%     windowData = dataframe(i:(i + windowSize - 1), 2); 
%     windowLabels = labels(:, 1);
% 
%     % Assign to training or testing set based on the partition
%     if training(cp, i + windowSize - 1)
%         trainingData = [trainingData; windowData];
%         trainingLabels = [trainingLabels; windowLabels];
%     else
%         testingData = [testingData; windowData];
%         testingLabels = [testingLabels; windowLabels];
%     end
% end
% 
% % Train the neural network
% net.trainParam.epochs = 25;
% net.trainParam.showWindow = true; % Display training progress
% net.trainParam.max_fail = 6; % Maximum validation failures
% net.trainParam.lr = 0.01; % learning rate
% net.trainParam.goal = 1e-5; % weight decay
% 
% net = train(net, trainingData(:, :)', trainingLabels');
% 
% % Test the neural network
% predictions = net(testingData(:, :)'');
% 
% % Convert predictions to indices
% [~, predictedLabels] = max(predictions, [], 1);
% 
% % Convert testingLabels to indices
% [~, trueLabels] = max(testingLabels', [], 1);
% 
% % Evaluate performance
% accuracy = sum(predictedLabels == trueLabels) / numel(trueLabels);
% disp(['Accuracy: ', num2str(accuracy)]);