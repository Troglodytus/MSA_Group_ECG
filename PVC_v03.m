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

FileID = 101
FileID = num2str(FileID);

HEADERFILE= [FileID,'.hea'];      % header-file in text format
ATRFILE= [FileID,'.atr'];         % attributes-file in binary format
DATAFILE=[FileID,'.dat'];         % data-file
SAMPLES2READ=30000;         % number of samples to be read
                            % in case of more than one signal:
                            % 2*SAMPLES2READ samples are read

%------ LOAD HEADER DATA --------------------------------------------------
fprintf(1,'\\n$> WORKING ON %s ...\n', HEADERFILE);
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
fprintf(1,'\\n$> LOADING DATA FINISHED \n');

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
ATRTIME= (cumsum(ATRTIME))/sfreq;
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
fprintf(1,'\\n$> DISPLAYING DATA FINISHED \n');

% ------HIGHPASS FILTERING-------------------------------------------------------------------

load ECG_FIR_Highpass_01.mat
Filter_FIR=ECG_FIR_Highpass_01;
Numerator = Filter_FIR.Numerator;

sig_filter = filter(Numerator,1,M(:,1));



figure(2); clf, box on, hold on
plot(TIME, sig_filter,'r');
for k=1:length(ATRTIMED)
    text(ATRTIMED(k),0,num2str(ANNOTD(k)));
end;
xlim([TIME(1), TIME(end)]);
xlabel('Time / s'); ylabel('Voltage / mV');
string=['ECG signal ',DATAFILE, ' - highpass filtered'];
title(string);
fprintf(1,'\\n$> DISPLAYING DATA FINISHED \n');


%Creating a dataframe
dataframe = ones(length(TIME), 3) * 50;


dataframe(:, 1) = TIME';
dataframe(:, 2) = sig_filter(:,1);
% Find the corresponding indices in the new array for each ATRTIME
[~, idx] = ismember(ATRTIMED, TIME);
windowSize = 10;

for i=1:size(idx,1)
    temp = idx(i,1);
    for j=1:windowSize
        if idx(1,1) > -windowSize/2 | windowSize/2 < idx(end,1)
            dataframe(temp+j,3) = ANNOTD(temp,1);
        else
            
        end
    end
end


% Fill the second column with the corresponding ANNOT values
dataframe(idx, 2) = ANNOTD(idx,1);





svm = [];

for j = 1:size(board)
    svm(j,1) = board(j,1);
    svm(j,2) = sqrt(board(j,2)*board(j,2)+board(j,3)*board(j,3)+board(j,4)*board(j,4));
end

aamv = [];

windowsize = 35;
for j = 1:size(svm)-windowsize
    for i=j:j+windowsize
        tempaamv(i,1) = abs((svm(j+1,2)-svm(j,2)));
    end
    aamv(j,1) = svm(j,1);
    aamv(j,2) = sum(tempaamv(:,1)) / windowsize;
    tempaamvboard = [];
end

figure;
sgtitle('SVM');
raise = false;
subplot(5, 1, 1);
plot(svmboard(:,1) / 1000, svmboard(:,2), 'r-', aamvboard(:,1) / 1000, aamvboard(:,2), 'b-');  % Plot
title('Board SVM', 'FontSize', 12);
legend('SVM', 'AAMV');


% Classification Approach

newlabel = [];
for i = 1:10
    newlabel{i,1} = 'board';
end


% Normalize X, Y, Z columns
data(:, 2:4, :) = normalize(data(:, 2:4, :), 'norm');

% Flatten the data for neural network input
flattenedData = permute(data, [1, 3, 2]); % Permute to have time as the second dimension
flattenedData = reshape(flattenedData, [], size(data, 1));

% Convert cell array of labels to numerical indices
labelStrings = string(newlabel);
[~, ~, labelIndices] = unique(labelStrings);

% Create the target matrix
labels = dummyvar(labelIndices);

% Define the neural network architecture
net = patternnet([20, 10, 5]); % 3 layers of hidden layers with 20, 10, 5 nodes


% Split the data into training and testing sets
trainingPercentage = 80;
n = size(data, 3);
seed = 40;
rng(seed);

cp = cvpartition(n, 'HoldOut', (100 - trainingPercentage)/100);

% Use logical indexing to select training data
trainingData = flattenedData(training(cp), :);
trainingLabels = labels(training(cp), :);

% Use logical indexing to select testing data
testingData = flattenedData(test(cp), :);
testingLabels = labels(test(cp), :);

% Train the neural network
net.trainParam.epochs = 40;
net.trainParam.showWindow = true; % Display training progress
net.trainParam.max_fail = 6; % Maximum validation failures
net.trainParam.lr = 0.01; %  learning rate
% net.trainParam.mc = 0.9;  %  momentum
net.trainParam.goal = 1e-5;  % weight decay

net = train(net, trainingData', trainingLabels');

% Test the neural network
predictions = net(testingData');

% Convert predictions to indices
[~, predictedLabels] = max(predictions, [], 1);

% Convert testingLabels to indices
[~, trueLabels] = max(testingLabels', [], 1);

% Evaluate performance
accuracy = sum(predictedLabels == trueLabels) / numel(trueLabels);
disp(['Accuracy: ', num2str(accuracy)]);