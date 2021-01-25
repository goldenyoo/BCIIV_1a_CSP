%% Call raw data
close all
clear all

% Ask user for input parameters
prompt = {'Data label: ', 'Feature vector length: ', 'Re-referencing: 0 (Non),1 (CAR), 2 (LAP)', 'BFP order'};
dlgtitle = 'Input';
dims = [1 50];
definput = {'a', '3', '1','20'};
answer = inputdlg(prompt,dlgtitle,dims,definput);


% Error detection
if isempty(answer), error("Not enough input parameters."); end

% Input parameters
data_label = string(answer(1,1));   % Calib_ds1 + "data_label"
m = double(string(answer(2,1))); % feature vector will have length (2m)
referencing = double(string(answer(3,1)));
order = double(string(answer(4,1)));

% Load file
FILENAME = strcat('C:\Users\유승재\Desktop\Motor Imagery EEG data\BCICIV_1_mat\BCICIV_calib_ds1',data_label,'.mat');
load(FILENAME);

% Data rescale
cnt= 0.1*double(cnt);
cnt = cnt';

% Exclude electrode (AF3, AF4, O1, O2, PO1, PO2)
cnt_c = cnt(3:55,:);
%% Preprocessing
if referencing ~= 0
    %%% Calculate differential voltage
    for i = 1 : size(cnt_c,1)
        cnt_c(i,:) = cnt_c(i,:) - cnt(29,:);
    end

    
    if referencing == 1 % common average
        Means = (1/size(cnt_c,1))*sum(cnt_c);
        for i = 1 : size(cnt_c,1)
            cnt_c(i,:) = cnt_c(i,:) - Means;
        end
    elseif referencing == 2 % LAP
        %%%
    end
end
%% 
%BPF Design
bpFilt = designfilt('bandpassfir','FilterOrder',order, ...
         'CutoffFrequency1',8,'CutoffFrequency2',30, ...
         'SampleRate',100);

% Apply BPF
for i = 1:size(cnt_c,1)
    cnt_c(i,:) = filtfilt(bpFilt, cnt_c(i,:));
end