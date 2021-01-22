% ----------------------------------------------------------------------- %
%    File_name: Evaluation.m
%    Programmer: Seungjae Yoo                             
%                                           
%    Last Modified: 2020_01_22                            
%                                                            
 % ----------------------------------------------------------------------- %
 %% 
% close all
% clear all

load('C:\Users\유승재\Desktop\true_labels\feature.mat');
data_label = string(answer(1,1));
m = double(string(answer(2,1)));
referencing = double(string(answer(3,1)));
order = double(string(answer(4,1)));

FILENAME = strcat('C:\Users\유승재\Desktop\true_labels\BCICIV_eval_ds1',data_label,'_1000Hz_true_y.mat');
load(FILENAME);
%% 

A=[];
B=[];
C=[];
D=[];
%NaN, -1, 0, 1
tmp = 1;
for i=1:length(true_y)-1
    if isnan(true_y(i,1))
        if isnan(true_y(i+1,1))
            continue
        else
            A = [A; [tmp i]];
            tmp = i+1;
        end
    elseif true_y(i,1) ~= true_y(i+1,1)
        if true_y(i,1)== -1
            B = [B; [tmp i]];
        elseif true_y(i,1)== 0
            C = [C; [tmp i]];
        elseif true_y(i,1)== 1
            D = [D; [tmp i]];
        end
        tmp = i+1;
    end
end

%% 
% Load file
FILENAME = strcat('C:\Users\유승재\Desktop\Motor Imagery EEG data\BCICIV_1_mat\BCICIV_eval_ds1',data_label,'.mat');
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


%% 
f1 = figure;
% f2 = figure;
score = [];
predictions = [];
checks = [];
for j = 1 : length(B)
    tmp1 = round(B(j,1)/10);
    tmp2 = round(B(j,2)/10);
    E = cnt_c(:, tmp1:tmp2);
    Z = P'*E;
    
     % Feature vector
    tmp_ind = size(Z,1);
    Z_reduce = [Z(1:m,:); Z(tmp_ind-(m-1):tmp_ind,:)];
   
    var_vector = var(Z_reduce,0,2)';
    var_vector = (1/sum(var_vector))*var_vector;
        
    fp = log(var_vector);
    fp = fp';
    
    % Graphical represent
     figure(f1)
     scatter3(Z(1,:), Z(size(cnt_c,1),:),Z(2,:),'b'); hold on;
%      figure(f2)
%      scatter3(fp(1),fp(2),fp(6),'b'); hold on;
       
    % Run classifier
    [check, prediction] = myClassifier(fp);
    if prediction == -1
        score = [score 1];
        
    else
        score = [score 0];
        
    end
    predictions = [predictions prediction];
    checks = [checks check];
    
end

for j = 1 : length(D)
    tmp1 = round(D(j,1)/10);
    tmp2 = round(D(j,2)/10);
    E = cnt_c(:, tmp1:tmp2);
    Z = P'*E;
    
    % Feature vector
    tmp_ind = size(Z,1);
    Z_reduce = [Z(1:m,:); Z(tmp_ind-(m-1):tmp_ind,:)];
   
    var_vector = var(Z_reduce,0,2)';
    var_vector = (1/sum(var_vector))*var_vector;
        
    fp = log(var_vector);
    fp = fp';
    
    % Graphical represent
     figure(f1)
     scatter3(Z(1,:), Z(size(cnt_c,1),:),Z(2,:),'r'); hold on;
%      figure(f2)
%      scatter3(fp(1),fp(2),fp(6),'r'); hold on;
       
    % Run classifier
    [check, prediction] = myClassifier(fp);
    if prediction == 1
        score = [score 1];
        
    else
        score = [score 0];
        
    end
    predictions = [predictions prediction];
    checks = [checks check];
    
end


fprintf('Score: %g\n',100*sum(score)/length(score));
