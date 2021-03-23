% ----------------------------------------------------------------------- %
%    File_name: Eval.m
%    Programmer: Seungjae Yoo                             
%                                           
%    Last Modified: 2020_01_27                           
%                                                            
 % ----------------------------------------------------------------------- %
function predictions = Eval(answer,Mr,Ml,Qr,Ql,P,ref)
data_label = string(answer(1,1));   
m = double(string(answer(2,1))); % feature vector will have length (2m)
low_f = double(string(answer(3,1))); % Low cutoff freq
high_f = double(string(answer(4,1))); % High cutoff freq
sampling_rate = double(string(answer(5,1)));
referencing = double(string(answer(6,1))); % Non(0), CAR(1), LAP(2)
order = double(string(answer(7,1))); % Filter order




%% 
% Load file
if sampling_rate == 0
    FILENAME = strcat('C:\Users\유승재\Desktop\Motor Imagery EEG data\BCICIV_1_mat\BCICIV_eval_ds1',data_label,'.mat');
    chunk = 350;
    fs = 100;
else
    FILENAME = strcat('C:\Users\유승재\Desktop\Motor Imagery EEG data\BCICIV_1eval_1000Hz_mat\BCICIV_eval_ds1',data_label,'_1000Hz.mat');
    chunk = 3500;
    fs = 1000;
end
load(FILENAME);

% Data rescale
cnt= 0.1*double(cnt);
cnt = cnt';
%% Preprocessing
if referencing ~= 0
    %%% Calculate differential voltage
    for i = 1 : size(cnt,1)
        cnt(i,:) = cnt(i,:) - cnt(ref,:);
    end
    
    % common average
    if referencing == 1        
        cnt_c = cnt; % Exclude electrode (AF3, AF4, O1, O2, PO1, PO2)     
        Means = (1/size(cnt,1))*sum(cnt);
        for i = 1 : size(cnt_c,1)
            cnt_c(i,:) = cnt_c(i,:) - Means; % CAR
        end
        cnt_c = cnt_c([27 29 31 44 46 50 52 54],:);
    % LAP
    elseif referencing == 2
        cnt_n = myLAP(cnt,nfo); % Laplacian
        cnt_c = cnt_n([27 29 31 44 46 50 52 54],:); % Exclude electrode (AF3, AF4, O1, O2, PO1, PO2)
    end
else
        %%% Calculate differential voltage
    for i = 1 : size(cnt,1)
        cnt(i,:) = cnt(i,:) - cnt(ref,:);
    end
    
    cnt_c = cnt([27 29 31 44 46 50 52 54],:); % Exclude electrode (AF3, AF4, O1, O2, PO1, PO2)
end

clear cnt cnt_n

%% 
%BPF Design
bpFilt = designfilt('bandpassiir','SampleRate',fs,'PassbandFrequency1',low_f, ...
        'PassbandFrequency2',high_f,'StopbandFrequency1',low_f-2,'StopbandFrequency2',high_f+2, ...
        'StopbandAttenuation1',40,'StopbandAttenuation2',40, 'PassbandRipple',1,'DesignMethod','cheby2');
    

% Apply BPF
for i = 1:size(cnt_c,1)
    cnt_c(i,:) = filtfilt(bpFilt, cnt_c(i,:));
end
%%
FILENAME = strcat('C:\Users\유승재\Desktop\Motor Imagery EEG data\true_labels\BCICIV_eval_ds1',data_label,'_1000Hz_true_y.mat');
load(FILENAME);

true_y = downsample(true_y,10);
%% 
% f1 = figure;
% f2 = figure;
score = [];
predictions = [];
checks = [];

% For class 1
iter = 1;
chunk = 150;
predictions = zeros(1,size(cnt_c,2));

while iter + chunk <= size(cnt_c,2)
        
    E = cnt_c(:, iter:iter+chunk-1);
    
    for k = 1:3
        Z = P'*E;
        % Feature vector
        tmp_ind = size(Z,1);
        Z_reduce = [Z(1:m,:); Z(tmp_ind-(m-1):tmp_ind,:)];
        var_vector = diag(Z_reduce*Z_reduce')/trace(Z_reduce*Z_reduce');
        fp(:,k) = log(var_vector);
    end
    
    [check, prediction] = myClassifier(fp,Mr,Ml,Qr,Ql);
    
    prediction0 = true_y(iter);
    
    for tmpp = iter:iter+chunk-1
        predictions(:,tmpp) = [prediction];
    end
    iter = iter + chunk*0.8;
end

end
% ----------------------------------------------------------------------- %
%                               EOF
% ----------------------------------------------------------------------- %
