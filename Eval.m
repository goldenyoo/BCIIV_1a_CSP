% ----------------------------------------------------------------------- %
%    File_name: Eval.m
%    Programmer: Seungjae Yoo                             
%                                           
%    Last Modified: 2020_01_27                           
%                                                            
 % ----------------------------------------------------------------------- %
function predictions = Eval(answer,fp_l,fp_r,P,ref)
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
        cnt_c = cnt([27 29 31 44 46 50 52 54],:); % Exclude electrode (AF3, AF4, O1, O2, PO1, PO2)        
        Means = (1/size(cnt,1))*sum(cnt);
        for i = 1 : size(cnt_c,1)
            cnt_c(i,:) = cnt_c(i,:) - Means; % CAR
        end
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
X_train = [fp_l';fp_r'];
Y_train = [repmat(-1,length(fp_l),1); repmat(1,length(fp_l),1)];

Mdl = fitcnb(X_train, Y_train);
% Mdl = fitcdiscr(X_train, Y_train);
% Mdl = fitcsvm(X_train, Y_train);
% Mdl = fitcsvm(X_train, Y_train,'Crossval','on');
% CVMdl = crossval(Mdl);


%%
% f1 = figure;
% f2 = figure;
score = [];
predictions = [];
checks = [];


iter = 1;
chunk = 200;
while iter + chunk <= size(cnt_c,2)
%     if rem(iter,10000)== 0 
%         fprintf("%d",iter)
%     end
    E = cnt_c(:, iter:iter+chunk);
    Z = P'*E;
    
     % Feature vector
    tmp_ind = size(Z,1);
    Z_reduce = [Z(1:m,:); Z(tmp_ind-(m-1):tmp_ind,:)];
   
    var_vector = var(Z_reduce,0,2)';
    var_vector = (1/sum(var_vector))*var_vector;
        
    fp = log(var_vector);
       
       
    % Run classifier
    prediction = predict(Mdl,fp);
    
    predictions = [predictions prediction];
    
    iter = iter+1;
       
end
predictions = [predictions repmat(0,1,size(cnt_c,2)- iter + 1)];


end
% ----------------------------------------------------------------------- %
%                               EOF
% ----------------------------------------------------------------------- %
