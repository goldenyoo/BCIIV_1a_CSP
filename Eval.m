% ----------------------------------------------------------------------- %
%    File_name: Eval.m
%    Programmer: Seungjae Yoo                             
%                                           
%    Last Modified: 2020_02_04                           
%                                                            
 % ----------------------------------------------------------------------- %
function myLabel = Eval(answer,M0,M12,M1,M2,Q0,Q12,Q1,Q2,P_0_vs_12,P_1_vs_2,ref)
data_label = string(answer(1,1));   
m = double(string(answer(2,1))); % feature vector will have length (2m)
low_f = double(string(answer(3,1))); % Low cutoff freq
high_f = double(string(answer(4,1))); % High cutoff freq
referencing = double(string(answer(5,1))); % Non(0), CAR(1), LAP(2)
order = double(string(answer(6,1))); % Filter order



%% 
% Load file

FILENAME = strcat('C:\Users\유승재\Desktop\Motor Imagery EEG data\BCICIV_1_mat\BCICIV_eval_ds1',data_label,'.mat');
chunk = 400;
fs = 100;

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
        cnt_c = cnt(3:55,:); % Exclude electrode (AF3, AF4, O1, O2, PO1, PO2)        
        Means = (1/size(cnt_c,1))*sum(cnt_c);
        for i = 1 : size(cnt_c,1)
            cnt_c(i,:) = cnt_c(i,:) - Means; % CAR
        end
     % LAP   
    elseif referencing == 2 
        cnt_n = myLAP(cnt,nfo); % Laplacian
        cnt_c = cnt_n(3:55,:); % Exclude electrode (AF3, AF4, O1, O2, PO1, PO2)
    end
else
        %%% Calculate differential voltage
    for i = 1 : size(cnt,1)
        cnt(i,:) = cnt(i,:) - cnt(ref,:);
    end
    
    cnt_c = cnt(3:55,:); % Exclude electrode (AF3, AF4, O1, O2, PO1, PO2)
end

clear cnt cnt_n

%% 
%BPF Design
bpFilt = designfilt('bandpassfir','FilterOrder',order, ...
         'CutoffFrequency1',low_f,'CutoffFrequency2',high_f, ...
         'SampleRate',fs);

% Apply BPF
for i = 1:size(cnt_c,1)
    cnt_c(i,:) = filtfilt(bpFilt, cnt_c(i,:));
end

%% 
myLabel = zeros(1,size(cnt_c,2));

i=1;
while i+150 < size(cnt_c,2)
    E = cnt_c(:,i:i+150);
    Z = P_0_vs_12'*E;
    
    % Feature vector
    tmp_ind = size(Z,1);
    Z_reduce = [Z(1:m,:); Z(tmp_ind-(m-1):tmp_ind,:)];
   
    var_vector = var(Z_reduce,0,2)';
    var_vector = (1/sum(var_vector))*var_vector;
        
    fp = log(var_vector);
    fp = fp';
    
    [check, prediction] = myClassifier(fp,M0,M12,Q0,Q12);
    if prediction == 1
        myLabel(1,i:i+150) = 0;
    else
        Z = P_1_vs_2'*E;
        
        % Feature vector
        tmp_ind = size(Z,1);
        Z_reduce = [Z(1:m,:); Z(tmp_ind-(m-1):tmp_ind,:)];
        
        var_vector = var(Z_reduce,0,2)';
        var_vector = (1/sum(var_vector))*var_vector;
        
        fp = log(var_vector);
        fp = fp';
        
        [check, prediction1] = myClassifier(fp,M1,M2,Q1,Q2);
        if prediction1 == 1
            myLabel(1,i:i+150) = 1;
        else
            myLabel(1,i:i+150) = -1;
        end
    end
    i = i+120;
end

end
% ----------------------------------------------------------------------- %
%                               EOF
% ----------------------------------------------------------------------- %
