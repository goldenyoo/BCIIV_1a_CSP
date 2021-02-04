% ----------------------------------------------------------------------- %
%    File_name: Calib.m
%    Programmer: Seungjae Yoo                             
%                                           
%    Last Modified: 2020_02_04                           
%                                                            
 % ----------------------------------------------------------------------- %
function [M0,M12,M1,M2,Q0,Q12,Q1,Q2,P_0_vs_12,P_1_vs_2] = Calib(answer,ref)

% Input parameters
data_label = string(answer(1,1));   
m = double(string(answer(2,1))); % feature vector will have length (2m)
low_f = double(string(answer(3,1))); % Low cutoff freq
high_f = double(string(answer(4,1))); % High cutoff freq

referencing = double(string(answer(5,1))); % Non(0), CAR(1), LAP(2)
order = double(string(answer(6,1))); % Filter order

% Load file

FILENAME = strcat('C:\Users\유승재\Desktop\Motor Imagery EEG data\BCICIV_1_mat\BCICIV_calib_ds1',data_label,'.mat');
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

%% Calculate spatial filter
a = 1; b = 1; d = 1;
C_0 = zeros(size(cnt_c,1)); C_1 = zeros(size(cnt_c,1)); C_2 = zeros(size(cnt_c,1));

% Training only for training data set
for i = 1:length(mrk.pos)
    
    % One trial data
    E = cnt_c(:,mrk.pos(1,i):mrk.pos(1,i)+chunk);
        
    % Covariance 연산
    C = E*E'/ trace( E*E');
    
    % According to its class, divide calculated covariance
    if mrk.y(1,i) == 1
        C_1 = C_1+C;
        a = a+1;
    else
        C_2 = C_2+C;
        b = b+1;
    end
    
     % One trial data
    E = cnt_c(:,mrk.pos(1,i)+chunk:mrk.pos(1,i)+2*chunk);
        
    % Covariance 연산
    C = E*E'/ trace( E*E');
    
    C_0 = C_0 + C;
    d = d+1;    
end

% Average covariance of each class
C_0 = C_0/(d-1);
C_1 = C_1/(a-1);
C_2 = C_2/(b-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Class 1 vs Class 2%%%%%%%%%%%%%%%%%%%%%%%%%
% composite covariace
C_c = C_1 + C_2;

% EVD for composite covariance
[V, D] = eig(C_c);

% sort eigen vector with descend manner
[d, ind] = sort(abs(diag(D)),'descend');
D_new = diag(d);
V_new = V(:,ind);

% whitening transformation
whiten_tf = V_new*D_new^(-0.5);
W = whiten_tf';

% Apply whitening to each averaged covariances
Sa = W*C_1*W';
Sb = W*C_2*W';

% EVD for transformed covariance
[U, phsi] = eig(Sa,Sb);

% sort
[d, ind] = sort(abs(diag(phsi)),'descend');
phsi_new = diag(d);
U_new = U(:,ind);

% Total Projection matrix,   Z = P'*X
P_1_vs_2 = (U_new'*W)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Class 0 vs Class 1,2%%%%%%%%%%%%%%%%%%%%%%%
% composite covariace
C_c = C_0 + 0.5*(C_1+C_2);

% EVD for composite covariance
[V, D] = eig(C_c);

% sort eigen vector with descend manner
[d, ind] = sort(abs(diag(D)),'descend');
D_new = diag(d);
V_new = V(:,ind);

% whitening transformation
whiten_tf = V_new*D_new^(-0.5);
W = whiten_tf';

% Apply whitening to each averaged covariances
Sa = W*C_0*W';
Sb = W*(0.5*(C_1+C_2))*W';

% EVD for transformed covariance
[U, phsi] = eig(Sa);

% sort
[d, ind] = sort(abs(diag(phsi)),'descend');
phsi_new = diag(d);
U_new = U(:,ind);

% Total Projection matrix,   Z = P'*X
P_0_vs_12 = (U_new'*W)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate feature vector
fp_0 = [];
fp_12 = [];
fp_1 = [];
fp_2 = [];

for i = 1:length(mrk.pos)
    
    % One trial data
    E = cnt_c(:,mrk.pos(1,i):mrk.pos(1,i)+chunk);
    
    % Project data using calculated spatial filter
    Z = P_1_vs_2'*E;
    
    % Feature vector
    tmp_ind = size(Z,1);
    Z_reduce = [Z(1:m,:); Z(tmp_ind-(m-1):tmp_ind,:)];
    var_vector = var(Z_reduce,0,2)';
    var_vector = (1/sum(var_vector))*var_vector;
    
    fp = log(var_vector);
    fp = fp';
    
    if mrk.y(1,i) == 1
        fp_1 = [fp_1 fp];
    else
        fp_2 = [fp_2 fp];
    end
    
    Z = P_0_vs_12'*E;
    
    % Feature vector
    tmp_ind = size(Z,1);
    Z_reduce = [Z(1:m,:); Z(tmp_ind-(m-1):tmp_ind,:)];
    var_vector = var(Z_reduce,0,2)';
    var_vector = (1/sum(var_vector))*var_vector;
    
    fp = log(var_vector);
    fp = fp';
    
    fp_12 = [fp_12 fp];
    
    % One trial data
    E = cnt_c(:,mrk.pos(1,i)+chunk:mrk.pos(1,i)+2*chunk);
    
    % Project data using calculated spatial filter
    Z = P_0_vs_12'*E;
    
    % Feature vector
    tmp_ind = size(Z,1);
    Z_reduce = [Z(1:m,:); Z(tmp_ind-(m-1):tmp_ind,:)];
    var_vector = var(Z_reduce,0,2)';
    var_vector = (1/sum(var_vector))*var_vector;
    
    fp = log(var_vector);
    fp = fp';
    
    fp_0 = [fp_0 fp];    
end
M0 = mean(fp_0,2);
M12 = mean(fp_12,2);
M1 = mean(fp_1,2);
M2 = mean(fp_2,2);

Q0 = zeros(2*m);
for i = 1:length(fp_0)
    tmp = (fp_0(:,i) - M0)*(fp_0(:,i) - M0)';
    Q0 = Q0 + tmp;
end
Q0 = (1/(length(fp_0)-1))*Q0;

Q12 = zeros(2*m);
for i = 1:length(fp_12)
    tmp = (fp_12(:,i) - M12)*(fp_12(:,i) - M12)';
    Q12 = Q12 + tmp;
end
Q12 = (1/(length(fp_12)-1))*Q12;


Q1 = zeros(2*m);
for i = 1:length(fp_1)
    tmp = (fp_1(:,i) - M1)*(fp_1(:,i) - M1)';
    Q1 = Q1 + tmp;
end
Q1 = (1/(length(fp_1)-1))*Q1;

Q2 = zeros(2*m);
for i = 1:length(fp_2)
    tmp = (fp_2(:,i) - M2)*(fp_2(:,i) - M2)';
    Q2 = Q2 + tmp;
end
Q2 = (1/(length(fp_2)-1))*Q2;

end
% ----------------------------------------------------------------------- %
%                               EOF
% ----------------------------------------------------------------------- %
