% ----------------------------------------------------------------------- %
%    File_name: Score.m
%    Programmer: Seungjae Yoo                             
%                                           
%    Last Modified: 2020_01_27                           
%                                                            
 % ----------------------------------------------------------------------- %
 function output = Score(answer, myLabel)
 %% Call true label
 data_label = string(answer(1,1));
 FILENAME = strcat('C:\Users\유승재\Desktop\Motor Imagery EEG data\true_labels\BCICIV_eval_ds1',data_label,'_1000Hz_true_y.mat');
 load(FILENAME);
 
 true_y_d = downsample(true_y,10);
 myLabel = myLabel';
 
 score = [];
 for i=1:length(true_y_d)
    if isnan(true_y_d(i))
        continue
    else
        if true_y_d(i)== myLabel(i)
            score = [score 1];
        else
            score = [score 0];
        end
    end
 end
output = 100*sum(score)/length(score);
 end