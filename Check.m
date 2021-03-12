function [score, total , tmp] = Check(answer,predictions)
data_label = string(answer(1,1));   

%% Call true label
FILENAME = strcat('C:\Users\유승재\Desktop\Motor Imagery EEG data\true_labels\BCICIV_eval_ds1',data_label,'_1000Hz_true_y.mat');
load(FILENAME);

true_y = downsample(true_y,10);

total = 0;
score = 0;
mse = 0;

score_0 = 0;
total_0 = 0;

score_12 = 0;
total_12 = 0;
for i = 1:size(true_y,1)
    if isnan(true_y(i))        
        continue
    elseif true_y(i)==0
%         total = total + 1;
%         total_0 = total_0 +1;
%         if predictions(i) == 0
%             score = score + 1;
%             score_0 = score_0 + 1;
%         else
%             mse = mse + (true_y(i)-predictions(i))^2;
%         end
    continue
        
    elseif true_y(i)== -1
        total = total + 1;
       
        if predictions(i) == -1
            score = score + 1;
                 
        else
%             mse = mse + (true_y(i)-predictions(i))^2;
        end
    elseif true_y(i)== 1
         total = total + 1;
         
        if predictions(i) == 1
            score = score + 1;
            
        else
%             mse = mse + (true_y(i)-predictions(i))^2;
        end
    end
end
% mse = mse/total;

tmp = [true_y'; predictions];

end