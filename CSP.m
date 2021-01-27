% ----------------------------------------------------------------------- %
%    File_name: CSP.m
%    Programmer: Seungjae Yoo                             
%                                           
%    Last Modified: 2020_01_26                            
%                                                            
 % ----------------------------------------------------------------------- %

%% Call raw data
close all
clear all

% Ask user for input parameters
prompt = {'Data label: ', 'Feature vector length: ', 'low cutoff freq', 'high cutoff freq'};
dlgtitle = 'Input';
dims = [1 50];
definput = {'a', '3','8','30'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
% Error detection
if isempty(answer), error("Not enough input parameters."); end



%% 
% Input parameters
ref_method = [0 1 2];
filt_ord = [10 15 20 25 30];

%% 


for i = 1:length(ref_method)
    for j = 1:length(filt_ord)
        answer(5,1) = {ref_method(i)};
        answer(6,1) = {filt_ord(j)};
        fprintf('Re-referencing: %d',ref_method(i));
        fprintf(' filter_order: %d',filt_ord(j));
        [Mr,Ml,Qr,Ql,P] = Calib(answer);
        tmp = Eval(answer,Mr,Ml,Qr,Ql,P);
        output(i,j) = tmp;
        fprintf('----> score: %f\n',tmp);
    end
    fprintf('\n');
end

fprintf('Data_Label: %s\n',string(answer(1,1)));
fprintf('BPF cutoff freq: %s ~ %s\n',string(answer(3,1)),string(answer(4,1)));
disp(output);
% ----------------------------------------------------------------------- %
%                               EOF
% ----------------------------------------------------------------------- %
