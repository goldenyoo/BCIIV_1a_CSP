% ----------------------------------------------------------------------- %
%    File_name: CSP.m
%    Programmer: Seungjae Yoo                             
%                                           
%    Last Modified: 2020_01_27                            
%                                                            
 % ----------------------------------------------------------------------- %

%% Get input parameter from user
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

%% Conditions
% Rereferencing method 
ref_method = [0 ]; % Non(0), CAR(1), LAP(2)

% Filter order
filt_ord = [20 25 30];

% Reference electrode number
ref = 29;        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Change

%% CSP 
for i = 1:length(ref_method)
    for j = 1:length(filt_ord)
        answer(5,1) = {ref_method(i)};
        answer(6,1) = {filt_ord(j)};
        fprintf('Re-referencing: %d',ref_method(i));
        fprintf(' filter_order: %d',filt_ord(j));
        [M0,M12,M1,M2,Q0,Q12,Q1,Q2,P_0_vs_12,P_1_vs_2] = Calib(answer,ref);
        myLabel = Eval(answer,M0,M12,M1,M2,Q0,Q12,Q1,Q2,P_0_vs_12,P_1_vs_2,ref);
        score = Score(answer, myLabel);
        fprintf(' ----> score: %f\n',score);
    end
    fprintf('\n');
end

%% Output present

fprintf('Data_Label: %s\n',string(answer(1,1)));
fprintf('BPF cutoff freq: %s ~ %s\n',string(answer(3,1)),string(answer(4,1)));

% ----------------------------------------------------------------------- %
%                               EOF
% ----------------------------------------------------------------------- %
