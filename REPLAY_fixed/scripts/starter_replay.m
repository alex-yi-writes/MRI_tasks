% STARTER SCRIPT : REPLAY
%  Makte sure you ru3n each cell individually to avoid MATLAB throwing
%  errors.
%  In emergencytt contact Alex: 
%              +49 162 713 20 4943431
%              yeo-jin.yi.15@ucl.ac.uk

%%  work log

%   23_09_2020      created the script


%% experiment preparation
%  always run this cell before starting an experiment.434

clear all; close all; clc;
rand('state',sum(100*clock));

% taskpath       = 'C:\Users\presentation\Desktop\REPLAY\';
taskpath       = 'C:\Users\direk\Documents\MATLAB\REPLAY\';
savepath_behav = [taskpath 'data\behav\'];
logpath        = [taskpath 'data\logs\'];
addpath(genpath([taskpath '\scripts\']))

% enter experiment detailst
input_prompt = {'ID (e.g. 4001)';'GROUP (1 = YAs, 2= OAs)';'0=encoding, 1=retrieval1, 2=retrieval2'; ...
    'block (1,2=enc / 1,2,3,4=ret)'};
defaults     = {'9999','1','1','0',' '};
input_answer = inputdlg(input_prompt, 'Input experiment details. NO SPACES!', 1, defaults);

ID           = str2num(input_answer{1,1});
AgeGroup     = str2num(input_answer{2,1});
TaskType     = str2num(input_answer{3,1});
BlockNum     = str2num(input_answer{4,1});
path_ECHO    = taskpath;

%% Start Experiment
% diary on
try
datetime
catch
datestr(now)
end
    

if TaskType == 0 % encoding
    
    filestrings         = 'REPLAY7Tenc';
    behavfilename       = [num2str(ID) '_enc' num2str(BlockNum) '.mat'];
    
    [dat] = runREPLAYencoding_cg(ID,AgeGroup,BlockNum,behavfilename,path_ECHO);
    
    % save results
    save([savepath_behav behavfilename],'dat')
    
%     if BlockNum==2
%         %%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE STIMLIST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         clear tmplist
%         makeRetrievalList(taskpath,ID)
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     end
    
    
elseif TaskType == 1 % retrieval 1
    
    filestrings         = 'REPLAY7Tretrieval1';
    behavfilename       = [num2str(ID) '_ret1' num2str(BlockNum) '.mat'];
    
    [dat] = runREPLAYretieval_1_ptb(ID,AgeGroup,BlockNum,behavfilename,path_ECHO);
    
    % save results
    save([savepath_behav behavfilename],'dat')
    
elseif TaskType == 2
    
    filestrings         = 'REPLAY7Tretrieval2';
    behavfilename       = [num2str(ID) '_ret2' num2str(BlockNum) '.mat'];
    
    [dat] = runREPLAYretieval_2_cg(ID,AgeGroup,BlockNum,behavfilename,path_ECHO);
    
    % save results
    save([savepath_behav behavfilename],'dat')
    
else
    error('Something''s wrong... better call Alex: +49 162 713 20 49')
    
end


% try
%     cd(savepath_behav)
%     mail = 'sfb779.behavdata@gmail.com';
%     password = '*caldzn3';
%     port = '465';
%     setpref('Internet','SMTP_Server','smtp.gmail.com');
%     setpref('Internet','E_mail',mail);
%     setpref('Internet','SMTP_Username',mail);
%     setpref('Internet','SMTP_Password',password);
%     props = java.lang.System.getProperties;
%     props.setProperty('mail.smtp.auth','true');
%     props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
%     props.setProperty('mail.smtp.socketFactory.port',port);
%     try
%         sendmail(mail,['Behav_' num2str(ID) '_' filestrings],char(datetime),behavfilename)
%     catch
%         sendmail(mail,['Behav_' num2str(ID) '_' filestrings],char(datestr(now)),behavfilename)
%     end
% catch
%     warning('data not sent! move them manually')
% end

cd(logpath)
try
datetime
catch
datestr(now)
end
    
disp('****************************************')

% diary off

%% troubleshoot

% document when something comes up