%% MEMORY TEST : immediate and delayed recall
%  In emergency contact Alex:
%              +49 162 713 20 49
%              yeo-jin.yi.15@ucl.ac.uk
% See also StarterScript_MRIpilot, runPracticeTask_MRIpilot,
%       runMainTask_MRIpilot

%%  work log

%   10_10_2018 created the script

%% START

function [dat] = runMemoryTest_3Tpilot_cg_mid(ID,TaskType,EyeFlag,MRIFlag,CnterBal,eyetracker_filename,Schedules,DayNumber)

global dat eyeflag eyetracker_filename

%% basic setups

% initialize eyetracking
eyeflag=EyeFlag;
fname_eyetracker = eyetracker_filename;

if eyeflag==1
    
    Eyelink('Initialize')
    
    %%%%%%%%%%%% INITIALISE CALIBRATION AND VALIDATION %%%%%%%%%%%%%%%
    PsychDefaultSetup(1);
    try
        fprintf('EyelinkToolbox Example\n\n\t');
        dummymode=0;       % set to 1 to initialize in dummymode (rather pointless for this example though)
        
        % STEP 1
        % Open a graphics window on the main screen
        % using the PsychToolbox's Screen function.
        screenNumber=max(Screen('Screens'));
        window=Screen('OpenWindow', screenNumber);
        
        % STEP 2
        % Provide Eyelink with details about the graphics environment
        % and perform some initializations. The information is returned
        % in a structure that also contains useful defaults
        % and control codes (e.g. tracker state bit and Eyelink key values).
        el=EyelinkInitDefaults(window);
        
        % Disable key output to Matlab window:
        ListenChar(2);
        
        % STEP 3
        % Initialization of the connection with the Eyelink Gazetracker.
        % exit program if this fails.
        if ~EyelinkInit(dummymode, 1)
            fprintf('Eyelink Init aborted.\n');
            %         cleanup;  % cleanup function
            return;
        end
        
        [v vs]=Eyelink('GetTrackerVersion');
        fprintf('Running experiment on a ''%s'' tracker.\n', vs );
        
        % make sure that we get gaze data from the Eyelink
        Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');
        
        % open file to record data to
%         edffile = 
        edfFile=fname_eyetracker;
        Eyelink('Openfile', edfFile);
        
        % STEP 4
        % Calibrate the eye tracker
        EyelinkDoTrackerSetup(el); % changed this bit to see what happens to the data
        
        
        % STEP 7
        % finish up: stop recording eye-movements,
        % close graphics window, close data file and shut down tracker
        %     Eyelink('StopRecording');
        %     Eyelink('CloseFile');
        % download data file
        %     try
        %         fprintf('Receiving data file ''%s''\n', edfFile );
        %         status=Eyelink('ReceiveFile');
        %         if status > 0
        %             fprintf('ReceiveFile status %d\n', status);
        %         end
        %         if 2==exist(edfFile, 'file')
        %             fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
        %         end
        %     catch rdf
        %         fprintf('Problem receiving data file ''%s''\n', edfFile );
        %         rdf;
        %     end
        
        %     cleanup;
        
    catch
        %this "catch" section executes in case of an error in the "try" section
        %above.  Importantly, it closes the onscreen window if its open.
        %     cleanup;
        psychrethrow(psychlasterror);
    end %try..catch.
    
    Screen('CloseAll'); sca;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     eyelink_calval(dat.eyetracker_filename)
    %     Eyelink('Openfile',   dat.eyetracker_filename);
    
end

% set pathss
% paths.parent   = '\\fs-md\users\yiy\Documents\DOPE\';
paths.parent   = '\\fs-md\users\yiy\Documents\DOPE\';
paths.stim     = [ paths.parent 'pilot_3T\stim_mirrored\' ];  %%%%% changed
paths.fb       = [ paths.parent 'pilot_3T\feedback\' ];
paths.pupil    = [ paths.parent 'pilot_3T\raw\pupil\' ];
paths.behav    = [ paths.parent 'pilot_3T\raw\behav\' ];
paths.cogent   = 'C:\Studies\Grid_experiment2\Cogent2000v1.32';
load([paths.stim 'Memorability_API_midrange_20190320.mat' ]) %%%%% changed
load([paths.behav num2str(ID) '_' num2str(DayNumber) '.mat'])
addpath(genpath(paths.cogent))
% addpath(genpath('\\fs-md\users\yiy\Documents\MATLAB\Cogent2000v1.33\'))
addpath(genpath('\\fs-md\users\yiy\Documents\MATLAB\Cogent2000v1.33\'))


% create data structure
if TaskType == 3 % if immediate recall
    if DayNumber == 1
        dat.day1.memorytest  = [];
        dat.day1.memorytest.immediate = [];
        dat.day1.memorytest.immediate.ETinfo  = eyetracker_filename;
        if Schedules == 1
            load([ paths.stim 'ImmMem_GA_designstruct_D1_ver1.mat' ])
        elseif Schedules == 2
            load([ paths.stim 'ImmMem_GA_designstruct_D1_ver2.mat' ])
        end
    elseif DayNumber == 2
        dat.day2.memorytest  = [];
        dat.day2.memorytest.immediate = [];
        dat.day2.memorytest.immediate.ETinfo  = eyetracker_filename;
        if Schedules == 1
            load([ paths.stim 'ImmMem_GA_designstruct_D2_ver1.mat' ])
        elseif Schedules == 2
            load([ paths.stim 'ImmMem_GA_designstruct_D2_ver2.mat' ])
        end
    end
elseif TaskType == 4 % if delayed recall
    if DayNumber == 1
        dat.day1.memorytest.delayed = [];
        dat.day1.memorytest.delayed.ETinfo  = eyetracker_filename;
    elseif DayNumber == 2
        dat.day2.memorytest.delayed = [];
        dat.day2.memorytest.delayed.ETinfo  = eyetracker_filename;
    end
end

% scanner preparation
if MRIFlag == 1
    
    if DayNumber == 1
        dat.day1.maintask.MRIinfo   = [];
        dat.day1.maintask.MRIinfo.scanPort = 1;
        dat.day1.maintask.MRIinfo.dummy    = 5;        % no. of dummy vols
        dat.day1.maintask.MRIinfo.nslice   = 51;       % no. of slices
        dat.day1.maintask.MRIinfo.TE       = 32;       % time for each slice in msec
        dat.day1.maintask.MRIinfo.TR       = 3.6; % MRIinfo.nslice * MRIinfo.TE/1000; % time for volume in msec
%         dat.day1.maintask.MRIinfo.nvols    = ceil(dat.day1.maintask.MRIinfo.dummy + 5 + ...  % no. of vols per run (5 dummy + 400 + 5 post)
%             (options.scanblocklength * options.meanTrialLength) / dat.day1.maintask.MRIinfo.TR);
%         dat.day1.maintask.MRIinfo.total_slices = dat.day1.maintask.MRIinfo.nvols*dat.day1.maintask.MRIinfo.nslice;    % per run!
        %     [MRIinfo] = runScannerSetup_MRIpilot; % redundant
        
    elseif DayNumber == 2
        dat.day2.maintask.MRIinfo   = [];
        dat.day2.maintask.MRIinfo.scanPort = 1;
        dat.day2.maintask.MRIinfo.dummy    = 5;        % no. of dummy vols
        dat.day2.maintask.MRIinfo.nslice   = 51;       % no. of slices
        dat.day2.maintask.MRIinfo.TE       = 32;       % time for each slice in msec
        dat.day2.maintask.MRIinfo.TR       = 3.6; % MRIinfo.nslice * MRIinfo.TE/1000; % time for volume in msec
%         dat.day2.maintask.MRIinfo.nvols    = ceil(dat.day2.maintask.MRIinfo.dummy + 5 + ...  % no. of vols per run (5 dummy + 400 + 5 post)
%             (options.scanblocklength * options.meanTrialLength) / dat.day2.maintask.MRIinfo.TR);
%         dat.day2.maintask.MRIinfo.total_slices = dat.day1.maintask.MRIinfo.nvols*dat.day2.maintask.MRIinfo.nslice;    % per run!
        %     [MRIinfo] = runScannerSetup_MRIpilot; % redundant
        
    end
end


%% experimental setups

% counterbalancing
if CnterBal == 1 % if land/private = reward
    if DayNumber == 1 % indoors dataset
        rewcond = 3; puncond = 4;
    elseif DayNumber ==2 % outdoors
        rewcond = 1; puncond = 2;
    end;
elseif CnterBal == 2 % if city/public == reward
    if DayNumber == 1 % indoors
        rewcond = 4; puncond = 3;
    elseif DayNumber ==2 % outdoors
        rewcond = 2; puncond = 1;
    end;
end

% timings
config.timing.outcome  = 2500;
tmp = [(design_struct.eventlist(:,4)+design_struct.eventlist(:,5)-2)']*1000;
config.timing.fixation1= horzcat(tmp(1:60), [NaN], tmp(61:120), [NaN], ...
    tmp(121:end)); clear tmp
config.timing.cue      = 2500;
config.timing.response = 2000;

% config.timing.outcome  = 2500;
% fixjit1 = repmat(1000:250:2000,1,48); fixjit2 = repmat(1000:250:2000,1,48);
% config.timing.fixation1= fixjit1(randperm(length(fixjit1)));
% config.timing.fixation2= fixjit2(randperm(length(fixjit2)));
% config.timing.cue      = 2500;
% config.timing.resp     = 2000;

% --- for breezeing through --- %
% config.timing.fixation1= zeros(1,240);
% config.timing.fixation2= zeros(1,240);
% config.timing.cue      = 1;
% config.timing.outcome  = 1;
% ----------------------------- %

% number of trials
config.numtrials.total      = 176;

if DayNumber == 1
    config.numtrials.old_total= 88;
    config.numtrials.old_Rew  = 66;
    config.numtrials.old_Pun  = 22;

    config.numtrials.new_total= 88;
    config.numtrials.new_Rew  = 66;
    config.numtrials.new_Pun  = 22;
    
elseif DayNumber == 2
    config.numtrials.old_total= 88;
    config.numtrials.old_Rew  = 22;
    config.numtrials.old_Pun  = 66;

    config.numtrials.new_total= 88;
    config.numtrials.new_Rew  = 22;
    config.numtrials.new_Pun  = 66;
    
end


%% roster

if DayNumber == 1
    maintask_fname          = dat.day1.maintask.config.stim.fname;
%     immediate_stimlist_old  = dat.day1.memorytest.immediate.config.stim.stimlist_old;
%     immediate_stimlist_new  = dat.day1.memorytest.immediate.config.stim.stimlist_new;
elseif DayNumber == 2
    maintask_fname          = dat.day2.maintask.config.stim.fname;
%     immediate_stimlist_old  = dat.day2.memorytest.immediate.config.stim.stimlist_old;
%     immediate_stimlist_new  = dat.day2.memorytest.immediate.config.stim.stimlist_new;
end

% ------------------- filenames setup: OLD stim list ------------------- %

if TaskType == 3 % if immediate recall
    
    c1 = 0; c2 = 0; c3 = 0; c4 = 0;
    stim_pun_o = []; stim_rew_o = [];
    for i = 1:size(maintask_fname(:,1),1)
        
        % if punished
        if cell2mat(maintask_fname(i,4)) == puncond
            c1 = c1+1;
                stim_pun_o{c1,1} = maintask_fname{i,1}; % filename
                stim_pun_o{c1,2} = maintask_fname{i,4}; % category
                stim_pun_o{c1,3} = maintask_fname{i,3}; % memorability 
            
        % if rewarded
        elseif cell2mat(maintask_fname(i,4)) == rewcond
                c3 = c3+1;
                stim_rew_o{c3,1} = maintask_fname{i,1};
                stim_rew_o{c3,2} = maintask_fname{i,4};
                stim_rew_o{c3,3} = maintask_fname{i,3};
        end
    end % close immediate recall stimlist loop
    
    % immediate recall old stimuli list (1 = old, 2 = new)
    fname_ImOld      = cellstr(vertcat(stim_rew_o{1:config.numtrials.old_Rew,1}, stim_pun_o{1:config.numtrials.old_Pun,1})); % filenames
    fname_ImOld(:,2) = vertcat(stim_rew_o(1:config.numtrials.old_Rew,2), stim_pun_o(1:config.numtrials.old_Pun,2)); % categories
    fname_ImOld(:,3) = vertcat(stim_rew_o(1:config.numtrials.old_Rew,3), stim_pun_o(1:config.numtrials.old_Pun,3)); % memorability scores
    fname_ImOld(:,4) = num2cell(ones(size(fname_ImOld,1),1)); % old/new marker: OLD
    config.stim.stimlist_old = fname_ImOld;
    
    % extra
    fname_ImOld_rew      = cellstr(vertcat(stim_rew_o{1:config.numtrials.new_Rew,1}));
    fname_ImOld_rew(:,2) = vertcat(stim_rew_o(1:config.numtrials.new_Rew,2));
    fname_ImOld_rew(:,3) = vertcat(stim_rew_o(1:config.numtrials.new_Rew,3));
    fname_ImOld_rew(:,4) = num2cell(ones(size(fname_ImOld_rew,1),1));
    
    fname_ImOld_pun      = cellstr(vertcat(stim_pun_o{1:config.numtrials.new_Pun,1}));
    fname_ImOld_pun(:,2) = vertcat(stim_pun_o(1:config.numtrials.new_Pun,2));
    fname_ImOld_pun(:,3) = vertcat(stim_pun_o(1:config.numtrials.new_Pun,3));
    fname_ImOld_pun(:,4) = num2cell(ones(size(fname_ImOld_pun,1),1));
    
    config.stim.stimlist_old_rew = fname_ImOld_rew;
    config.stim.stimlist_old_pun = fname_ImOld_pun;
    
elseif TaskType == 4 % if delayed recall
    
    if DayNumber == 1
        maintask_fname          = dat.day1.maintask.config.stim.fname;
        immediate_stimlist_old  = dat.day1.memorytest.immediate.config.stim.stimlist_old;
        immediate_stimlist_new  = dat.day1.memorytest.immediate.config.stim.stimlist_new;
    elseif DayNumber == 2
        maintask_fname          = dat.day2.maintask.config.stim.fname;
        immediate_stimlist_old  = dat.day2.memorytest.immediate.config.stim.stimlist_old;
        immediate_stimlist_new  = dat.day2.memorytest.immediate.config.stim.stimlist_new;
    end
    
    c1 = 0; c3 = 0; stim_pun_n = []; stim_rew_n = [];
    for i = 1:size(maintask_fname(:,1),1)
        
        % if punished & not used for the immediate test
        if cell2mat(maintask_fname(i,4)) == puncond && ...
                ~ismember(maintask_fname{i,1}, immediate_stimlist_old(:,1))
                c1 = c1+1;
                stim_pun_o{c1,1} = maintask_fname{i,1};
                stim_pun_o{c1,2} = maintask_fname{i,4};
                stim_pun_o{c1,3} = maintask_fname{i,3};            
        % if rewarded & not used for the immediate test
        elseif cell2mat(maintask_fname(i,2)) == rewcond && ...
                ~ismember(maintask_fname{i,1}, immediate_stimlist_old(:,1))
                c3 = c3+1;
                stim_rew_o{c3,1} = maintask_fname{i,1};
                stim_rew_o{c3,2} = maintask_fname{i,4};
                stim_rew_o{c3,3} = maintask_fname{i,3};
        end
        
    end % close delayed recall stimlist loop
    
    % delayed recall old stimuli list (1 = old, 2 = new)
    fname_DeOld      = cellstr(vertcat(stim_rew_o{1:config.numtrials.old_Rew,1}, stim_pun_o{1:config.numtrials.old_Pun,1})); % filenames
    fname_DeOld(:,2) = vertcat(stim_rew_o(1:config.numtrials.old_Rew,2), stim_pun_o(1:config.numtrials.old_Pun,2)); % categories
    fname_DeOld(:,3) = vertcat(stim_rew_o(1:config.numtrials.old_Rew,3), stim_pun_o(1:config.numtrials.old_Pun,3)); % memorability scores
    fname_DeOld(:,4) = num2cell(ones(size(fname_DeOld,1),1)); % old/new marker: OLD
    config.stim.stimlist_old = fname_DeOld;
    
    % extra
    fname_DeOld_rew      = cellstr(vertcat(stim_rew_o{1:config.numtrials.new_Rew,1}));
    fname_DeOld_rew(:,2) = vertcat(stim_rew_o(1:config.numtrials.new_Rew,2));
    fname_DeOld_rew(:,3) = vertcat(stim_rew_o(1:config.numtrials.new_Rew,3));
    fname_DeOld_rew(:,4) = num2cell(ones(size(fname_DeOld_rew,1),1));
    
    fname_DeOld_pun      = cellstr(vertcat(stim_pun_o{1:config.numtrials.new_Pun,1}));
    fname_DeOld_pun(:,2) = vertcat(stim_pun_o(1:config.numtrials.new_Pun,2));
    fname_DeOld_pun(:,3) = vertcat(stim_pun_o(1:config.numtrials.new_Pun,3));
    fname_DeOld_pun(:,4) = num2cell(ones(size(fname_DeOld_pun,1),1));
    
    config.stim.stimlist_old_rew = fname_DeOld_rew;
    config.stim.stimlist_old_pun = fname_DeOld_pun;
    
end % close immediate vs. delayed conditional



% ------------------- filenames setup: NEW stim list ------------------- %

if TaskType == 3 % if immediate recall
    
    c1 = 0; c3 = 0;
    stim_pun_n = []; stim_rew_n = []; %%%%% changed
    for i = 1:size(memorability_API_midrange_20190320,1) %%%%% changed
        
        % if punished & not used in the maintask
        if cell2mat(memorability_API_midrange_20190320(i,3)) == puncond && ...  %%%%% changed
                ~ismember(memorability_API_midrange_20190320{i,1},maintask_fname(:,1))  %%%%% changed
                c1 = c1+1;
                stim_pun_n{c1,1} = memorability_API_midrange_20190320{i,1}; %%%%% changed
                stim_pun_n{c1,2} = memorability_API_midrange_20190320{i,3}; %%%%% changed
                stim_pun_n{c1,3} = memorability_API_midrange_20190320{i,2}; %%%%% changed
                stim_pun_n{c1,4} = memorability_API_midrange_20190320{i,4}; %%%%% changed & added
        % if rewarded & not used in the maintask
        elseif cell2mat(memorability_API_midrange_20190320(i,3)) == rewcond && ...
                ~ismember(memorability_API_midrange_20190320{i,1},maintask_fname(:,1)) %%%%% changed
                c3 = c3+1;
                stim_rew_n{c3,1} = memorability_API_midrange_20190320{i,1}; %%%%% changed
                stim_rew_n{c3,2} = memorability_API_midrange_20190320{i,3}; %%%%% changed
                stim_rew_n{c3,3} = memorability_API_midrange_20190320{i,2}; %%%%% changed
                stim_rew_n{c3,4} = memorability_API_midrange_20190320{i,4}; %%%%% changed & added
        end
    end
    
    % immediate recall new stimuli list %%%%% changed
    fname_ImNew      = cellstr(vertcat(stim_rew_n{1:config.numtrials.new_Rew,1}, stim_pun_n{1:config.numtrials.new_Pun,1}));
    fname_ImNew(:,2) = vertcat(stim_rew_n(1:config.numtrials.new_Rew,2), stim_pun_n(1:config.numtrials.new_Pun,2));
    fname_ImNew(:,3) = vertcat(stim_rew_n(1:config.numtrials.new_Rew,3), stim_pun_n(1:config.numtrials.new_Pun,3));
    fname_ImNew(:,4) = num2cell(ones(size(fname_ImNew,1),1).*2);
    config.stim.stimlist_new = fname_ImNew;
    
    % extra %%%%% changed
    fname_ImNew_rew      = cellstr(vertcat(stim_rew_n{1:config.numtrials.new_Rew,1}));
    fname_ImNew_rew(:,2) = vertcat(stim_rew_n(1:config.numtrials.new_Rew,2));
    fname_ImNew_rew(:,3) = vertcat(stim_rew_n(1:config.numtrials.new_Rew,3));
    fname_ImNew_rew(:,4) = num2cell(ones(size(fname_ImNew_rew,1),1).*2);
    
    fname_ImNew_pun      = cellstr(vertcat(stim_pun_n{1:config.numtrials.new_Pun,1}));
    fname_ImNew_pun(:,2) = vertcat(stim_pun_n(1:config.numtrials.new_Pun,2));
    fname_ImNew_pun(:,3) = vertcat(stim_pun_n(1:config.numtrials.new_Pun,3));
    fname_ImNew_pun(:,4) = num2cell(ones(size(fname_ImNew_pun,1),1).*2);
    
    config.stim.stimlist_new_rew = fname_ImNew_rew;
    config.stim.stimlist_new_pun = fname_ImNew_pun;
    
    
elseif TaskType == 4 % if delayed recall
    
    c1 = 0; c3 = 0;
    stim_pun_n = []; stim_rew_n = []; %%%%% changed
    for i = 1:size(memorability_API_midrange_20190320,1) %%%%% changed
        
        if cell2mat(memorability_API_midrange_20190320(i,3)) == puncond && ... %%%%% changed
                ~ismember(memorability_API_midrange_20190320{i,1},maintask_fname(:,1)) && ... %%%%% changed
                ~ismember(memorability_API_midrange_20190320{i,1},immediate_stimlist_new(:,1)) %%%%% changed
                c1 = c1+1;
                stim_pun_n{c1,1} = memorability_API_midrange_20190320{i,1}; %%%%% changed
                stim_pun_n{c1,2} = memorability_API_midrange_20190320{i,3}; %%%%% changed
                stim_pun_n{c1,3} = memorability_API_midrange_20190320{i,2}; %%%%% changed
                stim_pun_n{c1,4} = memorability_API_midrange_20190320{i,4}; %%%%% changed & added
        elseif cell2mat(memorability_API_midrange_20190320(i,3)) == rewcond && ... %%%%% changed
                ~ismember(memorability_API_midrange_20190320{i,1},maintask_fname(:,1)) && ... %%%%% changed
                ~ismember(memorability_API_midrange_20190320{i,1},immediate_stimlist_new(:,1)) %%%%% changed
                c3 = c3+1;
                stim_rew_n{c3,1} = memorability_API_midrange_20190320{i,1}; %%%%% changed
                stim_rew_n{c3,2} = memorability_API_midrange_20190320{i,3}; %%%%% changed
                stim_rew_n{c3,3} = memorability_API_midrange_20190320{i,2}; %%%%% changed
                stim_rew_n{c3,4} = memorability_API_midrange_20190320{i,4}; %%%%% changed & added
        end
    end
    
    % delayed recall lure stimuli list %%%%% changed
    fname_DeNew      = cellstr(vertcat(stim_rew_n(1:config.numtrials.new_Rew,1), stim_pun_n(1:config.numtrials.new_Pun,1)));
    fname_DeNew(:,2) = vertcat(stim_rew_n(1:config.numtrials.new_Rew,2), stim_pun_n(1:config.numtrials.new_Pun,2));
    fname_DeNew(:,3) = vertcat(stim_rew_n(1:config.numtrials.new_Rew,3), stim_pun_n(1:config.numtrials.new_Pun,3));
    fname_DeNew(:,4) = num2cell(ones(size(fname_DeNew,1),1).*2);
    config.stim.stimlist_new = fname_DeNew;
    
    % extra %%%%% changed
    fname_DeNew_rew      = cellstr(vertcat(stim_rew_n{1:config.numtrials.new_Rew,1}));
    fname_DeNew_rew(:,2) = vertcat(stim_rew_n(1:config.numtrials.new_Rew,2));
    fname_DeNew_rew(:,3) = vertcat(stim_rew_n(1:config.numtrials.new_Rew,3));
    fname_DeNew_rew(:,4) = num2cell(ones(size(fname_DeNew_rew,1),1).*2);
    
    fname_DeNew_pun      = cellstr(vertcat(stim_pun_n{1:config.numtrials.new_Pun,1}));
    fname_DeNew_pun(:,2) = vertcat(stim_pun_n(1:config.numtrials.new_Pun,2));
    fname_DeNew_pun(:,3) = vertcat(stim_pun_n(1:config.numtrials.new_Pun,3));
    fname_DeNew_pun(:,4) = num2cell(ones(size(fname_DeNew_pun,1),1).*2);
    
    config.stim.stimlist_new_rew = fname_DeNew_rew;
    config.stim.stimlist_new_pun = fname_DeNew_pun;
end

% generate merged stimlist in presentation order
tmpvec       = 1:config.numtrials.total;
rewmat       = design_struct.eventlist(:,3)';
stimvec_temp = [];
orc = 0; nrc = 0; opc = 0; npc = 0;
for sc = 1:length(rewmat)
    if TaskType == 3 % immediate recall
        if rewmat(sc) == 11 % if old and rewarded
            orc = orc+1;
            stimvec_temp(1,sc) = tmpvec(sc); % trial order
            stimvec_temp(2,sc) = rewmat(sc); % reward/oldnew contingency
            stimvec_temp(3,sc) = fname_ImOld_rew{orc,2}; % category
            stimvec_flist{1,sc}= fname_ImOld_rew{orc,1};
            stimvec_flist{2,sc}= fname_ImOld_rew{orc,2};
            stimvec_flist{3,sc}= fname_ImOld_rew{orc,3};
        elseif rewmat(sc) == 21 % if old and punished
            opc = opc+1;
            stimvec_temp(1,sc) = tmpvec(sc);
            stimvec_temp(2,sc) = rewmat(sc);
            stimvec_temp(3,sc) = fname_ImOld_pun{opc,2};
            stimvec_flist{1,sc}= fname_ImOld_pun{opc,1};
            stimvec_flist{2,sc}= fname_ImOld_pun{opc,2};
            stimvec_flist{3,sc}= fname_ImOld_pun{opc,3};
        elseif rewmat(sc) == 12 % if new and rewarded
            nrc = nrc+1;
            stimvec_temp(1,sc) = tmpvec(sc);
            stimvec_temp(2,sc) = rewmat(sc);
            stimvec_temp(3,sc) = fname_ImNew_rew{nrc,2};
            stimvec_flist{1,sc}= fname_ImNew_rew{nrc,1};
            stimvec_flist{2,sc}= fname_ImNew_rew{nrc,2};
            stimvec_flist{3,sc}= fname_ImNew_rew{nrc,3};
        elseif rewmat(sc) == 22 % if new and punished
            npc = npc+1;
            stimvec_temp(1,sc) = tmpvec(sc);
            stimvec_temp(2,sc) = rewmat(sc);
            stimvec_temp(3,sc) = fname_ImNew_pun{npc,2};
            stimvec_flist{1,sc}= fname_ImNew_pun{npc,1};
            stimvec_flist{2,sc}= fname_ImNew_pun{npc,2};
            stimvec_flist{3,sc}= fname_ImNew_pun{npc,3};
        end
        
        
    elseif TaskType == 4 % delayed recall
        if rewmat(sc) == 11 % if old and rewarded
            orc = orc+1;
            stimvec_temp(1,sc) = tmpvec(sc); % trial order
            stimvec_temp(2,sc) = rewmat(sc); % reward/oldnew contingency
            stimvec_temp(3,sc) = fname_DeOld_rew{orc,2}; % category
            stimvec_flist{1,sc}= fname_DeOld_rew{orc,1};
            stimvec_flist{2,sc}= fname_DeOld_rew{orc,2};
            stimvec_flist{3,sc}= fname_DeOld_rew{orc,3};
        elseif rewmat(sc) == 21 % if old and punished
            opc = opc+1;
            stimvec_temp(1,sc) = tmpvec(sc);
            stimvec_temp(2,sc) = rewmat(sc);
            stimvec_temp(3,sc) = fname_DeOld_pun{opc,2};
            stimvec_flist{1,sc}= fname_DeOld_pun{opc,1};
            stimvec_flist{2,sc}= fname_DeOld_pun{opc,2};
            stimvec_flist{3,sc}= fname_DeOld_pun{opc,3};
        elseif rewmat(sc) == 12 % if new and rewarded
            nrc = nrc+1;
            stimvec_temp(1,sc) = tmpvec(sc);
            stimvec_temp(2,sc) = rewmat(sc);
            stimvec_temp(3,sc) = fname_DeNew_rew{nrc,2};
            stimvec_flist{1,sc}= fname_DeNew_rew{nrc,1};
            stimvec_flist{2,sc}= fname_DeNew_rew{nrc,2};
            stimvec_flist{3,sc}= fname_DeNew_rew{nrc,3};
        elseif rewmat(sc) == 22 % if new and punished
            npc = npc+1;
            stimvec_temp(1,sc) = tmpvec(sc);
            stimvec_temp(2,sc) = rewmat(sc);
            stimvec_temp(3,sc) = fname_DeNew_pun{npc,2};
            stimvec_flist{1,sc}= fname_DeNew_pun{npc,1};
            stimvec_flist{2,sc}= fname_DeNew_pun{npc,2};
            stimvec_flist{3,sc}= fname_DeNew_pun{npc,3};
        end
    end
end




stimlist_temp = vertcat(config.stim.stimlist_new,config.stim.stimlist_old);
stimvec       = randperm(config.numtrials.total);
for i2 = 1:length(stimvec)
    stimlist(i2,:) = stimlist_temp(stimvec(i2),:);
end

% insert breaks
clear tmp
tmp = horzcat(stimvec(1:60),777,stimvec(61:122),777,stimvec(123:end));
stimvec = tmp; clear tmp
config.stim.stimlist_all = stimlist;
config.stim.stimvec      = stimvec;

fprintf('\n\nTASK SETUP DONE')
fprintf('\nINITIALISING MEMORY TEST\n\n')
if TaskType == 3
    dat.day1.memorytest.immediate.config = config;
    fprintf('\n****************\nIMMEDIATE RECALL \n****************\n\n\n')
elseif TaskType == 4
    dat.day1.memorytest.delayed.config = config;
    fprintf('\n**************\nDELAYED RECALL \n**************\n\n\n')
end

%% run

% configs
config_display(0, 4, [0 0 0], [1 1 1], 'Arial', 20, 10, 0); % 1024x768
% config_display(screenMode, screenRes, bg ,fontcol, fontName, fontSize, number_of_buffers,0);   % open graphics window
config_sound;
config_keyboard;

monitor = get(0,'screensize');
bg_w = monitor(3); bg_h = monitor(4);

% scanner input
if MRIFlag == 1
    config_serial(dat.day1.memorytest.immediate.MRIinfo.scanPort); % links to scanner for waitslice
end

%%%%%%%%%%%%
start_cogent
%%%%%%%%%%%%

if MRIFlag == 1
    cgopen(4,0,0,2)
else
    cgopen(4,0,0,1)
end

map     = getkeymap; % define keyboard IDs
OldKey  = 28;
NewKey  = 29;
Surekey = 28;
Guesskey= 29;
buff    = 2;

Picture = 2;% name buffers
Background = 3;
Questions = 4;
Answers = 5;
Instruction = 6;
Breaks = 7;
Perfgr = 8;

% print instructions
cd(paths.fb);
cgloadbmp(1,'memorytest_instruction1.bmp');
cgsetsprite(0);
cgdrawsprite(1,0,0);
cgflip(0,0,0); waitkeydown(inf,[OldKey,NewKey]); clearkeys

if eyeflag==1
    Eyelinkerror = Eyelink('StartRecording'); % start recording to the file
    if Eyelinkerror ~= 0
        return
        error('Eyetracker failed')
    end
    Eyelink('Message',num2str([0 100]))
end

% let the scanner start the task, allow n dummy volumes
if MRIFlag == 1
    
    % hail the operator
    cd(paths.fb)
    cgloadbmp(1,'operator_hail.bmp');
    cgsetsprite(0);
    cgdrawsprite(1,0,0);
    cgflip(0,0,0);
    waitkeydown(inf,map.Space)
    
    % start the dummy scan
    cd(paths.fb)
    cgloadbmp(1,'ready.bmp');
    cgsetsprite(0);
    cgdrawsprite(1,0,0);
    t0_standby=cgflip(0,0,0);
    scannerinput = 33;
    
    %%%%%%%%%%%%%%
    tic % ding ding ding
    %%%%%%%%%%%%%%
    
    [k, dat.day1.memorytest.immediate.MRIinfo.dummyscanfirst] = waitkeydown(inf,scannerinput); trig_1st = toc; % the first dummy scan
    for dum = 1:(5 - 2)
        waitkeydown(inf,scannerinput);
    end
    [k, dat.day1.memorytest.immediate.MRIinfo.dummyscanlast] = waitkeydown(inf,scannerinput); trig_5th = toc; % the last dummy scan
    
    cd(paths.fb); cgloadbmp(1,'checker_b.bmp',1152,864);
    cgsetsprite(0);
    cgdrawsprite(1,0,0);
    cgfont('Arial',300);
    cgpencol(.5,.5,.5)
    cgtext('+',0,0);  % Prepare buffer 1 (fixation cross)
    t0_fix0 = cgflip(0,0,0);  t0_fix0_raw = toc;    % load fixation cross
    if eyeflag==1
        Eyelink('Message',[num2str(0) 'f0'])  % send trialnumber and stimulustype to eyetracker
    end
    wait(dat.day1.memorytest.immediate.MRIinfo.TR*1000)
    
else
    
    %%%%%%%%%%%%%%
    tic % ding ding ding
    %%%%%%%%%%%%%%
    
    cd(paths.fb); cgloadbmp(1,'checker_b.bmp',1152,864);
    cgsetsprite(0);
    cgdrawsprite(1,0,0);
    cgfont('Arial',300)
    cgpencol(.5,.5,.5)
    cgtext('+',0,0);
    cgflip(0,0,0);
    if eyeflag==1
        Eyelink('Message',[num2str(0) 'f0'])  % send trialnumber and stimulustype to eyetracker
    end
    wait(1000)
    
end


clear trl c1 c2
c1 = 0; c2 = 0; results = []; results.SOT = [];
if MRIFlag == 1
results.SOT.cgflip.t0_fix0 = t0_fix0;
results.SOT.cgflip.t0_standby = t0_standby;
results.SOT.raw.t0_fix0 = t0_fix0_raw;
results.SOT.raw.trig_1st = trig_1st;
results.SOT.raw.trig_5th = trig_5th;
end
resp = []; accu = []; confi = [];
f1c = 0; sc = 0;
SOT_f1=[]; SOT_s=[]; SOT_cg_s=[]; SOT_cg_f1=[];
stimcount = 0;
for trl = 1:length(stimvec)
    
    if stimvec(trl) ~= 777
        
        stimcount = stimcount+1;
        fprintf('\ntrial %d\n',stimcount)
                
        %%%%%% load background %%%%%%
        cd(paths.fb); 
        cgloadbmp(1,'checker_b.bmp',1152,864);
        
        %%%%%% load stimulus %%%%%%
        clearkeys;
        cd(paths.stim);
        cgloadbmp(2,stimlist{stimcount,1},500,360); % prepare buffer 2 (image)
        
        %%%%%% draw stimulus screen %%%%%%
        sc=sc+1; SOT_s(sc,1)=toc;
        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgdrawsprite(2,0,0);
        SOT_cg_s(trl,1) = cgflip(0,0,0);
        if eyeflag==1
            Eyelink('Message',[num2str(stimcount) 's'])  % send trialnumber and stimulustype to eyetracker
        end
        wait(config.timing.cue);
        
        %%%%%% load new-old screen %%%%%%
        cd(paths.fb);
        cgloadbmp(1,'checker_b.bmp',1152,864);
        cgloadbmp(2,'newold.bmp');
        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgdrawsprite(2,0,0); tmp = toc;
        
        newoldonset = cgflip(0,0,0); % newold screen
        if eyeflag==1
            Eyelink('Message',[num2str(stimcount) 'r'])  % send trialnumber and stimulustype to eyetracker
        end
        waitkeydown(config.timing.response,[OldKey,NewKey]);
        
        [key_M, t_M, n_M] = getkeydown;
        if n_M == 0 %if they didn't press anything
            wait(0);
        else % otherwise
            %%%%%% load new-old answer screen %%%%%%
            cgloadbmp(1,'checker_b.bmp',1152,864);
            cgloadbmp(2,'newold_received.bmp');
            cgsetsprite(0);
            cgdrawsprite(1,0,0); cgdrawsprite(2,0,0); cgflip(0,0,0);
            rt = t_M(1)/1000 - newoldonset; rt = rt*1000; % and their reaction time is the time they pressed the button-the time the stimulus apprered
            wait(2000-rt);
        end
        
        
        
        
        
        
        %%%%%% load confidence screen %%%%%%
        cd(paths.fb);
        cgloadbmp(1,'checker_b.bmp',1152,864);
        cgloadbmp(2,'memorytest_confidence.bmp');
        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgdrawsprite(2,0,0); tmp = toc;
        
        confonset = cgflip(0,0,0); % conf screen
        if eyeflag==1
            Eyelink('Message',[num2str(stimcount) 'c'])  % send trialnumber and stimulustype to eyetracker
        end
        waitkeydown(config.timing.response,[Surekey,Guesskey])
        
        [key_SG, t_SG, n_SG] = getkeydown;
        if n_M == 0 %if they didn't press anything
            wait(0);
        else % otherwise
            %%%%%% load confidence answer screen %%%%%%
            cgloadbmp(1,'checker_b.bmp',1152,864);
            cgloadbmp(2,'memorytest_confidence_received.bmp');
            cgsetsprite(0);
            cgdrawsprite(1,0,0); cgdrawsprite(2,0,0); cgflip(0,0,0);
            try
                rt_cnf = t_SG(1)/1000 - confonset; rt_cnf = rt_cnf*1000; % and their reaction time is the time they pressed the button-the time the stimulus apprered
                wait(2000-rt_cnf);
            catch
                wait(1500)
            end
        end
        
        % record confidence
        if key_SG == Surekey
            confi = 1;
        elseif key_SG == Guesskey
            confi = 0;
        end
        
        % accuracy
        if key_M == OldKey
            resp = 1;
            if stimlist{stimcount,4} == 1
                fprintf('\ncorrect\n')
                accu = 1;
            else
                fprintf('\nincorrect\n')
                accu = 0;
            end
        elseif key_M == NewKey
            resp = 2;
            if stimlist{stimcount,4} == 2
                fprintf('\ncorrect\n')
                accu = 1;
            else
                fprintf('\nincorrect\n')
                accu = 0;
            end
        end
        
        % RT
        if n_M == 0 %if they didn't press anything
            rt = NaN; % mark their reaction time as nan
        else % otherwise
            rt = t_M(1) - newoldonset; % and their reaction time is the time they pressed the button-the time the stimulus apprered
        end
        
        
        %%%%%% load and present fix x %%%%%%
        cd(paths.fb); cgloadbmp(1,'checker_b.bmp',1152,864);
        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgfont('Arial',300)
        cgpencol(.5,.5,.5)
        cgtext('+',0,0);
        f1c=f1c+1; SOT_f1(f1c,1)=toc;
        SOT_cg_f1(trl,1) = cgflip(0,0,0); % present fixation cross
        if eyeflag==1
            Eyelink('Message',[num2str(stimcount) 'f'])  % send trialnumber and stimulustype to eyetracker
        end
        wait(config.timing.fixation1(stimcount));
        
        
        % record
        if n_M == 0 || n_SG == 0
%             results.keypress(stimcount,1:2)   = [NaN NaN];
            results.oldnew_resp(stimcount,1)  = NaN;
            results.accuracy(stimcount,1)     = NaN;
            results.confidence(stimcount,1)   = NaN;
            results.rt(stimcount,1)           = NaN;
        else
%             try
%             results.keypress(stimcount,1:2)   = [key_M key_SG];
%             catch
%                 results.keypress(stimcount,1:2)   = [NaN NaN];
%             end
            results.oldnew_resp(stimcount,1)  = resp;
            results.accuracy(stimcount,1)     = accu;
            results.confidence(stimcount,1)   = confi;
            results.rt(stimcount,1)           = rt;
        end
        %     results.oldnew_resp         = resp;
        %     results.accuracy(trl,1)     = accu;
        %     results.confidence(trl,1)   = confi;
        %     results.rt(trl,1)           = rt;
        results.trl(stimcount,:)          = stimlist(stimcount,:);
        results.all = horzcat(results.oldnew_resp,results.accuracy,results.confidence,results.rt);
        
        % calculate cumulative SOT
        results.SOT.cumulative.fix  = cumsum(SOT_f1);
        results.SOT.cumulative.stim = cumsum(SOT_s);
        results.SOT.raw.fix  = SOT_f1;
        results.SOT.raw.stim = SOT_s;
        results.SOT.cgflip.stim = SOT_cg_s;
        results.SOT.cgflip.fix  = SOT_cg_f1;
        
        % intermediate save
        if DayNumber == 1
            dat.day1.memorytest.labels = {'stimlist_all:' 'filenames' 'Nature(1)/Urban(2)/Private(3)/Public(4)' 'Memorability' 'Old(1)/New(2)'; ...
                'results:' 'response(old-1/new-2)' 'accuracy' 'confidence rating' 'rt'};
            if TaskType == 3 % immediate
                dat.day1.memorytest.immediate.results.SOT    = results.SOT;
                dat.day1.memorytest.immediate.results.resp   = results.oldnew_resp;
                dat.day1.memorytest.immediate.results.accu   = results.accuracy;
                dat.day1.memorytest.immediate.results.confi  = results.confidence;
                dat.day1.memorytest.immediate.results.rt     = results.rt;
                dat.day1.memorytest.immediate.results.all    = results.all;
                dat.day1.memorytest.immediate.results.trl    = results.trl;
                dat.day1.memorytest.immediate.ETinfo         = eyetracker_filename;
                dat.day1.memorytest.immediate.config.keymap  = map;
            elseif TaskType == 4 % delayed
                dat.day1.memorytest.delayed.results.SOT      = results.SOT;
                dat.day1.memorytest.delayed.results.resp     = results.oldnew_resp;
                dat.day1.memorytest.delayed.results.accu     = results.accuracy;
                dat.day1.memorytest.delayed.results.confi    = results.confidence;
                dat.day1.memorytest.delayed.results.rt       = results.rt;
                dat.day1.memorytest.delayed.results.all      = results.all;
                dat.day1.memorytest.delayed.results.trl      = results.trl;
                dat.day1.memorytest.delayed.ETinfo           = eyetracker_filename;
                dat.day1.memorytest.delayed.config.keymap    = map;
            end
            save([paths.behav num2str(ID) '_' num2str(DayNumber) '_backup.mat'],'dat')
        elseif DayNumber == 2
            dat.day2.memorytest.labels = {'stimlist_all:' 'filenames' 'Nature(1)/Urban(2)/Private(3)/Public(4)' 'Memorability' 'Old(1)/New(2)'; ...
                'results:' 'response(old-1/new-2)' 'accuracy' 'confidence rating' 'rt'};
            if TaskType == 3 % immediate
                dat.day2.memorytest.immediate.results.SOT    = results.SOT;
                dat.day2.memorytest.immediate.results.resp   = results.oldnew_resp;
                dat.day2.memorytest.immediate.results.accu   = results.accuracy;
                dat.day2.memorytest.immediate.results.confi  = results.confidence;
                dat.day2.memorytest.immediate.results.rt     = results.rt;
                dat.day2.memorytest.immediate.results.all    = results.all;
                dat.day2.memorytest.immediate.results.trl    = results.trl;
                dat.day2.memorytest.immediate.ETinfo         = eyetracker_filename;
                dat.day2.memorytest.immediate.config.keymap  = map;
            elseif TaskType == 4 % delayed
                dat.day2.memorytest.delayed.results.SOT      = results.SOT;
                dat.day2.memorytest.delayed.results.resp     = results.oldnew_resp;
                dat.day2.memorytest.delayed.results.accu     = results.accuracy;
                dat.day2.memorytest.delayed.results.confi    = results.confidence;
                dat.day2.memorytest.delayed.results.rt       = results.rt;
                dat.day2.memorytest.delayed.results.all      = results.all;
                dat.day2.memorytest.delayed.results.trl      = results.trl;
                dat.day2.memorytest.delayed.ETinfo           = eyetracker_filename;
                dat.day2.memorytest.delayed.config.keymap    = map;
            end
            save([paths.behav num2str(ID) '_' num2str(DayNumber) '_backup.mat'],'dat')
        end
        
    elseif stimvec(trl) == 777
        
        cd(paths.fb);
        cgloadbmp(1, 'instruction_break.bmp');      
        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgflip(0,0,0);
        waitkeydown(inf,[OldKey,NewKey]);
        
    end
    
end % close the trial loop


%% wrap up

% calculate final cumulative SOT
results.SOT.cumulative.fix  = cumsum(SOT_f1);
results.SOT.cumulative.stim = cumsum(SOT_s);
results.SOT.raw.fix  = SOT_f1;
results.SOT.raw.stim = SOT_s;
results.SOT.cgflip.stim = SOT_cg_s;
results.SOT.cgflip.fix  = SOT_cg_f1;

% save data
if DayNumber == 1
    dat.day1.memorytest.labels = {'stimlist_all:' 'filenames' 'Nature(1)/Urban(2)/Private(3)/Public(4)' 'Memorability' 'Old(1)/New(2)'; ...
        'results:' 'response(old-1/new-2)' 'accuracy' 'confidence rating' 'rt'};
    if TaskType == 3 % immediate
        dat.day1.memorytest.immediate.ETinfo         = eyetracker_filename;
        dat.day1.memorytest.immediate.config.keymap  = map;
        dat.day1.memorytest.immediate.SOT            = results.SOT;
        dat.day1.memorytest.immediate.results.resp   = results.oldnew_resp;
        dat.day1.memorytest.immediate.results.accu   = results.accuracy;
        dat.day1.memorytest.immediate.results.confi  = results.confidence;
        dat.day1.memorytest.immediate.results.rt     = results.rt;
        dat.day1.memorytest.immediate.results.all    = results.all;
        dat.day1.memorytest.immediate.results.trl    = results.trl;
        
    elseif TaskType == 4 % delayed
        dat.day1.memorytest.delayed.ETinfo           = eyetracker_filename;
        dat.day1.memorytest.delayed.config.keymap    = map;
        dat.day1.memorytest.delayed.results.SOT      = results.SOT;
        dat.day1.memorytest.delayed.results.resp     = results.oldnew_resp;
        dat.day1.memorytest.delayed.results.accu     = results.accuracy;
        dat.day1.memorytest.delayed.results.confi    = results.confidence;
        dat.day1.memorytest.delayed.results.rt       = results.rt;
        dat.day1.memorytest.delayed.results.all      = results.all;
        dat.day1.memorytest.delayed.results.trl      = results.trl;
        
    end
    save([paths.behav num2str(ID) '_' num2str(DayNumber) '.mat'],'dat')
elseif DayNumber == 2
    
    dat.day2.memorytest.labels = {'stimlist_all:' 'filenames' 'Nature(1)/Urban(2)/Private(3)/Public(4)' 'Memorability' 'Old(1)/New(2)'; ...
        'results:' 'response(old-1/new-2)' 'accuracy' 'confidence rating' 'rt'};
    if TaskType == 3 % immediate
        dat.day2.memorytest.immediate.ETinfo         = eyetracker_filename;
        dat.day2.memorytest.immediate.config.keymap  = map;
        dat.day2.memorytest.immediate.SOT            = results.SOT;
        dat.day2.memorytest.immediate.results.resp   = results.oldnew_resp;
        dat.day2.memorytest.immediate.results.accu   = results.accuracy;
        dat.day2.memorytest.immediate.results.confi  = results.confidence;
        dat.day2.memorytest.immediate.results.rt     = results.rt;
        dat.day2.memorytest.immediate.results.all    = results.all;
        dat.day2.memorytest.immediate.results.trl    = results.trl;
        
    elseif TaskType == 4 % delayed
        dat.day2.memorytest.delayed.ETinfo           = eyetracker_filename;
        dat.day2.memorytest.delayed.config.keymap    = map;
        dat.day2.memorytest.delayed.results.SOT      = results.SOT;
        dat.day2.memorytest.delayed.results.resp     = results.oldnew_resp;
        dat.day2.memorytest.delayed.results.accu     = results.accuracy;
        dat.day2.memorytest.delayed.results.confi    = results.confidence;
        dat.day2.memorytest.delayed.results.rt       = results.rt;
        dat.day2.memorytest.delayed.results.all      = results.all;
        dat.day2.memorytest.delayed.results.trl      = results.trl;
        
    end
    save([paths.behav num2str(ID) '_' num2str(DayNumber) '.mat'],'dat')
end

% terminate protocol
cd(paths.fb)
cgloadbmp(1,'instruction_ende.bmp');
cgsetsprite(0);
cgdrawsprite(1,0,0);
cgflip(0,0,0); wait(3000);
cgshut;

fprintf('\n%%%%%%%%%%%%%%%%%%%%%%\n COMPLETED  \n%%%%%%%%%%%%%%%%%%%%%%\n')
fprintf('\nMean Accuracy: %3.3f\n', nanmean(results.accuracy))

end
