 %% MAIN TASK : Scene classification task
%  In emergency contact Alex:
%              +49 162 713 20 49
%              yeo-jin.yi.15@ucl.ac.uk
% See also StarterScript_MRPET, runPracticeTask_MRPET,
%       runMemoryTest_MRPET

%%  work log

%   29_10_2018 created the script

%% START

function [dat] = runMainTask_MRPET(ID,EyeFlag,DayNumber,MRIFlag,CnterBal,eyetracker_filename,PATH_echo,Schedules)

global eyeflag eyetracker_filename paths

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
%         EyelinkDoTrackerSetup(el); % changed this bit to see what happens to the data
        
        
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

% set paths
paths.parent   = PATH_echo;
if MRIFlag == 1
paths.stim     = [ paths.parent 'MRPET\stim_mirrored\' ];  %%%%% changed
else
paths.stim     = [ paths.parent 'MRPET\stim_original\' ];  %%%%% changed
end
paths.fb       = [ paths.parent 'MRPET\feedback\' ];
paths.pupil    = [ paths.parent 'MRPET\raw\pupil\' ];
paths.behav    = [ paths.parent 'MRPET\raw\behav\' ];
% paths.cogent   = 'C:\Studies\Grid_experiment2\Cogent2000v1.32';
% paths.cogent   = '\\fs-md\users\yiy\Documents\MATLAB\Cogent2000v1.33';
paths.cogent   = '\\fs-md\users\yiy\Documents\MATLAB\Cogent2000v1.33';
load([ paths.stim 'Memorability_API_midrange_20190611.mat' ]) %%%%% changed
load([ paths.behav num2str(ID) '_' num2str(DayNumber) '.mat' ])
addpath(genpath(paths.cogent))

% create data structure
if DayNumber == 1
    dat.day1.maintask  = [];
    dat.day1.maintask.ETinfo  = eyetracker_filename;
    dat.day1.maintask.StimSchedule = Schedules;
    if Schedules == 1
        load([ paths.stim 'Main_GA_designstruct_D1_ver1.mat' ])
    elseif Schedules == 2
        load([ paths.stim 'Main_GA_designstruct_D1_ver2.mat' ])
    end
elseif DayNumber == 2
    dat.day2.maintask  = [];
    dat.day2.maintask.ETinfo  = eyetracker_filename;
    dat.day2.maintask.StimSchedule = Schedules;
    if Schedules == 1
        load([ paths.stim 'Main_GA_designstruct_D2_ver1.mat' ])
    elseif Schedules == 2
        load([ paths.stim 'Main_GA_designstruct_D2_ver2.mat' ])
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
%---------- for MR-PET pilot -----------%
if Schedules == 1
    
    % INDOORS  dataset : Nature(1) / Urban(2)
    % OUTDOORS dataset : Private(3) / Public(4)
    
    if CnterBal == 1 % if land/private = reward
        if DayNumber == 2 % small rew
            rewcond = 1; puncond = 2;
            remind_scr = 'nature_small';
        elseif DayNumber == 1 % big rew
            rewcond = 3; puncond = 4;
            remind_scr = 'private_big';
        end
        
    elseif CnterBal == 2 % if city/public == reward
        if DayNumber == 2 % small rew
            rewcond = 2; puncond = 1;
            remind_scr = 'urban_small';
        elseif DayNumber == 1 % big rew
            rewcond = 4; puncond = 3;
            remind_scr = 'public_big';
        end
    end
    
elseif Schedules == 2 % OUTDOORS dataset is Day 1
    
    % INDOORS  dataset : Nature(1) / Urban(2)
    % OUTDOORS dataset : Private(3) / Public(4)
    
    if CnterBal == 1 % if land/private = reward
        if DayNumber == 1 % big rew
            rewcond = 1; puncond = 2;
            remind_scr = 'nature_big';
        elseif DayNumber == 2 % small rew
            rewcond = 3; puncond = 4;
            remind_scr = 'private_small';
        end
        
    elseif CnterBal == 2 % if city/public == reward
        if DayNumber == 1 % big rew
            rewcond = 2; puncond = 1;
            remind_scr = 'urban_big';
        elseif DayNumber == 2 % small rew
            rewcond = 4; puncond = 3;
            remind_scr = 'public_small';
        end
    end
end

%----------------------------------------%

%------------- for 3T pilot -------------%
% if CnterBal == 1 % if indoor = reward
%     rewcond = 1; puncond = 2;
%     cat1_1 = 3; cat1_2 = 4; cat2_1 = 1; cat2_2 = 2; %%%%% changed & added
% elseif CnterBal == 2 % if outdoor == reward
%     rewcond = 2; puncond = 1;
%     cat1_1 = 1; cat1_2 = 2; cat2_1 = 3; cat2_2 = 4; %%%%% changed & added
% end
%----------------------------------------%

% timings
config.timing.outcome  = 2500;
% fixjit1 = repmat(500:50:700,1,48); fixjit2 = repmat(500:50:700,1,48);
% config.timing.fixation1= fixjit1(randperm(length(fixjit1)));
% config.timing.fixation2= fixjit2(randperm(length(fixjit2)));
tmp = [(design_struct.eventlist(:,4)+design_struct.eventlist(:,5)-3.5)']*1000;
config.timing.fixation1= horzcat(tmp(1:45), [NaN], tmp(46:90), [NaN], ...
    tmp(91:135), [NaN], tmp(136:180)); clear tmp
tmp = repmat(1000:500:2000,1,60);
config.timing.dot      = tmp(randperm(length(tmp))); clear tmp
config.timing.cue      = 1500;
config.timing.stim     = 2500;
config.timing.response = 2000;

% --- for breezeing through --- %
% config.timing.fixation1= zeros(1,240);
% config.timing.fixation2= zeros(1,240);
% config.timing.stim     = 1;
% config.timing.cue      = 1;
% config.timing.response = 1;
% config.timing.dot      = zeros(1,300);
% ----------------------------- %

% trials
config.numtrials.total = 180;
if DayNumber == 1 % the first day, more rewards
    config.numtrials.rew   = 80;
    config.numtrials.pun   = 80;
    config.numtrials.null  = 20;
elseif DayNumber == 2 % the second day, fewer rewards
    config.numtrials.rew   = 80;
    config.numtrials.pun   = 80;
    config.numtrials.null  = 20;
end

% contingency vector
trlvec  = 1:config.numtrials.total;
rewmat  = design_struct.eventlist(:,3)';
rewmat_rewarded = find(rewmat==1); 
contingency_rew = horzcat(ones(1,length(rewmat_rewarded))); % rewarded: 1
contingency_rew = contingency_rew(randperm(length(contingency_rew)));
rewmat_punished = find(rewmat==2);
contingency_pun = horzcat(ones(1,length(rewmat_punished))*2); % neutral trials all 2
contingency_pun = contingency_pun(randperm(length(contingency_pun)));
rewmat_nulls    = find(rewmat==0);
contingency_null= horzcat(ones(1,length(rewmat_nulls))*0);
contingency_null= contingency_null(randperm(length(contingency_null)));
config.stim.rewmat_rewarded = vertcat(rewmat_rewarded,contingency_rew);
config.stim.rewmat_punished = vertcat(rewmat_punished,contingency_pun);
config.stim.rewmat_nulls    = vertcat(rewmat_nulls,contingency_null);

% roster
c1 = 0; c3 = 0; stim_pun_1 = []; stim_rew_1 = [];
for i = 1:size(memorability_API_midrange_20190611,1) %%%%% changed
    if cell2mat(memorability_API_midrange_20190611(i,3)) == puncond %%%%% changed
        c1 = c1+1;
        stim_pun_1{c1,1} = memorability_API_midrange_20190611{i,1}; %%%%% changed
        stim_pun_1{c1,2} = memorability_API_midrange_20190611{i,4}; %%%%% changed
        stim_pun_1{c1,3} = memorability_API_midrange_20190611{i,2}; %%%%% changed
        stim_pun_1{c1,4} = memorability_API_midrange_20190611{i,3}; %%%%% changed & added
    elseif cell2mat(memorability_API_midrange_20190611(i,3)) == rewcond %%%%% changed
        c3 = c3+1;
        stim_rew_1{c3,1} = memorability_API_midrange_20190611{i,1}; %%%%% changed
        stim_rew_1{c3,2} = memorability_API_midrange_20190611{i,4}; %%%%% changed
        stim_rew_1{c3,3} = memorability_API_midrange_20190611{i,2}; %%%%% changed
        stim_rew_1{c3,4} = memorability_API_midrange_20190611{i,3}; %%%%% changed & added
    end % close the contingency conditional %%%%% changed
end % close the roster loop

config.stim.fname      = cellstr(vertcat(stim_rew_1{1:config.numtrials.rew,1}, stim_pun_1{1:config.numtrials.pun,1}));
config.stim.fname(:,2) = num2cell(vertcat(stim_rew_1{1:config.numtrials.rew,2}, stim_pun_1{1:config.numtrials.pun,2}));
config.stim.fname(:,3) = num2cell(vertcat(stim_rew_1{1:config.numtrials.rew,3}, stim_pun_1{1:config.numtrials.pun,3}));
config.stim.fname(:,4) = num2cell(vertcat(stim_rew_1{1:config.numtrials.rew,4}, stim_pun_1{1:config.numtrials.pun,4}));

% schedule presentation
tmpvec       = trlvec;
stimvec_temp = [];
rhc = 0; phc = 0; nhc = 0; rwc = 0; pnc = 0; nwc = 0;
for sc = 1:length(rewmat)
    if rewmat(sc) == 1
        rhc = rhc+1; rwc = rwc+1;
        stimvec_temp(1,sc) = tmpvec(sc); % trial order
        stimvec_temp(2,sc) = rewcond; % reward contingency
        stimvec_temp(3,sc) = stim_rew_1{rhc,2}; % outdoor/indoor
        stimvec_temp(4,sc) = contingency_rew(rwc); % reward/punishment actual contingency
        stimvec_flist{1,sc}= stim_rew_1{rhc,1};
        stimvec_flist{2,sc}= stim_rew_1{rhc,2};
        stimvec_flist{3,sc}= stim_rew_1{rhc,3};
        stimvec_flist{4,sc}= stim_rew_1{rhc,4};
    elseif rewmat(sc) == 2
        phc = phc+1; pnc = pnc+1;
        stimvec_temp(1,sc) = tmpvec(sc);
        stimvec_temp(2,sc) = puncond;
        stimvec_temp(3,sc) = stim_pun_1{phc,2};
        stimvec_temp(4,sc) = contingency_pun(pnc);
        stimvec_flist{1,sc}= stim_pun_1{phc,1};
        stimvec_flist{2,sc}= stim_pun_1{phc,2};
        stimvec_flist{3,sc}= stim_pun_1{phc,3};
        stimvec_flist{4,sc}= stim_pun_1{phc,4};
    elseif rewmat(sc) == 0
        nhc = nhc+1; nwc = nwc+1;
        stimvec_temp(1,sc) = tmpvec(sc);
        stimvec_temp(2,sc) = 0;
        stimvec_temp(3,sc) = 0;
        stimvec_temp(4,sc) = contingency_null(nwc);
        stimvec_flist{1,sc}= 'Null';
        stimvec_flist{2,sc}= 0;
        stimvec_flist{3,sc}= 'Null';
        stimvec_flist{4,sc}= 0;
        
    end
end


% for sc = 1:length(tmpvec)
%     try
%         if config.stim.fname{tmpvec(sc),2} == rewcond % if rewarded category
%             rc = rc+1;
%             stimvec_temp(1,sc) = tmpvec(sc);
%             stimvec_temp(2,sc) = config.stim.rewmat_rewarded(rc);
%             stimvec_temp(3,sc) = config.stim.fname{tmpvec(sc),2};
%         elseif config.stim.fname{tmpvec(sc),2} == puncond % if punished category
%             nc = nc+1;
%             stimvec_temp(1,sc) = tmpvec(sc);
%             stimvec_temp(2,sc) = config.stim.rewmat_punished(nc);
%             stimvec_temp(3,sc) = config.stim.fname{tmpvec(sc),2};
%         end
%     catch
%         stimvec_temp(1,sc) = tmpvec(sc);
%         stimvec_temp(2,sc) = 0;
%         stimvec_temp(3,sc) = 0;
%     end
% end

stimvec_mod     = horzcat(stimvec_temp(:,1:45), [777;-1;-1;-1], stimvec_temp(:,46:90), [777;-1;-1;-1], ...
    stimvec_temp(:,91:135), [777;-1;-1;-1], stimvec_temp(:,136:180));
stimvec_compl   = stimvec_mod;
config.stim.stimvec = stimvec_mod;
config.stim.fname   = stimvec_flist';

if DayNumber == 1
    dat.day1.maintask.config = config;
elseif DayNumber == 2
    dat.day2.maintask.config = config;
end

fprintf('\n\nTASK SETUP DONE')
fprintf('\nINITIALISING MAIN TASK\n\n')
if DayNumber == 1
    fprintf('\n**** DAY 1 ****\n')
else
    fprintf('\n**** DAY 2 ****\n')
end

%% run

% configs
config_display(1, 4, [0 0 0], [1 1 1], 'Arial', 20, 10, 0); % 1024x768
% config_display(screenMode, screenRes, bg ,fontcol, fontName, fontSize, number_of_buffers,0);   % open graphics window
% config_sound;
config_log( [ paths.behav 'Cogent_' num2str(datevec(now),'-%02.0f') '_ID' num2str(ID) '_maintask.log'] );
config_keyboard;

% scanner input
if MRIFlag == 1
    if DayNumber == 1
        config_serial(dat.day1.maintask.MRIinfo.scanPort); % links to scanner for waitslice
    elseif DayNumber == 2
        config_serial(dat.day2.maintask.MRIinfo.scanPort); % links to scanner for waitslice
    end
end

monitor = get(0,'screensize');
bg_w = monitor(3); bg_h = monitor(4);

%%%%%%%%%%%%
start_cogent
%%%%%%%%%%%%

if MRIFlag == 1
    cgloadlib
    cgopen(4,0,0,2)
else
    cgloadlib
%     cgopen(5,0,0,0) % at alex's laptopcg
    cgopen(4,0,0,1)
end

map     = getkeymap; % define keyboard IDs
if CnterBal == 1  % if land/private reward / inside
    RewKey  = 28; % yellow / inside
    PunKey  = 29; % blue / outside
elseif CnterBal == 2 % if urban/private / outside
    RewKey  = 29;
    PunKey  = 28;
end

Picture = 2;% name buffers
Background = 3;
Questions = 4;
Answers = 5;
Instruction = 6;
Breaks = 7;
Perfgr = 8;

% print instructions
cd(paths.fb)

% ----------------- keytest ----------------- %

% left
tic;
if MRIFlag == 1
cgloadbmp(1,'keytest_left.bmp');
else
cgloadbmp(1,'keytest_left_orig.bmp');
end
cgsetsprite(0);
cgdrawsprite(1,0,40); tmp = toc;
tmpkeyonset = cgflip(0,0,0);
waitkeydown(5000,28);

[tmpkey, tmpt, tmpn] = getkeydown;
if tmpn == 0
else
if MRIFlag == 1
cgloadbmp(1,'keytest_left_done.bmp');
else
cgloadbmp(1,'keytest_left_done_orig.bmp');
end
cgsetsprite(0);
cgdrawsprite(1,0,40); cgflip(0,0,0);
tempkeyrt = tmpt(1)/1000 - tmpkeyonset; tempkeyrt = tempkeyrt*1000;
wait(5000-tempkeyrt);
end


% right
tic;
if MRIFlag == 1
cgloadbmp(1,'keytest_right.bmp');
else
cgloadbmp(1,'keytest_right_orig.bmp');
end
cgsetsprite(0);
cgdrawsprite(1,0,40); tmp = toc;
tmpkeyonset = cgflip(0,0,0);
waitkeydown(5000,29);

[tmpkey, tmpt, tmpn] = getkeydown;
if tmpn == 0
else
if MRIFlag == 1
cgloadbmp(1,'keytest_right_done.bmp');
else
cgloadbmp(1,'keytest_right_done_orig.bmp');
end
cgsetsprite(0);
cgdrawsprite(1,0,40); cgflip(0,0,0);
tempkeyrt = tmpt(1)/1000 - tmpkeyonset; tempkeyrt = tempkeyrt*1000;
wait(5000-tempkeyrt);
end

% ------------------------------------------- %


% --------------- instruction --------------- %

if MRIFlag == 1
%     cgloadbmp(1,'emomem_inst_task.bmp');
    if rewcond == 1 || rewcond == 2
        cgloadbmp(1,'emomem_inst_task_natureurban.bmp');        
    elseif rewcond == 3 || rewcond == 4
        cgloadbmp(1,'emomem_inst_task_privatepublic.bmp');        
    end
else
%     cgloadbmp(1,'emomem_inst_task_orig.bmp');
    if rewcond == 1 || rewcond == 2
        cgloadbmp(1,'emomem_inst_task_natureurban_orig.bmp');        
    elseif rewcond == 3 || rewcond == 4
        cgloadbmp(1,'emomem_inst_task_privatepublic_orig.bmp');        
    end
end
cgsetsprite(0);
cgdrawsprite(1,0,50);
cgflip(0,0,0); waitkeydown(inf,[RewKey,PunKey]); % VP input
clearkeys

% ------------------------------------------- %


% -------------- rew reminders -------------- %

if MRIFlag == 1
    cgloadbmp(1,[remind_scr '.bmp']);
else
    cgloadbmp(1,[remind_scr '_orig.bmp']);
end
cgsetsprite(0);
cgdrawsprite(1,0,50);
cgflip(0,0,0); waitkeydown(inf,[RewKey,PunKey]); % VP input
clearkeys


if eyeflag==1
    Eyelinkerror = Eyelink('StartRecording'); % start recording to the file
    if Eyelinkerror ~= 0
        return
        error ('Eyetracker failed')
    end
    Eyelink('Message',num2str([0 100]))
end


% let the scanner start the task, allow n dummy volumes
if MRIFlag == 1
    
    scannerinput = 33;
    
    % hail the operator
    fprintf(['\n *************** \n please press SPACE to continue!! \n *************** \n '])
    cd(paths.fb)
    cgloadbmp(1,'standby.bmp');
    cgsetsprite(0);
    cgdrawsprite(1,0,0);
    cgflip(0,0,0);
    waitkeydown(inf,[map.Space]) 
    
    % start the dummy scan
    cd(paths.fb)
    cgloadbmp(1,'standby_done.bmp');
    cgsetsprite(0);
    cgdrawsprite(1,0,0);
    t0_standby=cgflip(0,0,0);    
    
    %%%%%%%%%%%%%%
    tic % ding ding ding
    %%%%%%%%%%%%%%
    
    if DayNumber == 1
        [k, dat.day1.maintask.MRIinfo.dummyscanfirst] = waitkeydown(inf,scannerinput); trig_1st = toc; % the first dummy scan
        fprintf('\n *** 1 *** \n')
        for dum = 1:(dat.day1.maintask.MRIinfo.dummy - 2)
            waitkeydown(inf,scannerinput);
            fprintf('\n *** %1.0f *** \n',dum+1)
        end
        [k, dat.day1.maintask.MRIinfo.dummyscanlast] = waitkeydown(inf,scannerinput); trig_5th = toc; % the last dummy scan
        fprintf('\n *** 5 *** \n')
    else
        [k, dat.day2.maintask.MRIinfo.dummyscanfirst] = waitkeydown(inf,scannerinput); trig_1st = toc; % the first dummy scan
        fprintf('\n *** 1 *** \n')
        for dum = 1:(dat.day2.maintask.MRIinfo.dummy - 2)
            waitkeydown(inf,scannerinput);
            fprintf('\n *** %1.0f *** \n',dum+1)
        end
        [k, dat.day2.maintask.MRIinfo.dummyscanlast] = waitkeydown(inf,scannerinput); trig_5th = toc; % the last dummy scan
        fprintf('\n *** 5 *** \n')
    end
    
    cd(paths.fb); cgloadbmp(1,'checker_b.bmp',1152,864);
    cgsetsprite(0);
    cgdrawsprite(1,0,0);
    cgfont('Arial',300);
    cgpencol(.5,.5,.5)
    cgtext('+',0,0);  % Prepare buffer 1 (fixation cross)
    t0_fix0 = cgflip(0,0,0); t0_fix0_raw = toc;       % load fixation cross
    if eyeflag==1
        Eyelink('Message',[num2str(0) 'f0'])  % send trialnumber and stimulustype to eyetracker
    end
    if DayNumber == 1
        wait(dat.day1.maintask.MRIinfo.TR*1000)
    else
        wait(dat.day2.maintask.MRIinfo.TR*1000)
    end
    
    %     clearpict(1);
    
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
    %     clearpict(1);
    
end


results = []; results.SOT = [];
if MRIFlag == 1
    results.SOT.cgflip.t0_fix0 = t0_fix0;
    results.SOT.cgflip.t0_standby = t0_standby;
    results.SOT.raw.t0_fix0 = t0_fix0_raw;
    results.SOT.raw.trig_1st = trig_1st;
    results.SOT.raw.trig_5th = trig_5th;
end
eventmarker = 0;
breakmarker = 0;
rewamarker = -1;
i = 0; f1c = 0; sc = 0; rc = 0; cc = 0; dc = 0; nc = 0;
SOT_f1=[]; SOT_s=[]; SOT_r=[]; SOT_cue=[]; SOT_dot=[]; SOT_null=[];
SOT_cg_s=[];SOT_cg_f1=[];SOT_cg_r=[]; SOT_cg_cue=[]; SOT_cg_dot=[]; SOT_cg_null=[];
for trl = 1:size(config.stim.stimvec,2)
    
    %     try
    
    %     buff = 2;
    %     clearpict(1); % Clear buffer 1 (for fixation cross)
    %     clearpict(2); % Clear buffer 2 (for image presentation)
    %     clearpict(3); % Clear buffer 3 (for feedback cues)
    
    
    % ------------------------ present fix x -------------------------- %
    
    %     cd(paths.fb); loadpict('checker.bmp',1,0,0,bg_w,bg_h);
    %     setforecolour(.5,.5,.5);
    %     settextstyle('Arial', 300);
    %     preparestring('+', 1);  % Prepare buffer 1 (fixation cross)
    %
    %     f1c=f1c+1; SOT_f1(f1c,1)=toc; tic  % reset stopwatch
    %     drawpict(1)
    %     if eyeflag==1
    %         Eyelink('Message',[num2str(trl) 'f1'])  % send trialnumber and stimulustype to eyetracker
    %     end
    %     wait(config.timing.fixation1)
    %     eventmarker = eventmarker+1;
    %     results.presentation{eventmarker,1} = '+';
    %     clearpict(1);
    
    % ----------------------------------------------------------------- %
    
    
    % ------------------------ present break -------------------------- %
    
    if config.stim.stimvec(2,trl) == -1 % if break
        
        breakmarker = breakmarker + 1;
        
        categories(1,trl) = -1;
        temptrlord = categories;
        tempacc = results.accuracy;
        tempstimvec = stimvec_compl(4,:)';
        
        % calc points and money
        monie = sum(tempstimvec(1:trl-1) == 1);
        points = sum(tempstimvec(1:trl-1,1)== 1); %...
        maxpoints = numel(find(tempstimvec == 1));
        
        % conjure text
        %         monietext = [ 'Bisher haben Sie verdient:  ' num2str(monie) ' Euro von +6 Euro' ]; % 5 cent per point
        %         pointstxt = [ '(' num2str(points) ' Punkten von ' num2str(maxpoints) ' Punkten)' ];
        
        try
        % make a performance graph
        currents = points; maxes = maxpoints; graphpoints = [currents maxes];
        bargrphlabel={'Verdient';'Max Punkte'};
        figure('visible', 'off');
        bar(graphpoints);
        set(gca,'xticklabel',bargrphlabel);
        cd(paths.fb);
        export_fig perfgraph.bmp
        
        if MRIFlag == 1
            I = imread('perfgraph.bmp'); % flip the graph
            Ir = flipdim(I,2);
            imwrite(Ir,'perfgraph_mirror.bmp','bmp'); % save it again
        end
        
        cgpencol([1 1 1]);
        cgfont('Arial', 35);
        cd(paths.fb);
        if MRIFlag == 1
            cgloadbmp(1, 'instruction_break.bmp');
            cgloadbmp(2, 'perfgraph_mirror.bmp');
        else
            cgloadbmp(1, 'instruction_break_orig.bmp');
            cgloadbmp(2, 'perfgraph.bmp');
        end
        %         cgloadbmp(2, ['perfgraph_' num2str(breakmarker) '_' num2str(DayNumber) '.bmp']);
        
        catch
            if MRIFlag == 1
            cgloadbmp(1, 'instruction_break.bmp');
            else
            cgloadbmp(1, 'instruction_break_orig.bmp');
            end
        end
        cgsetsprite(0);
        cgdrawsprite(1,0,200);
        cgdrawsprite(2,0,-60);
        cgflip(0,0,0);
        wait(60000);
%         waitkeydown(60000,[RewKey,PunKey]);
%         waitkeydown(inf,[RewKey,PunKey]);
        %         wait(1);
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = 'break';
        %         delete perfgraph.bmp perfgraph_mirror.bmp
        
        
        % ----------------------------------------------------------------- %
        
        
        % ------------------------ present null --------------------------- %
        
    elseif config.stim.stimvec(2,trl) == 0 % if null
        
        i = i+1
        fprintf('\nTrial %d, null trial\n', i)
        
        clearkeys;
        
        %%%%%% load background %%%%%%
        cd(paths.fb); cgloadbmp(1,'checker_b.bmp',1152,864);
        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgfont('Arial',300)
        cgpencol(.5,.5,.5)
        cgtext('+',0,40);
        
        nc=nc+1; SOT_null(nc,1)=toc;
        SOT_cg_null(i) = cgflip(0,0,0);   % load fixation cross
%         if eyeflag==1
%             Eyelink('Message',[num2str(i) 'n'])  % send trialnumber and stimulustype to eyetracker
%         end
        wait(config.timing.fixation1(trl))
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = 'null';
        
        % record
        results.keypress(i,1) = 0;
        results.rt(i,1)       = 0;
        results.accuracy(i,1) = NaN;
        results.trl{i,1}      = 'null'; % filename
        results.trl{i,2}      = NaN; % category: reward always 1
        results.trl{i,3}      = 0; % contingencies
        results.trl{i,4}      = NaN; % memorabiliy
        
        % ----------------------------------------------------------------- %

        
        % ------------------------ present image -------------------------- %
        
    elseif config.stim.stimvec(2,trl) > 0
        
        i = i+1
        fprintf('\nTrial %d\n', i)
        
        clearkeys;
        
        
        %%%%%% load background %%%%%%
        cd(paths.fb);
        cgloadbmp(1,'checker_b.bmp',1152,864);
        
        %%%%%% load stimulus %%%%%%
        cd(paths.stim);
        cgloadbmp(2,char(config.stim.fname{i,1}),500,360);
        categories(1,i) = config.stim.stimvec(2,trl);
        
        %%%%%% draw stimulus screen %%%%%%
        sc=sc+1;
        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgdrawsprite(2,0,40);
        SOT_s(sc,1)=toc;
        SOT_cg_s(sc) = cgflip(0,0,0);
        if eyeflag==1
            Eyelink('Message',[num2str(i) 's'])  % send trialnumber and stimulustype to eyetracker
        end
        wait(config.timing.stim);
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = char(config.stim.fname{i,1});
        
        %%%%%% load Question screen %%%%%%
        cd(paths.fb);
        cgloadbmp(1,'checker_b.bmp',1152,864);
        if Schedules == 1
            if DayNumber == 1 % if it's inside day
                if MRIFlag == 1
                    questionSc = 'privatepublic.bmp';
                    answerSc_L = 'privatepublic_private.bmp';
                    answerSc_R = 'privatepublic_public.bmp';
                else
                    questionSc = 'privatepublic_orig.bmp';
                    answerSc_L = 'privatepublic_private_orig.bmp';
                    answerSc_R = 'privatepublic_public_orig.bmp';
                end
            elseif DayNumber == 2
                if MRIFlag == 1
                    questionSc = 'natureurban.bmp';
                    answerSc_L = 'natureurban_nature.bmp';
                    answerSc_R = 'privatepublic_urban.bmp';
                else
                    questionSc = 'natureurban_orig.bmp';
                    answerSc_L = 'natureurban_nature_orig.bmp';
                    answerSc_R = 'privatepublic_urban_orig.bmp';
                end
            end
            
        elseif Schedules == 2
            if DayNumber == 2 % if it's inside day
                if MRIFlag == 1
                    questionSc = 'privatepublic.bmp';
                    answerSc_L = 'privatepublic_private.bmp';
                    answerSc_R = 'privatepublic_public.bmp';
                else
                    questionSc = 'privatepublic_orig.bmp';
                    answerSc_L = 'privatepublic_private_orig.bmp';
                    answerSc_R = 'privatepublic_public_orig.bmp';
                end
            elseif DayNumber == 1
                if MRIFlag == 1
                    questionSc = 'natureurban.bmp';
                    answerSc_L = 'natureurban_nature.bmp';
                    answerSc_R = 'natureurban_urban.bmp';
                else
                    questionSc = 'natureurban_orig.bmp';
                    answerSc_L = 'natureurban_nature_orig.bmp';
                    answerSc_R = 'natureurban_urban_orig.bmp';
                end
            end
        end
        cgloadbmp(2,questionSc);
        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgdrawsprite(2,0,40); tmp = toc;
        responset = cgflip(0,0,0);
        if eyeflag==1
            Eyelink('Message',[num2str(i) 'r'])  % send trialnumber and stimulustype to eyetracker
        end
        waitkeydown(config.timing.stim,[RewKey,PunKey]);
        
        % RT
        [key, t, n] = getkeydown;
        if n == 0 %if they didn't press anything
            keypress = NaN; % mark their response as nan
            rt = NaN; % mark their reaction time as nan
            wait(0);
        else % otherwise
            %%%%%% load answer screen %%%%%%
            cgloadbmp(1,'checker_b.bmp',1152,864);
            if key(1) == 28
                cgloadbmp(2,answerSc_L);
            elseif key(1) == 29
                cgloadbmp(2,answerSc_R);
            end
            cgsetsprite(0);
            cgdrawsprite(1,0,0); cgdrawsprite(2,0,40); SOT_cg_r(i) = cgflip(0,0,0);
            keypress = key(1); % their response is whatever key they pressed.
            rt = t(1)/1000 - responset; rt = rt*1000; % and their reaction time is the time they pressed the button-the time the stimulus apprered
            rc=rc+1; SOT_r(rc,1) = tmp+rt; clear tmp
            wait(2000-rt);
        end
        
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = 'response';
        
        % accuracy
        if (keypress == RewKey && config.stim.stimvec(2,trl) == rewcond)
            accuracy = 1;
            fprintf('\nCorrect \n')
        elseif (keypress == PunKey && config.stim.stimvec(2,trl) == puncond)
            accuracy = 1;
            fprintf('\nCorrect \n')
        elseif (keypress ~= RewKey  && config.stim.stimvec(2,trl) == rewcond && ~isnan(keypress)) ...
                || (keypress ~= PunKey && config.stim.stimvec(2,trl) == puncond && ~isnan(keypress))
            accuracy = 0;
            fprintf('\nIncorrect \n')
        elseif isnan(keypress)
            accuracy = 0;
            fprintf('\nMissed \n')
        end
        
        % ----------------------------------------------------------------- %
        
        
        % record
        results.keypress(i,1) = keypress;
        results.rt(i,1)       = rt;
        results.accuracy(i,1) = accuracy;
        results.trl{i,1}      = char(config.stim.fname{i,1}); % filename
        results.trl{i,2}      = config.stim.stimvec(2,trl); % category: reward always 1
        results.trl{i,3}      = config.stim.stimvec(4,trl); % contingencies
        results.trl{i,4}      = config.stim.fname{i,3}; % memorabiliy
        
        %         categories    = cell2mat(results.trl(:,2));
        %         memorability  = cell2mat(results.trl(:,3));
        %         contingencies = cell2mat(results.trl(:,4));
        %         accuracies    = results.accuracy;
        %         RTs           = results.rt;
        %
        %         results.infotable = table(categories,memorability,contingencies,accuracies,RTs);
        
        fprintf('Current accuracy: %2.2f\n', mean(results.accuracy))
        
        
        % ------------------------- present dot --------------------------- %
        
        cd(paths.fb); cgloadbmp(1,'checker_b.bmp',1152,864);
        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgfont('Arial',200)
        cgpencol(.3,.3,.3)
        cgtext('.',0,50);
        
        dc=dc+1; SOT_dot(dc,1)=toc;
        SOT_cg_dot(i) = cgflip(0,0,0); % load fixation cross
        if eyeflag==1
            Eyelink('Message',[num2str(i) 'd'])  % send trialnumber and stimulustype to eyetracker
        end
        wait(config.timing.dot(i))
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = '.';
        
        % ----------------------------------------------------------------- %
        
        
        % ------------------------ present fback -------------------------- %
        
        %%%%%% load background %%%%%%
        cd(paths.fb);
        cgloadbmp(1,'checker_b.bmp',1152,864);
        
        %%%%%% load feedback %%%%%%
        cd(paths.fb);
        if config.stim.stimvec(4,trl) == 1
            if DayNumber == 1
                if MRIFlag == 1
                    cgloadbmp(2,'1euro.bmp',150,150);
                else
                    cgloadbmp(2,'1euro_orig.bmp',150,150);
                end
                rewamarker = rewamarker+1;
                
                cues(1,i) = 1;
                
                %%%%%% draw progress bar %%%%%
                if MRIFlag == 1
                    cgloadbmp(3,'progress_bar.bmp',272,42);
                else
                    cgloadbmp(3,'progress_bar_orig.bmp',272,42);
                end
                cgloadbmp(4,'marker.bmp')
                %                 cgmakesprite(3,config.numtrials.rew*8*0.3,10)
                %                 cgsetsprite(3)
                %                 cgpencol(0.5,0.5,0.5) % pen colour
                %                 cgrect(0,0,config.numtrials.rew*8*0.3,10) % first draw max bar
                %                 cgpencol(0.5,0.8,0.5) % pen colour
                %                 prg_num = sum(cues==1)*10*0.3;
                %                 cgrect(0,0,prg_num,30) % now the progress
                %                 cgpencol(0.5,0.5,0.5) % pen colour
                
                %%%%%% draw feedback screen %%%%%%
                cc=cc+1;
                cgsetsprite(0);
                cgdrawsprite(1,0,0);
                cgdrawsprite(2,0,40);
                cgdrawsprite(3,0,-60)
                if rewamarker <= 60
                    if MRIFlag == 1
                        cgdrawsprite(4,123-4.725*rewamarker,-62) %%% +123 is the start position
%                         cgdrawsprite(4,123-6.3*rewamarker,-62) %%% +123 is the start position
                    else
                        cgdrawsprite(4,-123+4.725*rewamarker,-62) %%% -123 is the start position
%                         cgdrawsprite(4,-123+6.3*rewamarker,-62) %%% -123 is the start position
                    end
                elseif rewamarker > 60
                    if MRIFlag == 1
                        cgdrawsprite(4,123-6.3*40,-62)
                    else
                        cgdrawsprite(4,-123+6.3*40,-62) %%% -123 is the start position
                    end
                end
                SOT_cue(cc,1)=toc;
                SOT_cg_cue(cc) = cgflip(0,0,0);
                if eyeflag==1
                    Eyelink('Message',[num2str(i) 'c'])  % send trialnumber and stimulustype to eyetracker
                end
                wait(config.timing.cue);
                eventmarker = eventmarker+1;
                results.presentation{eventmarker,1} = char(['feedback_' num2str(cues(1,i))]);
                cgfreesprite(2)
                
            else
                cues(1,i) = 1;
                if MRIFlag == 1
                    cgloadbmp(2,'5cent.bmp',150,150);
                else
                    cgloadbmp(2,'5cent_orig.bmp',150,150);
                end
                cgloadbmp(3,'progress_bar_grey.bmp',272,42);
                %%%%%% draw progress bar %%%%%
                %                 cgmakesprite(3,config.numtrials.rew*8*0.3,30)
                %                 cgsetsprite(3)
                %                 cgpencol(0.4,0.4,0.4) % pen colour
                %                 cgrect(0,0,config.numtrials.rew*8*0.3,30) % greyed out bar
                
                %%%%%% draw feedback screen %%%%%%
                cc=cc+1;
                cgsetsprite(0);
                cgdrawsprite(1,0,0);
                cgdrawsprite(2,0,40);
                cgdrawsprite(3,0,-60)
                %         cgdrawsprite(4,-123,-62) %%% -123 is the start position
                SOT_cue(cc,1)=toc;
                SOT_cg_cue(cc) = cgflip(0,0,0);
                if eyeflag==1
                    Eyelink('Message',[num2str(i) 'c'])  % send trialnumber and stimulustype to eyetracker
                end
                wait(config.timing.cue);
                eventmarker = eventmarker+1;
                results.presentation{eventmarker,1} = char(['feedback_' num2str(cues(1,i))]);
                cgfreesprite(2)
                
            end
            
        else
            if MRIFlag == 1
                cgloadbmp(2,'zerocoin.bmp',150,150);
            else
                cgloadbmp(2,'zerocoin_orig.bmp',150,150);
            end
            cues(1,i) = 0;
            
            %%%%%% draw progress bar %%%%%
            cgloadbmp(3,'progress_bar_grey.bmp',272,42);
            %             cgmakesprite(3,config.numtrials.rew*8*0.3,30)
            %             cgsetsprite(3)
            %             cgpencol(0.4,0.4,0.4) % pen colour
            %             cgrect(0,0,config.numtrials.rew*8*0.3,30) % greyed out bar
            
            %%%%%% draw feedback screen %%%%%%
            cc=cc+1;
            cgsetsprite(0);
            cgdrawsprite(1,0,0);
            cgdrawsprite(2,0,40);
            cgdrawsprite(3,0,-60)
            %         cgdrawsprite(4,-123,-62) %%% -123 is the start position
            SOT_cue(cc,1)=toc;
            SOT_cg_cue(cc) = cgflip(0,0,0);
            if eyeflag==1
                Eyelink('Message',[num2str(i) 'c'])  % send trialnumber and stimulustype to eyetracker
            end
            wait(config.timing.cue);
            eventmarker = eventmarker+1;
            results.presentation{eventmarker,1} = char(['feedback_' num2str(cues(1,i))]);
            cgfreesprite(2)
            
        end
        
        % ----------------------------------------------------------------- %
        
        
        
    end % close the presentation conditional
    
    % ------------------------ present fix x -------------------------- %
    
    
    cd(paths.fb); cgloadbmp(1,'checker_b.bmp',1152,864);
    cgsetsprite(0);
    cgdrawsprite(1,0,0);
    cgfont('Arial',300)
    cgpencol(.5,.5,.5)
    cgtext('+',0,40);
    
    f1c=f1c+1; SOT_f1(f1c,1)=toc;
    SOT_cg_f1(i) = cgflip(0,0,0);             % load fixation cross
    if eyeflag==1
        Eyelink('Message',[num2str(i) 'f'])  % send trialnumber and stimulustype to eyetracker
    end
    wait(config.timing.fixation1(trl))
    eventmarker = eventmarker+1;
    results.presentation{eventmarker,1} = '+';
    
    % ----------------------------------------------------------------- %
    
    
    % calculate cumulative SOT
    %     results.SOT.cumulative.fix  = cumsum(SOT_f1);
    %     results.SOT.cumulative.stim = cumsum(SOT_s);
    %     results.SOT.cumulative.resp = cumsum(SOT_r);
    results.SOT.raw.fix  = SOT_f1;
    results.SOT.raw.stim = SOT_s;
    results.SOT.raw.resp = SOT_r;
    results.SOT.raw.cue  = SOT_cue;
    results.SOT.raw.null = SOT_null;
    results.SOT.cgflip.stim = SOT_cg_s;
    results.SOT.cgflip.resp = SOT_cg_r;
    results.SOT.cgflip.fix  = SOT_cg_f1;
    results.SOT.cgflip.cue  = SOT_cg_cue;
    results.SOT.cgflip.null = SOT_cg_null;
    
    % intermediate save
    if DayNumber == 1
        dat.day1.maintask.results = results;
        save([paths.behav num2str(ID) '_' num2str(DayNumber) '_backup.mat'],'dat')
    elseif DayNumber == 2
        dat.day2.maintask.results = results;
        save([paths.behav num2str(ID) '_' num2str(DayNumber) '_backup.mat'],'dat')
    end
    %     catch
    %         eventmarker = eventmarker+1;
    %         results.stimord{eventmarker,1} = 'ERROR';
    %     end
    
end % close the trial loop

fprintf('\nTask Done\n')

% compute the final amount
monie2 = (sum(results.accuracy==1 & config.stim.stimvec(2,find(stimvec_mod(1,:)~=777))'==1)*0.05) + ...
    (sum(results.accuracy==1 & config.stim.stimvec(2,find(stimvec_mod(1,:)~=777))' == 3)*-0.05) + ...
    (sum(results.accuracy==1)*0.0125);
monietext2 = ['Der endgueltige Geldbetrag: ' num2str(monie2)];

% terminate protocol
cd(paths.fb)
if MRIFlag == 1
    cgloadbmp(1,'instruction_ende.bmp')
else
    cgloadbmp(1,'instruction_ende_orig.bmp');
end
cgsetsprite(0);
cgdrawsprite(1,0,0);
cgflip(0,0,0); wait(3000);
if MRIFlag == 1
    cgloadbmp(1,'eyes.bmp');
    cgsetsprite(0);
    cgdrawsprite(1,0,0);
    cgflip(0,0,0); wait(6000);
end

cgshut;

%% wrap up
% calculate final cumulative SOT
% results.SOT.cumulative.fix  = cumsum(SOT_f1);
% results.SOT.cumulative.stim = cumsum(SOT_s);
% results.SOT.cumulative.resp = cumsum(SOT_r);
results.SOT.raw.fix  = SOT_f1;
results.SOT.raw.stim = SOT_s;
results.SOT.raw.resp = SOT_r;
results.SOT.raw.cue  = SOT_cue;
results.SOT.raw.null = SOT_null;
results.SOT.cgflip.stim = SOT_cg_s;
results.SOT.cgflip.resp = SOT_cg_r;
results.SOT.cgflip.fix  = SOT_cg_f1;
results.SOT.cgflip.cue  = SOT_cg_cue;
results.SOT.cgflip.null = SOT_cg_null;
    
fprintf('\nwrapped up\n')

% save results
if DayNumber == 1
    fprintf('\n\nBEHAVIOURAL DATA SAVING...PLEASE WAIT\n\n')
    dat.day1.maintask.results = results;
    dat.day1.maintask.labels  = {'categories->Nature(1)/Urban(2)/Private(3)/Public(4)';'contingencies->rewarded(1),neutral(2)'; ...
        'results.trl: filenames, reward(1)/neutral(2), contingencies(actually rewarded: 1, not rewarded 2, neutral 0), memorability'};
    dat.day1.maintask.config.keymap = map;
    save([paths.behav num2str(ID) '_' num2str(DayNumber) '.mat'],'dat')
    fprintf('\n data day 1 saved \n')
elseif DayNumber == 2
    fprintf('\n\nBEHAVIOURAL DATA SAVING...PLEASE WAIT\n\n')
    dat.day2.maintask.results = results;
    dat.day2.maintask.labels  = {'categories->Nature(1)/Urban(2)/Private(3)/Public(4)';'contingencies->rewarded(1),neutral(2)'; ...
        'results.trl: filenames, reward(1)/neutral(2), contingencies(actually rewarded: 1, not rewarded 2, neutral 0), memorability'};
    dat.day2.maintask.config.keymap = map;
    save([paths.behav num2str(ID) '_' num2str(DayNumber) '.mat'],'dat')
    fprintf('\n data day 2 saved \n')
    
end



fprintf(monietext2);
if EyeFlag == 1
    fprintf('\n\nEYETRACKING DATA SAVING...PLEASE WAIT\n\n')
end

% fprintf('\n%%%%%%%%%%%%%%%%%%%%%%\nTotal money earned: %2.3f Euro\n%%%%%%%%%%%%%%%%%%%%%%\n', monie2)
end