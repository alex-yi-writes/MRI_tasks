%% PRACTICE TASK : Scene (Indoor/Outdoor) classification task
%  In emergency contact Alex: 
%              +49 162 713 20 49
%              yeo-jin.yi.15@ucl.ac.uk
% See also StarterScript_MRIpilot, runMainTask_MRIpilot,
%       runMemoryTest_MRIpilot

%%  work log

%   09_10_2018 created the script

%% START

function [dat] = runPracticeTask_3Tpilot(ID,AgeGroup,EyeFlag,CnterBal,eyetracker_filename,DayNumber)

global dat eyeflag

%% basic setups

% initialize eyetracking
eyeflag=EyeFlag;
% fname_eyetracker = ['pr' num2str(ID) '.edf'];
fname_eyetracker = eyetracker_filename;

if eyeflag==1
    
    sca;
    
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


% set paths
if DayNumber == 1
    if CnterBal == 1 % if private = reward
        path.stim_rew = '\\fs-md\users\yiy\Documents\DOPE\pilot_3T\practicestim\practicestim_private\';
        path.stim_pun = '\\fs-md\users\yiy\Documents\DOPE\pilot_3T\practicestim\practicestim_public\';
        rewcond = 1; puncond = 2; rewchar = 'v'; punchar = 'b';
    elseif CnterBal == 2 % if public == reward
        path.stim_rew = '\\fs-md\users\yiy\Documents\DOPE\pilot_3T\practicestim\practicestim_public\';
        path.stim_pun = '\\fs-md\users\yiy\Documents\DOPE\pilot_3T\practicestim\practicestim_private\';
        rewcond = 2; puncond = 1; rewchar = 'b'; punchar = 'v';
    end
    
elseif DayNumber == 2
    if CnterBal == 1 % if nature = reward
        path.stim_rew = '\\fs-md\users\yiy\Documents\DOPE\pilot_3T\practicestim\practicestim_nature\';
        path.stim_pun = '\\fs-md\users\yiy\Documents\DOPE\pilot_3T\practicestim\practicestim_urban\';
        rewcond = 1; puncond = 2; rewchar = 'n'; punchar = 'u';
    elseif CnterBal == 2 % if urban == reward
        path.stim_rew = '\\fs-md\users\yiy\Documents\DOPE\pilot_3T\practicestim\practicestim_urban\';
        path.stim_pun = '\\fs-md\users\yiy\Documents\DOPE\pilot_3T\practicestim\practicestim_nature\';
        rewcond = 2; puncond = 1; rewchar = 'u'; punchar = 'n';
    end
end
path.parent   = '\\fs-md\users\yiy\Documents\DOPE\';
path.stim     = [ path.parent 'pilot_3T\stim_resized_IO_sml_new\' ];
path.fb       = [ path.parent 'pilot_3T\feedback\' ];
path.pupil    = [ path.parent 'pilot_3T\raw\pupil\' ];
path.behav    = [ path.parent 'pilot_3T\raw\behav\' ];
% paths.cogent   = 'C:\Studies\Grid_experiment2\Cogent2000v1.32';
path.cogent   = '\\fs-md\users\yiy\Documents\MATLAB\Cogent2000v1.33';

addpath(genpath(path.cogent))


% create the data structure
cd(path.behav)
if exist([num2str(ID) '_' num2str(DayNumber) '.mat']) == 2
    error('\nFILE EXISTS!!! Move the .mat file out of the folder and try again.\n')
else
    dat = [];
    dat.ID       = ID;
    dat.AgeGroup = AgeGroup;
    dat.practice = [];
    if DayNumber == 1
        if CnterBal == 1
            dat.day1.RewardCategory = 'private(3)';
            fprintf('\nCounterbalance marker: 1\n')
        else
            dat.day1.RewardCategory = 'public (4)';
            fprintf('\nCounterbalance marker: 2\n')
        end
    elseif DayNumber == 2
        if CnterBal == 1
            dat.day2.RewardCategory = 'nature (1)';
            fprintf('\nCounterbalance marker: 1\n')
        else
            dat.day2.RewardCategory = 'urban  (2)';
            fprintf('\nCounterbalance marker: 2\n')
        end
    end
    dat.practice.ETinfo  = fname_eyetracker;
    
end

%% experimental setups

% timings
config.timing.outcome  = 2500;
fixjit1 = repmat(500:50:700,1,48); fixjit2 = repmat(500:50:700,1,48);
config.timing.fixation1= fixjit1(randperm(length(fixjit1)));
config.timing.fixation2= fixjit2(randperm(length(fixjit2)));
config.timing.cue      = 2500;

% --- for breezeing through --- %
% config.timing.fixation1= zeros(1,240);
% config.timing.fixation2= zeros(1,240);
% config.timing.cue      = 1;
% config.timing.outcome  = 1;
% ----------------------------- %

% trials
config.numtrials.total = 40;
if DayNumber == 1
    config.numtrials.rew   = 30;
    config.numtrials.pun   = 10;
elseif DayNumber == 2
    config.numtrials.rew   = 10;
    config.numtrials.pun   = 30;
end

% contingency vector
% trlvec    = randperm(config.numtrials.total);
Rew = ones(1,config.numtrials.rew*.8); nRew = ones(1,config.numtrials.rew*.2).*2; rewmat_rewarded = horzcat(Rew,nRew); % 1 for reward 2 for nothing
Pun = ones(1,config.numtrials.pun*.8)*3; nPun = ones(1,config.numtrials.pun*.2).*4; rewmat_punished = horzcat(Pun,nPun); % 3 for punishment 4 for nothing
rewmat_rewarded = rewmat_rewarded(randperm(length(rewmat_rewarded)));
rewmat_punished = rewmat_punished(randperm(length(rewmat_punished)));

% roster
stimlist = [];
cd(path.stim_rew)
temp = dir('*.bmp'); [stimlist.rew{1:config.numtrials.rew}] = deal(temp(1:config.numtrials.rew).name); clear temp
cd(path.stim_pun)
temp = dir('*.bmp'); [stimlist.pun{1:config.numtrials.pun}] = deal(temp(1:config.numtrials.pun).name); clear temp
stimlist.all = horzcat(stimlist.rew,stimlist.pun);

flist = stimlist.all;
flist = flist(randperm(length(stimlist.all)));
stimlist.all = flist;

rc = 0; nc = 0; stimvec = [];
for i = 1:config.numtrials.total
    if strcmp(flist{i}(1,1),rewchar)
        rc = rc+1;
        stimvec(1,i) = i;
        stimvec(2,i) = rewmat_rewarded(rc);
    elseif strcmp(flist{i}(1,1),punchar)
        nc = nc+1;
        stimvec(1,i) = i;
        stimvec(2,i) = rewmat_punished(nc);
    end
end
config.trials = stimvec;

dat.practice.config = config; % save the configuration

fprintf('\n\nTASK SETUP DONE\n\n')
fprintf('\n\nINITIALISING PRACTICE TASK\n\n')

%% run

% configs
config_display(1, 6, [0 0 0], [1 1 1], 'Arial', 20, 4); % 1280*1024
config_sound;
config_keyboard;

monitor = get(0,'screensize');
bg_w = monitor(3); bg_h = monitor(4);

%%%%%%%%%%%%
start_cogent
%%%%%%%%%%%%

map         = getkeymap; % define keyboard IDs
if CnterBal == 1 % if private-nature reward
    RewKey  = map.A;
    PunKey  = map.L;
elseif CnterBal == 2
    RewKey  = map.L;
    PunKey  = map.A;
end
buff        = 2;

instr{1,1}  = 'Willkommen zu unserem Experiment!';
instr{2,1}  = 'Bitte klassifizieren Sie die Art des Bildes auf dem Bildschirm.';
instr{3,1}  = 'Wenn Sie bereit sind, druecken Sie die Leertaste, um zu beginnen.';

settextstyle('Arial', 35);
preparestring(instr{1,1},1)
preparestring(instr{2,1},2,0,50)
preparestring(instr{3,1},2,0,20)

if eyeflag==1
    Eyelinkerror = Eyelink('StartRecording'); % start recording to the file
    if Eyelinkerror ~= 0
        return
        error ('Eyetracker failed')
    end
    Eyelink('Message',num2str([0 100]))
end

drawpict(1); wait(500);
drawpict(2); waitkeydown(inf,map.Space);

for trl = 1:size(stimvec,2)
    
    buff = 2;
    clearpict(1); % Clear buffer 1 (for fixation cross)
    clearpict(2); % Clear buffer 2 (for image presentation)
    clearpict(3); % Clear buffer 3 (for feedback cues)
    
    
    % ------------------------ present fix x -------------------------- %
    
    cd(path.fb); loadpict('checker.bmp',1,0,0,bg_w,bg_h);
    setforecolour(.5,.5,.5);
    settextstyle('Arial', 300);
    preparestring('+', 1);  % Prepare buffer 1 (fixation cross)
    
    drawpict(1)             % load fixation cross
    if eyeflag==1
        Eyelink('Message',[num2str(trl) 'f1'])  % send trialnumber and stimulustype to eyetracker
    end
    wait(config.timing.fixation1)
    clearpict(1);
    
    % ----------------------------------------------------------------- %
    
    
    % ------------------------ present image -------------------------- %
    
    fprintf('\nTrial %d\n', trl)
    
    clearkeys;
    cd(path.fb); loadpict('checker.bmp',1,0,0,bg_w,bg_h);loadpict('checker.bmp',2,0,0,bg_w,bg_h);
    setforecolour(.9,.9,.9);
    settextstyle('Arial', 80);
    if DayNumber == 1
        preparestring('Privat        Oeffentlich',2);
    elseif DayNumber == 2
        preparestring('Natur        Staedtisch',2);
    end
    
    if strcmp(stimlist.all{trl}(1,1),rewchar)
        cd(path.stim_rew)
        loadpict(char(stimlist.all(trl)),1);
    else
        cd(path.stim_pun)
        loadpict(char(stimlist.all(trl)),1);
    end
    fname = char(stimlist.all(trl));
    
    drawpict(1);
    if eyeflag==1
        Eyelink('Message',[num2str(trl) 's'])  % send trialnumber and stimulustype to eyetracker
    end
    wait(config.timing.cue);
    
    stimonset = drawpict(2);
    if eyeflag==1
        Eyelink('Message',[num2str(trl) 'r'])  % send trialnumber and stimulustype to eyetracker
    end
    waitkeydown(config.timing.cue,[RewKey,PunKey]);
    clearpict(1);clearpict(2);
    
    [key, t, n] = getkeydown;
    if n == 0; %if they didn't press anything
        response = NaN; % mark their response as nan
        rt = NaN; % mark their reaction time as nan
    else % otherwise
        response = key(1); % their response is whatever key they pressed.
        rt = t(1) - stimonset; % and their reaction time is the time they pressed the button-the time the stimulus apprered
    end
    recordresponse.keypress(trl) = response;
    
    % accuracy
    if (response == RewKey &&  strcmp(fname(1,1),rewchar))
        accuracy = 1;
        fprintf('\nCorrect \n')
    elseif (response == PunKey && strcmp(fname(1,1),punchar))
        accuracy = 1;
        fprintf('\nCorrect \n')
    elseif (response == PunKey  && strcmp(fname(1,1),rewchar) && ~isnan(response)) ...
            || (response == RewKey && strcmp(fname(1,1),punchar) && ~isnan(response))
        accuracy = 0;
        fprintf('\nIncorrect \n')
    elseif isnan(response)
        accuracy = 0;
        fprintf('\nMissed \n')
    end
    recordresponse.accu(trl) = accuracy;
    
    % ----------------------------------------------------------------- %
    
    
    % ---------------------- present feedback -------------------------- %
    
    if  stimvec(2,trl) == 1 && accuracy == 1 % if it is indeed rewarded and gave correct answer
        
        cd(path.fb); loadpict('checker.bmp',1,0,0,bg_w,bg_h);loadpict('checker.bmp',2,0,0,bg_w,bg_h);
        setforecolour(.5,.5,.5);
        settextstyle('Arial', 300);
        preparestring('+', 1);  % Prepare buffer 1 (fixation cross)
        loadpict('50plus.bmp',2,0,0,bg_w/6,bg_h/4);
        
        drawpict(1);
        if eyeflag==1
            Eyelink('Message',[num2str(trl) 'f2'])  % send trialnumber and stimulustype to eyetracker
        end
        wait(config.timing.fixation2);
        
        drawpict(2);
        if eyeflag==1
            Eyelink('Message',[num2str(trl) 'f'])  % send trialnumber and stimulustype to eyetracker
        end
        wait(config.timing.outcome) ;
        
        clearpict(1); clearpict(2);
        
    elseif (stimvec(2,trl) == 4 && accuracy == 1) || (stimvec(2,trl) == 2 && accuracy == 1) % if bad luck with reward cat or good luck with norew cat
        
        cd(path.fb); loadpict('checker.bmp',1,0,0,bg_w,bg_h);loadpict('checker.bmp',2,0,0,bg_w,bg_h);
        
        setforecolour(.5,.5,.5);
        settextstyle('Arial', 300);
        preparestring('+', 1);  % Prepare buffer 1 (fixation cross)
        loadpict('zerocoin.bmp',2,0,0,bg_w/6,bg_h/4);
        
        drawpict(1);
        if eyeflag==1
            Eyelink('Message',[num2str(trl) 'f2'])  % send trialnumber and stimulustype to eyetracker
        end
        wait(config.timing.fixation2);
        
        drawpict(2);
        if eyeflag==1
            Eyelink('Message',[num2str(trl) 'f'])  % send trialnumber and stimulustype to eyetracker
        end
        wait(config.timing.outcome) ;
        
        clearpict(1); clearpict(2);
        
    elseif (stimvec(2,trl) == 3 && accuracy == 1) % if it is punished condition and correct
        
        cd(path.fb); loadpict('checker.bmp',1,0,0,bg_w,bg_h);loadpict('checker.bmp',2,0,0,bg_w,bg_h);
        setforecolour(.5,.5,.5);
        settextstyle('Arial', 300);
        preparestring('+', 1);  % Prepare buffer 1 (fixation cross)
        loadpict('zerocoin.bmp',2,0,0,bg_w/6,bg_h/4); %%%%% changed here!
        
        drawpict(1);
        if eyeflag==1
            Eyelink('Message',[num2str(trl) 'f2'])  % send trialnumber and stimulustype to eyetracker
        end
        wait(config.timing.fixation2);
        
        drawpict(2);
        if eyeflag==1
            Eyelink('Message',[num2str(trl) 'f'])  % send trialnumber and stimulustype to eyetracker
        end
        wait(config.timing.outcome) ;
        
        clearpict(1); clearpict(2);
        
    elseif accuracy == 0
        
        cd(path.fb); loadpict('checker.bmp',1,0,0,bg_w,bg_h);loadpict('checker.bmp',2,0,0,bg_w,bg_h);
        setforecolour(.5,.5,.5);
        settextstyle('Arial', 300);
        preparestring('+', 1);  % Prepare buffer 1 (fixation cross)
        loadpict('20minus2.bmp',2,0,0,bg_w/6,bg_h/4);
        
        drawpict(1);
        if eyeflag==1
            Eyelink('Message',[num2str(trl) 'f2'])  % send trialnumber and stimulustype to eyetracker
        end
        wait(config.timing.fixation2);
        
        drawpict(2);
        if eyeflag==1
            Eyelink('Message',[num2str(trl) 'f'])  % send trialnumber and stimulustype to eyetracker
        end
        wait(config.timing.outcome) ;
        
        clearpict(1); clearpict(2);
        
    end % close the feedback loop
    
    % ----------------------------------------------------------------- %

    
    recordstim.fname{trl,1} = fname;
    recordstim.reward(trl,1) = stimvec(2,trl);
    
    % intermediate save
    dat.practice.results = recordresponse;
    dat.practice.stim    = recordstim;
    save([path.behav num2str(ID) '_' num2str(DayNumber) '_backup.mat'],'dat')

end % close the trial loop

% wrap up
dat.practice.results = recordresponse;
dat.practice.stim    = recordstim;
save([path.behav num2str(ID) '_' num2str(DayNumber) '.mat'],'dat')

% terminate protocol
clearpict(1)
settextstyle('Arial', 35);
preparestring('Das Ende der Uebungsaufgabe',1,0,0)
drawpict(1); wait(3000);
stop_cogent;


end