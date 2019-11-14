%% EMOTIONAL MEMORY TASK : Scene classification task with emotional/neutral images
%  In emergency contact Alex:
%              +49 162 713 20 49
%              yeo-jin.yi.15@ucl.ac.uk

%%  work log

%   23_01_2019 created the script

%% START

function [dat] = runEmoMemTask(ID,AgeGroup,EyeFlag,MRIFlag,eyetracker_filename)

global eyeflag

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
paths.parent   = 'Z:\Documents\EmoMem\';
paths.supports = [paths.parent 'stim\supports\'];
paths.stim     = [paths.parent 'stim\old\'];
paths.pupil    = [paths.parent 'data\pupil\'];
paths.behav    = [paths.parent 'data\behav\'];
paths.cogent   = '\\fs-md\users\yiy\Documents\MATLAB\Cogent2000v1.33';
addpath(genpath(paths.cogent))

% create the data structure
cd(paths.behav)
if exist([num2str(ID) '_EM.mat']) == 2
else
    dat = [];
    dat.ID       = ID;
    dat.AgeGroup = AgeGroup;
    dat.emomem = [];
    dat.emomem.task.ETinfo  = fname_eyetracker;
end
load([paths.parent 'scripts\EM_GA_designstruct.mat'])


% scanner preparation
if MRIFlag == 1
    dat.emomem.task.MRIinfo   = [];
    dat.emomem.task.MRIinfo.scanPort = 1;
    dat.emomem.task.MRIinfo.dummy    = 5;        % no. of dummy vols
    dat.emomem.task.MRIinfo.nslice   = 51;       % no. of slices
    dat.emomem.task.MRIinfo.TE       = 32;       % time for each slice in msec
    dat.emomem.task.MRIinfo.TR       = 3.6; % MRIinfo.nslice * MRIinfo.TE/1000; % time for volume in msec
    dat.emomem.task.MRIinfo.nvols    = ceil(MRIinfo.dummy + 5 + ...  % no. of vols per run (5 dummy + 400 + 5 post)
        (options.scanblocklength * options.meanTrialLength) / MRIinfo.TR);
    dat.emomem.task.MRIinfo.total_slices = MRIinfo.nvols*MRIinfo.nslice;    % per run!
    %     [MRIinfo] = runScannerSetup_MRIpilot; % redundant
end

%% experimental setups

% condition markers
in_emo = 11; in_neu = 12;
out_emo = 21; out_neu = 22;

% timings
% fixjit = repmat(500:50:700,1,24);
config.timing.fixation = [(design_struct.eventlist(:,4)+design_struct.eventlist(:,5)-2)'];
config.timing.stim     = 2500;
config.timing.response = 2000; % maximum waiting time for the button press

% --- for breezeing through --- %
% config.timing.fixation1= zeros(1,240);
% config.timing.fixation2= zeros(1,240);
% config.timing.cue      = 1;
% config.timing.outcome  = 1;
% ----------------------------- %

% trials
config.numtrials.total   = 175;
config.numtrials.in_emo  = 35;
config.numtrials.in_neu  = 35;
config.numtrials.out_emo = 35;
config.numtrials.out_neu = 35;

% contingency vector
trlvec          = 1:config.numtrials.total;
contingencies   = design_struct.eventlist(:,3)';
indoors         = contingencies(find(contingencies<13));
outdoors        = contingencies(find(contingencies>20));
config.stim.indoors  = indoors;
config.stim.outdoors = outdoors;


% roster
% emotional-indoor
fname = [];
cd(paths.stim);cd emotional;cd in;
tmp = dir('*.bmp');
[fname.in_emo{1:config.numtrials.in_emo}] = deal(tmp(1:config.numtrials.in_emo).name); clear tmp;
% emotional-outdoor
cd(paths.stim);cd emotional;cd out;
tmp = dir('*.bmp');
[fname.out_emo{1:config.numtrials.out_emo}] = deal(tmp(1:config.numtrials.out_emo).name); clear tmp;
% neutral-indoor
cd(paths.stim);cd neutral;cd in;
tmp = dir('*.bmp');
[fname.in_neu{1:config.numtrials.in_neu}] = deal(tmp(1:config.numtrials.in_neu).name); clear tmp;
% neutral-outdoor
cd(paths.stim);cd neutral;cd out;
tmp = dir('*.bmp');
[fname.out_neu{1:config.numtrials.out_neu}] = deal(tmp(1:config.numtrials.out_neu).name); clear tmp;


% schedule presentation
tmpvec       = trlvec;
stimvec_temp = [];
inemo = 0; outemo = 0; inneu = 0; outneu = 0;
for sc = 1:length(contingencies)
    if contingencies(sc) == in_emo
        inemo = inemo+1;
        stimvec_temp(1,sc) = tmpvec(sc);
        stimvec_temp(2,sc) = contingencies(sc);
        stimvec_fname{1,sc}= fname.in_emo{inemo};
    elseif contingencies(sc) == out_emo
        outemo = outemo+1;
        stimvec_temp(1,sc) = tmpvec(sc);
        stimvec_temp(2,sc) = contingencies(sc);
        stimvec_fname{1,sc}= fname.out_emo{outemo};
    elseif contingencies(sc) == in_neu
        inneu = inneu+1;
        stimvec_temp(1,sc) = tmpvec(sc);
        stimvec_temp(2,sc) = contingencies(sc);
        stimvec_fname{1,sc}= fname.in_neu{inneu};
    elseif contingencies(sc) == out_neu
        outneu = outneu+1;
        stimvec_temp(1,sc) = tmpvec(sc);
        stimvec_temp(2,sc) = contingencies(sc);
        stimvec_fname{1,sc}= fname.out_neu{outneu};
    elseif contingencies(sc) == 0
        stimvec_temp(1,sc) = tmpvec(sc);
        stimvec_temp(2,sc) = contingencies(sc);
        stimvec_fname{1,sc}= 'null';
    end
end
stimvec = stimvec_temp; clear stimvec_temp

config.stim.stimvec  = stimvec;
config.stim.stimlist = stimvec_fname;
config.stim.filelist = fname;

fprintf('\n\nTASK SETUP DONE')
fprintf('\nINITIALISING EMOMEM TASK\n\n')


%% run

% configs
config_display(1, 5, [0 0 0], [1 1 1], 'Arial', 20, 4); % 1024x768
config_sound;
config_keyboard;

% scanner input
if MRIFlag == 1
    config_serial(MRIinfo.scanPort); % links to scanner for waitslice
end

monitor = get(0,'screensize');
bg_w = monitor(3); bg_h = monitor(4);

%%%%%%%%%%%%
start_cogent
%%%%%%%%%%%%

map     = getkeymap; % define keyboard IDs
InKey   = map.L;
OutKey  = map.A;
VeryKey = map.L;
NVeryKey= map.A;

buff    = 2;

% print instructions
cd(paths.supports)
loadpict('emomem_inst_task.png',1,0,0);
drawpict(1); waitkeydown(inf,map.Space);

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
    
    % hail the operator
    clearpict(1); clearpict(2);
    preparestring('The operator should initiate the task when ready.',1,0,0)
    drawpict(1)
    waitkeydown(inf,map.Space)
    
    % start the dummy scan
    scannerinput = 5;
    [k, dat.maintask.MRIinfo.dummyscanfirst] = waitkeydown(inf,scannerinput); % the first dummy scan
    for dum = 1:(dat.maintask.MRIinfo.dummy - 2)
        waitkeydown(inf,scannerinput);
    end
    [k, dat.maintask.MRIinfo.dummyscanlast] = waitkeydown(inf,scannerinput); % the last dummy scan
    
    cd(paths.stim); loadpict('checker.bmp',1,0,0,bg_w,bg_h);
    setforecolour(.5,.5,.5);
    settextstyle('Arial', 300);
    preparestring('+', 1);  % Prepare buffer 1 (fixation cross)
    drawpict(1)             % load fixation cross
    if eyeflag==1
        Eyelink('Message',[num2str(0) 'f0'])  % send trialnumber and stimulustype to eyetracker
    end
    wait(dat.emomem.task.MRIinfo.TR*1000)
    clearpict(1);
    
else
    cd(paths.parent); loadpict('checker.bmp',1,0,0,bg_w,bg_h);
    setforecolour(.5,.5,.5);
    settextstyle('Arial', 300);
    preparestring('+', 1);  % Prepare buffer 1 (fixation cross)
    drawpict(1)             % load fixation cross
    if eyeflag==1
        Eyelink('Message',[num2str(0) 'f0'])  % send trialnumber and stimulustype to eyetracker
    end
    wait(1000)
    %     clearpict(1);
    
end

results = []; results.SOT = [];
eventmarker = 0;
i = 0; f1c = 0; sc = 0; rc1 = 0; rc2 = 0; vc = 0;
SOT_f1=[]; SOT_s=[]; SOT_r=[]; SOT_v=[];

for trl = 1:size(stimvec,2)
    
    tic % ding ding ding
    
    buff = 5;
    clearpict(1); % pic
    clearpict(2); % resp
    clearpict(3); % resp received
    
    
    % ------------------------ present image -------------------------- %
    i = i+1;
    fprintf('\nTrial %d\n', trl)
    
    clearkeys;
    cd(paths.supports); loadpict('checker.bmp',1,0,0,bg_w,bg_h); loadpict('checker.bmp',2,0,0,bg_w,bg_h); loadpict('checker.bmp',3,0,0,bg_w,bg_h);
    %loadpict('checker.bmp',4,0,0,bg_w,bg_h);loadpict('checker.bmp',5,0,0,bg_w,bg_h);
    setforecolour(.9,.9,.9);
    settextstyle('Arial', 80);
    loadpict('inout.bmp',2,0,0);
    loadpict('inout_received.bmp',3,0,0);
    
    cd(paths.stim)
    if stimvec(2,trl) ~= 0
        stimname = stimvec_fname{trl};
        loadpict(stimname,1);
        
        categories(1,i) = stimvec(2,trl);
        
        sc=sc+1; SOT_s(sc,1)=toc; tic % reset stopwatch
        drawpict(1);
        if eyeflag==1
            Eyelink('Message',[num2str(trl) 's'])  % send trialnumber and stimulustype to eyetracker
        end
        wait(config.timing.stim);
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = stimname;
        
        % ----------------------------------------------------------------- %
        
        
        % ----------------------- response screen ------------------------- %
        
        tmp = toc; % clearkeys; readkeys;
        stimonset1 = drawpict(2);
        if eyeflag==1
            Eyelink('Message',[num2str(trl) 'r'])  % send trialnumber and stimulustype to eyetracker
        end
        waitkeydown(config.timing.response,[InKey,OutKey]);
        %     wait(config.timing.response);
        
        % RT
        [keypress1, t1, n1] = getkeydown;
        if n1 == 0; %if they didn't press anything
            response1 = NaN; % mark their response as nan
            keypress1 = NaN;
            rt1 = NaN; % mark their reaction time as nan
            wait(0);
        else % otherwise
            drawpict(3);
            response1 = keypress1(1); % their response is whatever key they pressed.
            rt1 = t1(1) - stimonset1; % and their reaction time is the time they pressed the button-the time the stimulus apprered
            rc1=rc1+1; SOT_r(rc1,1) = tmp+rt1; clear tmp
            wait(2000-rt1)
        end
        
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = 'response';
        
        % accuracy
        if (response1 == InKey && strcmp(stimname(5),'i'))
            accuracy = 1;
            fprintf('\nCorrect \n')
        elseif (response1 == OutKey && strcmp(stimname(5),'o'))
            accuracy = 1;
            fprintf('\nCorrect \n')
        elseif (response1 == OutKey  && strcmp(stimname(5),'i') && ~isnan(response1)) ...
                || (response1 == InKey && strcmp(stimname(5),'i') && ~isnan(response1))
            accuracy = 0;
            fprintf('\nIncorrect \n')
        elseif isnan(response1)
            accuracy = 0;
            fprintf('\nMissed \n')
        end
        
        
        % ----------------------------------------------------------------- %
        
        
        % ----------------------- valence ratings ------------------------- %
        
        clearpict(1);clearpict(2);
        loadpict('checker.bmp',1,0,0,bg_w,bg_h);loadpict('checker.bmp',2,0,0,bg_w,bg_h);
        loadpict('valence.bmp',1,0,0);
        loadpict('valence_received.bmp',2,0,0);
        
        tmp = toc; % clearkeys; readkeys;
        stimonset2 = drawpict(1);
        if eyeflag==1
            Eyelink('Message',[num2str(trl) 'v'])  % send trialnumber and stimulustype to eyetracker
        end
        waitkeydown(config.timing.response,[VeryKey,NVeryKey]);
        %     wait(config.timing.response);
        
        % RT
        [keypress2, t2, n2] = getkeydown;
        if n2 == 0 % if they didn't press anything
            response2 = NaN; % mark their response as nan
            keypress2 = NaN;
            rt2 = NaN; % mark their reaction time as nan
            wait(0);
            valenceR = NaN;
        else % otherwise
            drawpict(2);
            response2 = keypress2(1); % their response is whatever key they pressed.
            rt2 = t2(1) - stimonset2; % and their reaction time is the time they pressed the button-the time the stimulus apprered
            rc2=rc2+1; SOT_v(rc2,1) = tmp+rt2; clear tmp
            wait(2000-rt2)
            % record valence
            if keypress2 == VeryKey
                valenceR = 1;
            elseif keypress2 == NVeryKey
                valenceR = 0;
            end
        end
        
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = 'valence';
        

        
        % ----------------------------------------------------------------- %
        
        clearpict(1);clearpict(2); clearpict(3);
        
    elseif stimvec(2,trl) == 0 % if the stim type is null
        stimname = 'null';
        setforecolour(.5,.5,.5);
        settextstyle('Arial', 300);
        preparestring('+', 1);
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = 'Null';
        categories(1,i) = 0;
        response1 = NaN; % mark their response as nan
        keypress1 = NaN;
        rt1 = NaN; % mark their reaction time as nan
        response2 = NaN; % mark their response as nan
        keypress2 = NaN;
        rt2 = NaN; % mark their reaction time as nan
        valenceR = NaN;
    end
    
    % ------------------------ present fix x -------------------------- %
    
    cd(paths.parent); loadpict('checker.bmp',1,0,0,bg_w,bg_h);
    setforecolour(.5,.5,.5);
    settextstyle('Arial', 300);
    preparestring('+', 1);  % Prepare buffer 1 (fixation cross)
    
    f1c=f1c+1; SOT_f1(f1c,1)=toc; tic  % reset stopwatch
    drawpict(1)             % load fixation cross
    if eyeflag==1
        Eyelink('Message',[num2str(trl) 'f1'])  % send trialnumber and stimulustype to eyetracker
    end
    wait(config.timing.fixation(trl))
    eventmarker = eventmarker+1;
    results.presentation{eventmarker,1} = '+';
    
    clearpict(1);
    
    % ----------------------------------------------------------------- %
    
    
    % record
    results.keypress(i,1) = keypress1; % response
    results.keypress(i,2) = keypress2; % valence
    results.rt_resp(i,1)  = rt1;
    results.rt_valnce(i,1)= rt2;
    results.accuracy(i,1) = accuracy;
    results.trl{i,1}      = stimname;
    results.trl{i,2}      = categories(1,i);
    results.trl{i,3}      = valenceR;
    
    % calculate cumulative SOT
    results.SOT.cumulative.fix  = cumsum(SOT_f1);
    results.SOT.cumulative.stim = cumsum(SOT_s);
    results.SOT.cumulative.resp = cumsum(SOT_r);
    results.SOT.cumulative.val  = cumsum(SOT_v);
    results.SOT.raw.fix  = SOT_f1;
    results.SOT.raw.stim = SOT_s;
    results.SOT.raw.resp = SOT_r;
    results.SOT.raw.val  = SOT_v;
    
    % intermediate save
    dat.emomem.task.results = results;
    save([paths.behav num2str(ID) '_EM_backup.mat'],'dat')
    
    clear stimname
    
end

%% wrap up
% calculate final cumulative SOT
results.SOT.cumulative.fix  = cumsum(SOT_f1);
results.SOT.cumulative.stim = cumsum(SOT_s);
results.SOT.cumulative.resp = cumsum(SOT_r);
results.SOT.cumulative.val  = cumsum(SOT_v);
results.SOT.raw.fix  = SOT_f1;
results.SOT.raw.stim = SOT_s;
results.SOT.raw.resp = SOT_r;
results.SOT.raw.val  = SOT_v;

% save results
dat.emomem.task.results = results;
dat.emomem.task.labels  = {'contingencies-> 11=in_emo, 12=in_neu, 21=out_emo, 22=out_neu'; ...
    'results.trl: filenames, contingencies, valenceRating(1=very,0=notvery)'};
dat.emomem.task.config.keymap = map;
save([paths.behav num2str(ID) '_EM.mat'],'dat')

% terminate protocol
clearpict(1); cd(paths.supports)
loadpict('ending.png',1,0,0)
drawpict(1); wait(3000);
stop_cogent;

end