%% EMOTIONAL MEMORY TASK : Scene classification task with emotional/neutral images
%  In emergency contact Alex:
%              +49 162 713 20 49
%              yeo-jin.yi.15@ucl.ac.uk

%%  work log

%   23_01_2019 created the script

%% START

function [dat] = runEmoMemTask_cg(ID,AgeGroup,EyeFlag,MRIFlag,eyetracker_filename,Schedule,path_ECHO)

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
paths.parent   = path_ECHO;
paths.supports = [paths.parent 'stim\supports\'];
paths.stim     = [paths.parent 'stim\old\'];
paths.pupil    = [paths.parent 'data\pupil\'];
paths.behav    = [paths.parent 'data\behav\'];
% paths.cogent   = 'E:\Dorothea\Cogent2000v1.33';
paths.cogent   = [paths.parent 'Cogent2000v1.33\'];
addpath(genpath(paths.cogent))

% create the data structure
cd(paths.behav)
if exist([num2str(ID) '_EM.mat']) == 2
    load([paths.behav num2str(ID) '_EM.mat'])
else
    dat = [];
    dat.ID       = ID;
    dat.AgeGroup = AgeGroup;
    dat.emomem = [];
    dat.emomem.task.ETinfo  = fname_eyetracker;
end
if Schedule == 1
    load([paths.parent 'scripts\EMtask_GA_optimised_ver1.mat'])
else
    load([paths.parent 'scripts\EMtask_GA_optimised_ver2.mat'])
end

% scanner preparation
if MRIFlag == 1
    dat.emomem.task.MRIinfo   = [];
    dat.emomem.task.MRIinfo.scanPort = 1;
    dat.emomem.task.MRIinfo.dummy    = 5;        % no. of dummy vols
    dat.emomem.task.MRIinfo.nslice   = 51;       % no. of slices
    dat.emomem.task.MRIinfo.TE       = 32;       % time for each slice in msec
    dat.emomem.task.MRIinfo.TR       = 3.6; % MRIinfo.nslice * MRIinfo.TE/1000; % time for volume in msec
%     dat.emomem.task.MRIinfo.nvols    = ceil(MRIinfo.dummy + 5 + ...  % no. of vols per run (5 dummy + 400 + 5 post)
%         (options.scanblocklength * options.meanTrialLength) / MRIinfo.TR);
%     dat.emomem.task.MRIinfo.total_slices = MRIinfo.nvols*MRIinfo.nslice;    % per run!
    %     [MRIinfo] = runScannerSetup_MRIpilot; % redundant
end

%% experimental setups

% condition markers
in_emo = 1; in_neu = 3;
out_emo = 2; out_neu = 4;

% timings
% fixjit = repmat(500:50:700,1,24);
config.timing.fixation = [(design_struct.eventlist(:,4)+design_struct.eventlist(:,5))']*1000;
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
config.numtrials.null    = 35;

% contingency vector
trlvec          = 1:config.numtrials.total;
contingencies   = design_struct.eventlist(:,3)';
indoors         = contingencies(find(contingencies==1 | contingencies==3));
outdoors        = contingencies(find(contingencies==2 | contingencies==4));
nulls           = contingencies(find(contingencies==0));
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
config_display(1, 5, [0 0 0], [1 1 1], 'Arial', 20, 10, 0); % 1024x768
% config_display(screenMode, screenRes, bg ,fontcol, fontName, fontSize, number_of_buffers,0);   % open graphics window
config_sound;
config_keyboard;

% scanner input
if MRIFlag == 1
    config_serial(dat.emomem.task.MRIinfo.scanPort); % links to scanner for waitslice
end

monitor = get(0,'screensize');
bg_w = monitor(3); bg_h = monitor(4);

%%%%%%%%%%%%
start_cogent
%%%%%%%%%%%%

% cgopen(4,0,0,2) % at the MRI, choose the secondary screen
cgopen(5,0,0,1)

map     = getkeymap; % define keyboard IDs
InKey   = 28; %%%%% check again
OutKey  = 29;
VeryKey = 29;
NVeryKey= 28;

% print instructions
cd(paths.supports)
cgloadbmp(1,'emomem_inst_task_orig.bmp');
cgsetsprite(0);
cgdrawsprite(1,0,0);
t0_inst=cgflip(0,0,0); waitkeydown(inf,[InKey,OutKey]); clearkeys

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
    cd(paths.supports)
    cgloadbmp(1,'operator_hail.bmp');
    cgsetsprite(0);
    cgdrawsprite(1,0,0);
    cgflip(0,0,0);
    waitkeydown(inf,map.Space)
    
    % start the dummy scan
    cd(paths.supports)
    cgloadbmp(1,'ready.bmp');
    cgsetsprite(0);
    cgdrawsprite(1,0,0);
    t0_standby=cgflip(0,0,0);
    scannerinput = map.T;
    
    %%%%%%%%%%%%%%
    tic % ding ding ding
    %%%%%%%%%%%%%%
    
    [k, dat.emomem.task.MRIinfo.dummyscanfirst] = waitkeydown(inf,scannerinput); % the first dummy scan
    for dum = 1:(dat.emomem.task.MRIinfo.dummy - 2)
        waitkeydown(inf,scannerinput);
    end
    [k, dat.emomem.task.MRIinfo.dummyscanlast] = waitkeydown(inf,scannerinput); % the last dummy scan
    
%     %%%%%%%%%%%%%%
%     tic % ding ding ding
%     %%%%%%%%%%%%%%    
    
    cd(paths.supports); cgloadbmp(1,'checker_b.bmp',1280,1024);
    cgsetsprite(0);
    cgdrawsprite(1,0,0);
    cgfont('Arial',300);
    cgpencol(.5,.5,.5)
    cgtext('+',0,0);  % Prepare buffer 1 (fixation cross)
    t0_fix0 = cgflip(0,0,0);   t0_fix0_raw = toc;    % load fixation cross
%     if eyeflag==1
%         Eyelink('Message',[num2str(0) 'f0'])  % send trialnumber and stimulustype to eyetracker
%     end
    wait(dat.emomem.task.MRIinfo.TR*1000)
    clearpict(1);
    
else
    
    %%%%%%%%%%%%%%
    tic % ding ding ding
    %%%%%%%%%%%%%%
    
    cd(paths.supports); cgloadbmp(1,'checker_b.bmp',1280,1024);
    cgsetsprite(0);
    cgdrawsprite(1,0,0);
    cgfont('Arial',300);
    cgpencol(.5,.5,.5)
    cgtext('+',0,0);  % Prepare buffer 1 (fixation cross)
    cgflip(0,0,0);             % load fixation cross
%     if eyeflag==1
%         Eyelink('Message',[num2str(0) 'f0'])  % send trialnumber and stimulustype to eyetracker
%     end
    wait(1000)
    
end

results = []; results.SOT = [];
results.SOT.cgflip.t0_fix0 = t0_fix0;
results.SOT.cgflip.t0_standby = t0_standby;
results.SOT.cgflip.t0_inst = t0_inst;
results.SOT.raw.t0_fix0 = t0_fix0_raw;
eventmarker = 0;
i = 0; f1c = 0; sc = 0; rc1 = 0; rc2 = 0; vc = 0; nc = 0;
SOT_f1=[]; SOT_s=[]; SOT_r=[]; SOT_v = []; SOT_n = [];
SOT_cg_s=[];SOT_cg_f1=[];SOT_cg_r=[]; SOT_cg_v=[]; SOT_cg_n = [];
for trl = 1:size(stimvec,2)
        
    % ------------------------ present image -------------------------- %
    i = i+1;
    fprintf('\nTrial %d\n', trl)
    
    clearkeys;
    
    %%%%%% load background %%%%%%
    cd(paths.supports); cgloadbmp(1,'checker_b.bmp',1280,1024);
    
    
    if stimvec(2,trl) ~= 0
        %%%%%% load stimulus %%%%%%
        cd(paths.stim)
        stimname = stimvec_fname{trl};
        cgloadbmp(2,stimname,500,360);
        categories(1,i) = stimvec(2,trl);
        
        %%%%%% draw stimulus screen %%%%%%
        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgdrawsprite(2,0,0);
        sc=sc+1; SOT_s(sc,1)=toc;
        SOT_cg_s(sc) = cgflip(0,0,0);
        if eyeflag==1
            Eyelink('Message',[num2str(trl) 's'])  % send trialnumber and stimulustype to eyetracker
        end
        wait(config.timing.stim);
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = stimname;
        
        % ----------------------------------------------------------------- %
        
        
        % --------------------- IN/OUT response screen -------------------- %
        
        %%%%%% load IN/OUT Question screen %%%%%%
        cd(paths.supports);
        cgloadbmp(1,'checker_b.bmp',1280,1024);
        cgloadbmp(2,'inout_orig.bmp');
        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgdrawsprite(2,0,0); tmp = toc;
        responset1 = cgflip(0,0,0);
        if eyeflag==1
            Eyelink('Message',[num2str(i) 'r'])  % send trialnumber and stimulustype to eyetracker
        end
        waitkeydown(config.timing.response,[InKey,OutKey]);
        
        % RT
        [keypress1, t1, n1] = getkeydown;
        if n1 == 0; %if they didn't press anything
            response1 = NaN; % mark their response as nan
            keypress1 = NaN;
            rt1 = NaN; % mark their reaction time as nan
            wait(0);
        else % otherwise
            %%%%%% load answer screen %%%%%%
            cgloadbmp(1,'checker_b.bmp',1280,1024);
            cgloadbmp(2,'inout_received_orig.bmp');
            cgsetsprite(0);
            cgdrawsprite(1,0,0); cgdrawsprite(2,0,0); SOT_cg_r(i) = cgflip(0,0,0);
            response1 = keypress1(1); % their response is whatever key they pressed.
            rt1 = t1(1)/1000 - responset1; rt1 = rt1*1000; % and their reaction time is the time they pressed the button-the time the stimulus apprered
            rc1=rc1+1; SOT_r(rc1,1) = tmp+rt1; SOT_cg_r(rc1,1) = tmp+rt1; clear tmp
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
        else
            accuracy = NaN;
        end
        
        
        % ----------------------------------------------------------------- %
        
        
        % ----------------------- valence ratings ------------------------- %
        
        %%%%%% load VALENCE Question screen %%%%%%
        cgloadbmp(1,'checker_b.bmp',1280,1024);
        cgloadbmp(2,'valence_orig.bmp');
        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgdrawsprite(2,0,0); tmp = toc;
        responset2 = cgflip(0,0,0);
        if eyeflag==1
            Eyelink('Message',[num2str(trl) 'v'])  % send trialnumber and stimulustype to eyetracker
        end
        waitkeydown(config.timing.response,[VeryKey,NVeryKey]);
        
        % RT
        [keypress2, t2, n2] = getkeydown;
        if n2 == 0 % if they didn't press anything
            response2 = NaN; % mark their response as nan
            keypress2 = NaN;
            rt2 = NaN; % mark their reaction time as nan
            wait(0);
            valenceR = NaN;
        else % otherwise
            %%%%%% load answer screen %%%%%%
            cgloadbmp(1,'checker_b.bmp',1280,1024);
            cgloadbmp(2,'valence_received_orig.bmp');
            cgsetsprite(0);
            cgdrawsprite(1,0,0); cgdrawsprite(2,0,0); SOT_cg_v(i) = cgflip(0,0,0);
            response2 = keypress2(1); % their response is whatever key they pressed.
            rt2 = t2(1)/1000 - responset2; rt2 = rt2*1000; % and their reaction time is the time they pressed the button-the time the stimulus apprered
            rc2=rc2+1; SOT_v(rc2,1) = tmp+rt2; SOT_cg_v(rc2,1) = tmp+rt2; clear tmp
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
        
        
    elseif stimvec(2,trl) == 0 % if the stim type is null
        stimname = 'null';
        cd(paths.supports);
        cgloadbmp(1,'checker_b.bmp',1280,1024);
        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgfont('Arial',300)
        cgpencol(.5,.5,.5)
        cgtext('+',0,0);
        
        nc = nc+1; SOT_n(nc,1)=toc;
        SOT_cg_n(nc) = cgflip(0,0,0);
        wait(config.timing.stim)
        
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
    
    cd(paths.supports); cgloadbmp(1,'checker_b.bmp',1280,1024);
    cgsetsprite(0);
    cgdrawsprite(1,0,0);
    cgfont('Arial',300)
    cgpencol(.5,.5,.5)
    cgtext('+',0,0);
    f1c=f1c+1; SOT_f1(f1c,1)=toc;
    SOT_cg_f1(f1c) = cgflip(0,0,0);             % load fixation cross
    if eyeflag==1
        Eyelink('Message',[num2str(i) 'f'])  % send trialnumber and stimulustype to eyetracker
    end
    wait(config.timing.fixation(trl))
    eventmarker = eventmarker+1;
    results.presentation{eventmarker,1} = '+';
        
    % ----------------------------------------------------------------- %
    
    
    if stimvec(2,trl) == 0 % if the stim type is null
        results.keypress(i,1) = NaN; % response
        results.keypress(i,2) = NaN; % valence
        results.rt_resp(i,1)  = NaN;
        results.rt_valnce(i,1)= NaN;
        results.accuracy(i,1) = NaN;
        results.trl{i,1}      = 'null';
        results.trl{i,2}      = 0;
        results.trl{i,3}      = NaN;
    else
        % record
        
        try
            results.keypress(i,1) = keypress1; % response
            results.keypress(i,2) = keypress2; % valence
        catch
            try
                results.keypress(i,1) = keypress1(1); % response
                results.keypress(i,2) = keypress2(1); % valence
            catch
                results.keypress(i,1) = NaN; % response
                results.keypress(i,2) = NaN; % valence
            end
        end
        
        results.rt_resp(i,1)  = rt1;
        results.rt_valnce(i,1)= rt2;
        results.accuracy(i,1) = accuracy;
        results.trl{i,1}      = stimname;
        results.trl{i,2}      = categories(1,i);
        results.trl{i,3}      = valenceR;
        
        % calculate cumulative SOT
        results.SOT.raw.fix  = SOT_f1;
        results.SOT.raw.stim = SOT_s;
        results.SOT.raw.resp = SOT_r;
        results.SOT.raw.val  = SOT_v;
        results.SOT.raw.nullt= SOT_n;
        results.SOT.cgflip.stim = SOT_cg_s;
        results.SOT.cgflip.resp = SOT_cg_r;
        results.SOT.cgflip.fix  = SOT_cg_f1;
        results.SOT.cgflip.val  = SOT_cg_v;
        results.SOT.cgflip.nullt= SOT_cg_n;
    end
    
    % intermediate save
    dat.emomem.task.results = results;
    save([paths.behav num2str(ID) '_EM_backup.mat'],'dat')
    
    clear stimname
    
end

%% wrap up
% calculate final cumulative SOT
results.SOT.raw.fix  = SOT_f1;
results.SOT.raw.stim = SOT_s;
results.SOT.raw.resp = SOT_r;
results.SOT.raw.val  = SOT_v;
results.SOT.raw.nullt= SOT_n;
results.SOT.cgflip.stim = SOT_cg_s;
results.SOT.cgflip.resp = SOT_cg_r;
results.SOT.cgflip.fix  = SOT_cg_f1;
results.SOT.cgflip.val  = SOT_cg_v;
results.SOT.cgflip.nullt= SOT_cg_n;

% save results
dat.emomem.task.results = results;
dat.emomem.task.labels  = {'contingencies-> 11=in_emo, 12=in_neu, 21=out_emo, 22=out_neu'; ...
    'results.trl: filenames, contingencies, valenceRating(1=very,0=notvery)'};
dat.emomem.task.config.keymap = map;
save([paths.behav num2str(ID) '_EM.mat'],'dat')

% terminate protocol
cd(paths.supports)
cgloadbmp(1,'instruction_ende.bmp');
cgsetsprite(0);
cgdrawsprite(1,0,0); 
cgflip(0,0,0); wait(3000);
cgshut;

end