%% EMOTIONAL MEMORY TEST : Scene classification task with emotional/neutral images
%  In emergency contact Alex:
%              +49 162 713 20 49
%              yeo-jin.yi.15@ucl.ac.uk

%%  work log

%   23_01_2019 created the script

%% START
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
function [dat] = runEmoMemTest_cg(ID,AgeGroup,EyeFlag,MRIFlag,eyetracker_filename,Schedule,path_ECHO)

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
% paths.parent   = '/Users/yeojin/Desktop/B_scripts/BD_tasks/task_scripts/EmoMem/';
paths.supports = [paths.parent 'stim\supports\'];
paths.stim_old = [paths.parent 'stim\old\'];
paths.stim_new = [paths.parent 'stim\new\'];
paths.pupil    = [paths.parent 'data\pupil\'];
paths.behav    = [paths.parent 'data\behav\'];
paths.cogent   = [paths.parent 'Cogent2000v1.33\'];
addpath(genpath(paths.cogent))

% create the data structure
cd(paths.behav)
if exist([num2str(ID) '_EM.mat']) == 2
    load([paths.behav num2str(ID) '_EM.mat'])
    dat.emomem.test = [];
else
    error('Previous task results do not exist!!!')
end
if Schedule == 1
    load([paths.parent 'scripts\EMtest_GA_optimised_ver1.mat'])
else
    load([paths.parent 'scripts\EMtest_GA_optimised_ver2.mat'])
end

% scanner preparation
if MRIFlag == 1
    dat.emomem.test.MRIinfo   = [];
    dat.emomem.test.MRIinfo.scanPort = 1;
    dat.emomem.test.MRIinfo.dummy    = 5;        % no. of dummy vols
    dat.emomem.test.MRIinfo.nslice   = 51;       % no. of slices
    dat.emomem.test.MRIinfo.TE       = 32;       % time for each slice in msec
    dat.emomem.test.MRIinfo.TR       = 3.6; % MRIinfo.nslice * MRIinfo.TE/1000; % time for volume in msec
%     dat.emomem.test.MRIinfo.nvols    = ceil(MRIinfo.dummy + 5 + ...  % no. of vols per run (5 dummy + 400 + 5 post)
%         (options.scanblocklength * options.meanTrialLength) / MRIinfo.TR);
%     dat.emomem.test.MRIinfo.total_slices = MRIinfo.nvols*MRIinfo.nslice;    % per run!
    %     [MRIinfo] = runScannerSetup_MRIpilot; % redundant
end

%% experimental setups

config = [];

% number of trials
config.numtrials.total      = 160;

config.numtrials.old_total  = 80;
config.numtrials.old_in_emo = 20;
config.numtrials.old_out_emo= 20;
config.numtrials.old_in_neu = 20;
config.numtrials.old_out_neu= 20;

config.numtrials.new_total  = 80;
config.numtrials.new_in_emo = 20;
config.numtrials.new_out_emo= 20;
config.numtrials.new_in_neu = 20;
config.numtrials.new_out_neu= 20;

% condition markers
in_emo_old = 1; in_emo_new = 2; 
in_neu_old = 5; in_neu_new = 6;
out_emo_old = 3; out_emo_new = 4;
out_neu_old = 7; out_neu_new = 8;

% contingency vector
trlvec          = 1:config.numtrials.total;
contingencies   = design_struct.eventlist(:,3)';
olds            = contingencies(find(contingencies==1 | contingencies==3 | contingencies==5 | contingencies==7));
news            = contingencies(find(contingencies==2 | contingencies==4 | contingencies==6 | contingencies==8));

% timing
config.timing.outcome  = 2500;
config.timing.fixation = [(design_struct.eventlist(:,4)+design_struct.eventlist(:,5))']*1000;
config.timing.cue      = 2500;
config.timing.response = 2000;

% --- for breezeing through --- %
% config.timing.fixation= zeros(1,240);
% config.timing.fixation2= zeros(1,240);
% config.timing.cue      = 1;
% config.timing.outcome  = 1;
% ----------------------------- %


%% roster

% ------------------- filenames setup: OLD stim list ------------------- %

% indoor-emotional
cd(paths.stim_old); cd emotional; cd in;
tmp = dir('*.bmp');
[fname.old_in_emo{1:config.numtrials.old_in_emo}] = deal(tmp(1:config.numtrials.old_in_emo).name); clear tmp;
% outdoor-emotional
cd(paths.stim_old); cd emotional; cd out;
tmp = dir('*.bmp');
[fname.old_out_emo{1:config.numtrials.old_out_emo}] = deal(tmp(1:config.numtrials.old_out_emo).name); clear tmp;
% indoor-neutral
cd(paths.stim_old); cd neutral; cd in;
tmp = dir('*.bmp');
[fname.old_in_neu{1:config.numtrials.old_in_neu}] = deal(tmp(1:config.numtrials.old_in_neu).name); clear tmp;
% outdoor-neutral
cd(paths.stim_old); cd neutral; cd out;
tmp = dir('*.bmp');
[fname.old_out_neu{1:config.numtrials.old_out_neu}] = deal(tmp(1:config.numtrials.old_out_neu).name); clear tmp;



% ------------------- filenames setup: NEW stim list ------------------- %

% indoor-emotional
cd(paths.stim_new); cd emotional; cd in;
tmp = dir('*.bmp');
[fname.new_in_emo{1:config.numtrials.new_in_emo}] = deal(tmp(1:config.numtrials.new_in_emo).name); clear tmp;
% outdoor-emotional
cd(paths.stim_new); cd emotional; cd out;
tmp = dir('*.bmp');
[fname.new_out_emo{1:config.numtrials.new_out_emo}] = deal(tmp(1:config.numtrials.new_out_emo).name); clear tmp;
% indoor-neutral
cd(paths.stim_new); cd neutral; cd in;
tmp = dir('*.bmp');
[fname.new_in_neu{1:config.numtrials.new_in_neu}] = deal(tmp(1:config.numtrials.new_in_neu).name); clear tmp;
% outdoor-neutral
cd(paths.stim_new); cd neutral; cd out;
tmp = dir('*.bmp');
[fname.new_out_neu{1:config.numtrials.new_out_neu}] = deal(tmp(1:config.numtrials.new_out_neu).name); clear tmp;

% schedule presentation
tmpvec       = trlvec;
stimvec_temp = [];
inemoold = 0; outemoold = 0; inneuold = 0; outneuold = 0;
inemonew = 0; outemonew = 0; inneunew = 0; outneunew = 0;
for sc = 1:length(contingencies)
    if contingencies(sc) == in_emo_old
        inemoold = inemoold+1;
        stimvec_temp(1,sc) = tmpvec(sc);
        stimvec_temp(2,sc) = contingencies(sc);
        stimvec_fname{1,sc}= fname.old_in_emo{inemoold};
    elseif contingencies(sc) == in_emo_new
        inemonew = inemonew+1;
        stimvec_temp(1,sc) = tmpvec(sc);
        stimvec_temp(2,sc) = contingencies(sc);
        stimvec_fname{1,sc}= fname.new_in_emo{inemonew};
    elseif contingencies(sc) == out_emo_old
        outemoold = outemoold+1;
        stimvec_temp(1,sc) = tmpvec(sc);
        stimvec_temp(2,sc) = contingencies(sc);
        stimvec_fname{1,sc}= fname.old_out_emo{outemoold};
    elseif contingencies(sc) == out_emo_new
        outemonew = outemonew+1;
        stimvec_temp(1,sc) = tmpvec(sc);
        stimvec_temp(2,sc) = contingencies(sc);
        stimvec_fname{1,sc}= fname.new_out_emo{outemonew};
    elseif contingencies(sc) == in_neu_old
        inneuold = inneuold+1;
        stimvec_temp(1,sc) = tmpvec(sc);
        stimvec_temp(2,sc) = contingencies(sc);
        stimvec_fname{1,sc}= fname.old_in_neu{inneuold};
    elseif contingencies(sc) == in_neu_new
        inneunew = inneunew+1;
        stimvec_temp(1,sc) = tmpvec(sc);
        stimvec_temp(2,sc) = contingencies(sc);
        stimvec_fname{1,sc}= fname.new_in_neu{inneunew};
    elseif contingencies(sc) == out_neu_old
        outneuold = outneuold+1;
        stimvec_temp(1,sc) = tmpvec(sc);
        stimvec_temp(2,sc) = contingencies(sc);
        stimvec_fname{1,sc}= fname.old_out_neu{outneuold};
    elseif contingencies(sc) == out_neu_new
        outneunew = outneunew+1;
        stimvec_temp(1,sc) = tmpvec(sc);
        stimvec_temp(2,sc) = contingencies(sc);
        stimvec_fname{1,sc}= fname.new_out_neu{outneunew};
    end
    
end
stimvec = stimvec_temp; clear stimvec_temp
config.stim.stimvec  = stimvec;
config.stim.stimlist = stimvec_fname;
config.stim.filelist = fname;

% recall phase stimuli list
% config.stim.new = vertcat( horzcat(fname.new_in_emo,fname.new_in_neu,fname.new_out_emo,fname.new_out_neu),...
%     num2cell( ones(1,config.numtrials.new_total)*2 ),...
%     num2cell( horzcat( ones(1,config.numtrials.new_in_emo)*in_emo,ones(1,config.numtrials.new_in_neu)*in_neu,...
%     ones(1,config.numtrials.new_out_emo)*out_emo,ones(1,config.numtrials.new_out_neu)*out_neu ) ));
% config.stim.old = vertcat( horzcat(fname.old_in_emo,fname.old_in_neu,fname.old_out_emo,fname.old_out_neu),...
%     num2cell( ones(1,config.numtrials.old_total) ),...
%     num2cell( horzcat( ones(1,config.numtrials.old_in_emo)*in_emo,ones(1,config.numtrials.old_in_neu)*in_neu,...
%     ones(1,config.numtrials.old_out_emo)*out_emo,ones(1,config.numtrials.old_out_neu)*out_neu ) ));
% config.stim.all = horzcat(config.stim.new,config.stim.old);


% scramble
% stimvec = config.stim.all;
% stimvec = stimvec(:,randperm(size(stimvec,2)));
% 
% config.stim.stimvec = stimvec;
% dat.emomem.test.config = config;


%% run

% configs
config_display(1, 5, [0 0 0], [1 1 1], 'Arial', 35, 4); % 1280*1024
config_sound;
config_keyboard;

if MRIFlag == 1
    config_serial(dat.emomem.test.MRIinfo.scanPort); % links to scanner for waitslice
end

monitor = get(0,'screensize');
bg_w = monitor(3); bg_h = monitor(4);

%%%%%%%%%%%%
start_cogent
%%%%%%%%%%%%

% cgopen(4,0,0,2) % at the MRI, choose the secondary screen
cgopen(5,0,0,1)

map     = getkeymap; % define keyboard IDs
OldKey  = 28;
NewKey  = 29;
SureKey = 28;
GuessKey= 29;

% print instructions
cd(paths.supports)
cgloadbmp(1,'memorytest_instruction1_orig.bmp');
cgsetsprite(0);
cgdrawsprite(1,0,0);
cgflip(0,0,0); waitkeydown(inf,[OldKey,NewKey]); clearkeys

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
    cgflip(0,0,0);
    scannerinput = 33;
    [k, dat.emomem.test.MRIinfo.dummyscanfirst] = waitkeydown(inf,scannerinput); % the first dummy scan
    for dum = 1:(dat.emomem.test.MRIinfo.dummy - 2)
        waitkeydown(inf,scannerinput);
    end
    [k, dat.emomem.test.MRIinfo.dummyscanlast] = waitkeydown(inf,scannerinput); % the last dummy scan
    
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
    wait(dat.emomem.test.MRIinfo.TR*1000)
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


clear trl c1 c2
c1 = 0; c2 = 0; results = []; results.SOT = [];
eventmarker = 0;
i = 0; f1c = 0; sc = 0; rc = 0; cc = 0; vc = 0; nc = 0;
SOT_f1=[]; SOT_s=[]; SOT_r=[]; SOT_c = [];
SOT_cg_s=[];SOT_cg_f1=[]; SOT_cg_r=[]; SOT_cg_c=[];
resp = []; accu = []; confi = [];
for trl = 1:config.numtrials.total
        
    % ------------------------ present image -------------------------- %
    i = i+1;
    fprintf('\nTrial %d\n', trl)
    
    clearkeys;
    
    %%%%%% load background %%%%%%
    cd(paths.supports); cgloadbmp(1,'checker_b.bmp',1280,1024);
    
    %%%%%% load stimulus %%%%%%
    if stimvec(2,trl) == 1 || stimvec(2,trl) == 3 || stimvec(2,trl) == 5 || stimvec(2,trl) == 7
        cd(paths.stim_old)
    else
        cd(paths.stim_new)
    end
    cgloadbmp(2,stimvec_fname{1,trl},500,360);
    categories(1,i) = stimvec(2,trl);
    
    %%%%%% draw stimulus screen %%%%%%
    sc=sc+1; SOT_s(sc,1)=toc;
    cgsetsprite(0);
    cgdrawsprite(1,0,0);
    cgdrawsprite(2,0,0);
    SOT_cg_s(sc) = cgflip(0,0,0);
    if eyeflag==1
        Eyelink('Message',[num2str(trl) 's'])  % send trialnumber and stimulustype to eyetracker
    end
    wait(config.timing.cue);
    eventmarker = eventmarker+1;
    results.presentation{eventmarker,1} = stimvec_fname{1,trl};
    % ----------------------------------------------------------------- %

    
    % --------------------- OLD/NEW response screen -------------------- %
        
    %%%%%% load IN/OUT Question screen %%%%%%
    cd(paths.supports);
    cgloadbmp(1,'checker_b.bmp',1280,1024);
    cgloadbmp(2,'newold_orig.bmp');
    cgsetsprite(0);
    cgdrawsprite(1,0,0);
    cgdrawsprite(2,0,0); tmp = toc;
    responset_ON = cgflip(0,0,0);
    if eyeflag==1
        Eyelink('Message',[num2str(i) 'r'])  % send trialnumber and stimulustype to eyetracker
    end
    waitkeydown(config.timing.response,[NewKey,OldKey]);
    
    % RT
    [key_ON, t_ON, n_ON] = getkeydown;
    if n_ON == 0; %if they didn't press anything
        response1 = NaN; % mark their response as nan
        key_ON = NaN;
        rt_ON = NaN; % mark their reaction time as nan
        wait(0);
    else % otherwise
        %%%%%% load answer screen %%%%%%
        cgloadbmp(1,'checker_b.bmp',1280,1024);
        cgloadbmp(2,'newold_received_orig.bmp');
        cgsetsprite(0);
        cgdrawsprite(1,0,0); cgdrawsprite(2,0,0); SOT_cg_r(i) = cgflip(0,0,0);
        response1 = key_ON(1); % their response is whatever key they pressed.
        rt_ON = t_ON(1)/1000 - responset_ON; rt_ON = rt_ON*1000; % and their reaction time is the time they pressed the button-the time the stimulus apprered
        rc=rc+1; SOT_r(rc,1) = tmp+rt_ON; SOT_cg_r(rc,1) = tmp+rt_ON; clear tmp
        wait(2000-rt_ON)
    end
    
    eventmarker = eventmarker+1;
    results.presentation{eventmarker,1} = 'response_ON';
    
    % accuracy
    if key_ON == OldKey
        resp = 1;
        if stimvec(2,trl) == 1 || stimvec(2,trl) == 3 || stimvec(2,trl) == 5 || stimvec(2,trl) == 7
            accu = 1;
        else
            accu = 0;
        end
    elseif key_ON == NewKey
        resp = 2;
        if stimvec(2,trl) == 2 || stimvec(2,trl) == 4 || stimvec(2,trl) == 6 || stimvec(2,trl) == 8
            accu = 1;
        else
            accu = 0;
        end
    end
    % ----------------------------------------------------------------- %
    
    
    % --------------------- confidence ratings ------------------------ %
    
    %%%%%% load confidence Question screen %%%%%%
    cd(paths.supports);
    cgloadbmp(1,'checker_b.bmp',1280,1024);
    cgloadbmp(2,'memorytest_confidence_orig.bmp');
    cgsetsprite(0);
    cgdrawsprite(1,0,0);
    cgdrawsprite(2,0,0); tmp = toc;
    responset_C = cgflip(0,0,0);
    if eyeflag==1
        Eyelink('Message',[num2str(i) 'c'])  % send trialnumber and stimulustype to eyetracker
    end
    waitkeydown(config.timing.response,[SureKey,GuessKey]);
    
    % RT
    [key_SG, t_SG, n_SG] = getkeydown;
    if n_SG == 0; %if they didn't press anything
        response1 = NaN; % mark their response as nan
        key_SG = NaN;
        rt_SG = NaN; % mark their reaction time as nan
        wait(0);
    else % otherwise
        %%%%%% load answer screen %%%%%%
        cgloadbmp(1,'checker_b.bmp',1280,1024);
        cgloadbmp(2,'memorytest_confidence_received_orig.bmp');
        cgsetsprite(0);
        cgdrawsprite(1,0,0); cgdrawsprite(2,0,0); SOT_cg_c(i) = cgflip(0,0,0);
        response1 = key_SG(1); % their response is whatever key they pressed.
        rt_SG = t_SG(1)/1000 - responset_C; rt_SG = rt_SG*1000; % and their reaction time is the time they pressed the button-the time the stimulus apprered
        cc=cc+1; SOT_r(cc,1) = tmp+rt_SG; SOT_cg_c = tmp+rt_SG; clear tmp
        wait(2000-rt_SG)
    end
    
    eventmarker = eventmarker+1;
    results.presentation{eventmarker,1} = 'response_SG';
    
    % record confidence
    if key_SG == SureKey
        confi = 1;
    elseif key_SG == GuessKey
        confi = 0;
    end
    
    % ----------------------------------------------------------------- %
    
    
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
        
   
    % record
    if n_ON == 0 || n_SG ~= 0
        results.keypress(trl,1:2)   = [NaN key_SG];
        results.oldnew_resp(trl,1)  = NaN;
        results.accuracy(trl,1)     = NaN;
        results.confidence(trl,1)   = confi;
        results.rt(trl,1:2)           = [NaN rt_SG];
    elseif n_ON ~= 0 || n_SG == 0
        results.keypress(trl,1:2)   = [key_ON NaN];
        results.oldnew_resp(trl,1)  = resp;
        results.accuracy(trl,1)     = accu;
        results.confidence(trl,1)   = NaN;
        results.rt(trl,1:2)           = [rt_ON NaN];
    else
        results.keypress(trl,1:2)   = [key_ON key_SG];
        results.oldnew_resp(trl,1)  = resp;
        results.accuracy(trl,1)     = accu;
        results.confidence(trl,1)   = confi;
        results.rt(trl,1:2)           = [rt_ON rt_SG];
    end
    results.trl(trl,:)          = stimvec(:,trl)';
    results.all = horzcat(results.oldnew_resp,results.accuracy,results.confidence,results.rt);
    
    % calculate cumulative SOT
    results.SOT.raw.fix  = SOT_f1;
    results.SOT.raw.stim = SOT_s;
    results.SOT.raw.resp = SOT_r;
    results.SOT.raw.conf = SOT_c;
    results.SOT.cgflip.stim = SOT_cg_s;
    results.SOT.cgflip.resp = SOT_cg_r;
    results.SOT.cgflip.fix  = SOT_cg_f1;
    results.SOT.cgflip.conf = SOT_cg_c;   
    
    % intermediate save
    dat.emomem.test.labels = {'results:' 'response(old-1/new-2)' 'accuracy' 'confidence rating' 'rt'};
    dat.emomem.test.results.SOT    = results.SOT;
    dat.emomem.test.results.resp   = results.oldnew_resp;
    dat.emomem.test.results.accu   = results.accuracy;
    dat.emomem.test.results.confi  = results.confidence;
    dat.emomem.test.results.rt     = results.rt;
    dat.emomem.test.results.all    = results.all;
    dat.emomem.test.results.trl    = results.trl;
    dat.emomem.test.ETinfo         = eyetracker_filename;
    dat.emomem.test.config.keymap  = map;
    
    save([paths.behav num2str(ID) '_EM_backup.mat'],'dat')
    
end

% wrap up

% calculate final cumulative SOT
results.SOT.raw.fix  = SOT_f1;
results.SOT.raw.stim = SOT_s;
results.SOT.raw.resp = SOT_r;
results.SOT.raw.conf = SOT_c;
results.SOT.cgflip.stim = SOT_cg_s;
results.SOT.cgflip.resp = SOT_cg_r;
results.SOT.cgflip.fix  = SOT_cg_f1;
results.SOT.cgflip.conf = SOT_cg_c;

% save data
dat.emomem.test.labels = {'stimlist_all:' 'filenames' 'Indoor(1)/Outdoor(2)' 'Old(1,3,5,7)/New(2,4,6,8)' 'Emo(1,2,3,4)/Neu(5,6,7,8)'; ...
    'results:' 'response(old-1/new-2)' 'accuracy' 'confidence rating(1=sure,0=guess)' 'rt'};
dat.emomem.test.ETinfo         = eyetracker_filename;
dat.emomem.test.config.keymap  = map;
dat.emomem.test.SOT            = results.SOT;
dat.emomem.test.results.resp   = results.oldnew_resp;
dat.emomem.test.results.accu   = results.accuracy;
dat.emomem.test.results.confi  = results.confidence;
dat.emomem.test.results.rt     = results.rt;
dat.emomem.test.results.all    = results.all;
dat.emomem.test.results.trl    = results.trl;

% stats for future convenience
% % accuracies
% dat.emomem.test.results.quickstat.MeanAccuracy = nanmean(results.accuracy);
% dat.emomem.test.results.quickstat.MeanAccuracy_Old = ...
%     nanmean(results.accuracy(cell2mat(stimlist(:,4))==1));
% dat.emomem.test.results.quickstat.MeanAccuracy_New = ...
%     nanmean(results.accuracy(cell2mat(stimlist(:,4))==2));
% dat.emomem.test.results.quickstat.MeanAccuracy_Rew = ...
%     nanmean(results.accuracy(cell2mat(stimlist(:,2))==rewcond));
% dat.emomem.test.results.quickstat.MeanAccuracy_Pun = ...
%     nanmean(results.accuracy(cell2mat(stimlist(:,2))==puncond));
% % confidences
% dat.emomem.test.results.quickstat.MeanConfidence_Old = ...
%     nanmean(results.confidence(cell2mat(stimlist(:,4))==1));
% dat.emomem.test.results.quickstat.MeanConfidence_New = ...
%     nanmean(results.confidence(cell2mat(stimlist(:,4))==2));
% dat.emomem.test.results.quickstat.MeanConfidence_Rew = ...
%     nanmean(results.confidence(cell2mat(stimlist(:,2))==rewcond));
% dat.emomem.test.results.quickstat.MeanConfidence_Pun = ...
%     nanmean(results.confidence(cell2mat(stimlist(:,2))==puncond));
% dat.emomem.test.results.quickstat.MedianConfidence_Old = ...
%     nanmedian(results.confidence(cell2mat(stimlist(:,4))==1));
% dat.emomem.test.results.quickstat.MedianConfidence_New = ...
%     nanmedian(results.confidence(cell2mat(stimlist(:,4))==2));
% dat.emomem.test.results.quickstat.MedianConfidence_Rew = ...
%     nanmedian(results.confidence(cell2mat(stimlist(:,2))==rewcond));
% dat.emomem.test.results.quickstat.MedianConfidence_Pun = ...
%     nanmedian(results.confidence(cell2mat(stimlist(:,2))==puncond));
% % FA rates
% dat.emomem.test.results.quickstat.FA = ...
%     nanmean(double(results.accuracy==0 & cell2mat(stimlist(:,4))==1));
% dat.emomem.test.results.quickstat.FA_Rew = ...
%     nanmean(double(results.accuracy==0 & ...
%     cell2mat(stimlist(:,4))==1 & cell2mat(stimlist(:,2))==rewcond));
% dat.emomem.test.results.quickstat.FA_Pun = ...
%     nanmean(double(results.accuracy==0 & ...
%     cell2mat(stimlist(:,4))==1 & cell2mat(stimlist(:,2))==puncond));

save([paths.behav num2str(ID) '_EM.mat'],'dat')

% terminate protocol
cd(paths.supports)
cgloadbmp(1,'instruction_ende.bmp');
cgsetsprite(0);
cgdrawsprite(1,0,0); 
cgflip(0,0,0); wait(3000);
cgshut;


fprintf('\n%%%%%%%%%%%%%%%%%%%%%%\n COMPLETED  \n%%%%%%%%%%%%%%%%%%%%%%\n')
fprintf('\nMean Accuracy: %3.3f\n', nanmean(results.accuracy))


end
