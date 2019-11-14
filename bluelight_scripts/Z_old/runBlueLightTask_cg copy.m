%% BLUE LIGHT MEMORY TASK : Scene classification task with windows/no windows images
%  In emergency contact Alex:
%              +49 162 713 20 49
%              yeo-jin.yi.15@ucl.ac.uk

%%  work log

%   23_01_2019 created the script

%% START

function [dat] = runBlueLightTask_cg(ID,AgeGroup,EyeFlag,MRIFlag,eyetracker_filename,path_ECHO)

global eyeflag

%% basic setups

% rng(0)

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
paths.stim     = [paths.parent 'stim\stim\'];
paths.pupil    = [paths.parent 'data\pupil\'];
paths.behav    = [paths.parent 'data\behav\'];
% paths.cogent   = 'E:\Dorothea\Cogent2000v1.33';
paths.cogent   = [paths.parent 'Cogent2000v1.33\'];
% paths.cogent = 'C:\Users\ibmi\Downloads\Desktop\pilot_7T\Cogent2000v1.33';
addpath(genpath(paths.cogent))

% create the data structure
cd(paths.behav)
if exist([num2str(ID) '_BL.mat']) == 2
    load([paths.behav num2str(ID) '_BL.mat'])
else
    dat = [];
    dat.ID       = ID;
    dat.AgeGroup = AgeGroup;
    dat.bluelight = [];
    dat.bluelight.task.ETinfo  = fname_eyetracker;
end
load([paths.parent 'stiminfo_BL.mat']);

% scanner preparation
if MRIFlag == 1
    dat.bluelight.task.MRIinfo   = [];
    dat.bluelight.task.MRIinfo.scanPort = 1;
    dat.bluelight.task.MRIinfo.dummy    = 5;        % no. of dummy vols
    dat.bluelight.task.MRIinfo.nslice   = 51;       % no. of slices
    dat.bluelight.task.MRIinfo.TE       = 32;       % time for each slice in msec
    dat.bluelight.task.MRIinfo.TR       = 3.6; % MRIinfo.nslice * MRIinfo.TE/1000; % time for volume in msec
    %     dat.bluelight.task.MRIinfo.nvols    = ceil(MRIinfo.dummy + 5 + ...  % no. of vols per run (5 dummy + 400 + 5 post)
    %         (options.scanblocklength * options.meanTrialLength) / MRIinfo.TR);
    %     dat.bluelight.task.MRIinfo.total_slices = MRIinfo.nvols*MRIinfo.nslice;    % per run!
    %     [MRIinfo] = runScannerSetup_MRIpilot; % redundant
end

%% experimental setups

% condition markers
windows   = 1;  nowindows = 2;
bri_blue  = 11; bri_oran  = 21;
dar_blue  = 12; dar_oran  = 22;

% timings
config.timing.intermission = 5000;
templ = repmat([1:3],1,40);
config.timing.intermission_jitter = templ(randperm(length(templ)))*1000;clear templ
config.timing.cue      = 1500;
config.timing.stim     = 4500;
config.timing.response = 2000;

% config.timing.intermission = 0;
% config.timing.intermission_jitter = zeros(1,40);
% config.timing.cue      = 0;
% config.timing.stim     = 0;
% config.timing.response = 0;

% trials
config.numtrials.total          = 120;
config.numtrials.dark_blues     = 30;
config.numtrials.dark_oranges   = 30;
config.numtrials.bright_blues   = 30;
config.numtrials.bright_oranges = 30;
config.numtrials.windows        = 60;
config.numtrials.no_windows     = 60;

%% roster

% tmp1 = [ ones(1,config.numtrials.bright_blues)*11 ones(1,config.numtrials.dark_blues)*12 ]; tmp1 = tmp1(randperm(length(tmp1))); % blue markers
% tmp1_1 = [ ones(1,config.numtrials.windows*2) ones(1,config.numtirals.no_windows*2)*2 ]; tmp1_1 = tmp1_1(randperm(length(tmp1_1))); % windows marker
% tmp_blue = vertcat(tmp1,tmp1_1);

% tmp2 = [ ones(1,config.numtrials.bright_oranges)*21 ones(1,config.numtrials.dark_oranges)*22 ]; tmp2 = tmp2(randperm(length(tmp2)));
% tmp2_1 = [ ones(1,config.numtrials.windows*2) ones(1,config.numtirals.no_windows*2)*2 ]; tmp2_1 = tmp2_1(randperm(length(tmp2_1)));
% tmp_orange = vertcat(tmp2,tmp2_1);

% randomise cover task contingency
tmp1_1 = [ ones(1,config.numtrials.windows) ones(1,config.numtrials.no_windows)*2 ];
contingencies = tmp1_1(randperm(length(tmp1_1)));

% randomise colour blocks
tmp_base = [11 12 21 22];
tmp1 = tmp_base(randperm(length(tmp_base)));
tmp2 = tmp_base(randperm(length(tmp_base)));

tmp_colourblk = tmp1;
blkcount=1;
while blkcount < 10
    if tmp_colourblk(end) == tmp2(1)
        clear tmp2
        tmp2 = tmp_base(randperm(length(tmp_base)));
        continue
    else
        blkcount = blkcount+1;
        tmp_colourblk = horzcat(tmp_colourblk,tmp2);
    end
end

try % older versions of matlab doesn't have repelem
    colour_blk = repelem(tmp_colourblk,3);
catch
    multipliers = ones(1,3);
    colour_blk = kron(tmp_colourblk,multipliers);
end
stimvec = vertcat(colour_blk,contingencies);

temp_stimlist = [];
wc = 0; nwc = 0;
for sc = 1:length(contingencies)
    if contingencies(sc) == windows
        wc = wc+1;
        temp_stimlist{1,sc} = colour_blk(1,sc);
        temp_stimlist{2,sc} = contingencies(sc);
        temp_stimlist{3,sc} = stiminfo_BL.windows_old{wc};
    elseif contingencies(sc) == nowindows
        nwc = nwc+1;
        temp_stimlist{1,sc} = colour_blk(1,sc);
        temp_stimlist{2,sc} = contingencies(sc);
        temp_stimlist{3,sc} = stiminfo_BL.no_windows_old{nwc};
    end
end

% trial vectors with intermission screens
interm_index  = repmat([1 1 1 0],1,40);
intermvec1 = zeros(1,40);
intermvec2 = repmat([1 -1],1,20);
intermissions = vertcat(intermvec2,intermvec1);
tmptrlvec = zeros(2,160);
tmptrlvec(:,interm_index==1) = stimvec;
tmptrlvec(:,interm_index==0) = intermissions;
trlvec = tmptrlvec; clear tmptrlvec

config.stim.trialvec = trlvec';
config.stim.stimlist = temp_stimlist';
% config.stim.intermissions = (-1).^(0:140);


% backgrounds
grey_bright = 'grey_bright.bmp';
grey_dark   = 'grey_dark.bmp';
blue_bright = 'blue_bright.bmp';
blue_dark   = 'blue_dark.bmp';
oran_bright = 'orange_bright.bmp';
oran_dark   = 'orange_dark.bmp';


fprintf('\n\nTASK SETUP DONE')

%% run

% configs
config_display(0, 5, [0 0 0], [1 1 1], 'Arial', 20, 10, 0); % 1024x768
config_sound;
config_keyboard;

% scanner input
if MRIFlag == 1
    config_serial(dat.bluelight.task.MRIinfo.scanPort); % links to scanner for waitslice
end

monitor = get(0,'screensize');
bg_w = monitor(3); bg_h = monitor(4);

%%%%%%%%%%%%
start_cogent
%%%%%%%%%%%%

% cgopen(5,0,0,2) % at the MRI, choose the secondary screen
cgopen(5,0,0,1)

map         = getkeymap; % define keyboard IDs
% WinKey      = map.A;
% noWinKey    = map.L;
WinKey      = 28; %%%%% check again
noWinKey    = 29;

% print instructions
cd(paths.supports)
cgloadbmp(1,'emomem_inst_task_orig.bmp');
cgsetsprite(0);
cgdrawsprite(1,0,0);
t0_inst = cgflip(0,0,0);
% waitkeydown(inf)
waitkeydown(inf,[WinKey,noWinKey]); clearkeys

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
    
    %%%%%%%%%%%%%%
    tic % ding ding ding
    %%%%%%%%%%%%%%
    
    scannerinput = map.T;
    [k, dat.bluelight.task.MRIinfo.dummyscanfirst] = waitkeydown(inf,scannerinput); % the first dummy scan
    for dum = 1:(dat.bluelight.task.MRIinfo.dummy - 2)
        waitkeydown(inf,scannerinput);
    end
    [k, dat.bluelight.task.MRIinfo.dummyscanlast] = waitkeydown(inf,scannerinput); % the last dummy scan
    
%     %%%%%%%%%%%%%%
%     tic % ding ding ding
%     %%%%%%%%%%%%%%
    
    cd(paths.supports); cgloadbmp(1,grey_bright,1280,1024);
    cgsetsprite(0);
    cgdrawsprite(1,0,0);
    cgfont('Arial',300);
    cgpencol(.2,.2,.2)
    cgtext('+',0,0);  % Prepare buffer 1 (fixation cross)
    t0_fix0 = cgflip(0,0,0);  t0_fix0_raw = toc; % load fixation cross
    %     if eyeflag==1
    %         Eyelink('Message',[num2str(0) 'f0'])  % send trialnumber and stimulustype to eyetracker
    %     end
    wait(dat.bluelight.task.MRIinfo.TR*1000)
    clearpict(1);
    
else
    
    %%%%%%%%%%%%%%
    tic % ding ding ding
    %%%%%%%%%%%%%%
    
    cd(paths.supports); cgloadbmp(1,grey_bright,1280,1024);
    cgsetsprite(0);
    cgdrawsprite(1,0,0);
    cgfont('Arial',300);
    cgpencol(.2,.2,.2)
    cgtext('+',0,0);  % Prepare buffer 1 (fixation cross)
    cgflip(0,0,0);             % load fixation cross
    %     if eyeflag==1
    %         Eyelink('Message',[num2str(0) 'f0'])  % send trialnumber and stimulustype to eyetracker
    %     end
    wait(1000)
    
end

results = []; results.SOT = [];
% results.SOT.cgflip.t0_fix0 = t0_fix0;
% results.SOT.cgflip.t0_standby = t0_standby;
% results.SOT.cgflip.t0_inst = t0_inst;
% results.SOT.raw.t0_fix0 = t0_fix0_raw;
eventmarker = 0; interm_cnt = 0;
i = 0; f1c = 0; sc = 0; rc1 = 0; rc2 = 0; vc = 0; cc = 0; ic = 0;
SOT_f1=[]; SOT_s=[]; SOT_r=[]; SOT_v = []; SOT_c = []; SOT_i = [];
SOT_cg_s=[];SOT_cg_f1=[];SOT_cg_r=[]; SOT_cg_v=[]; SOT_cg_c = []; SOT_cg_i = [];
for trl = 1:size(trlvec,2)
    
    if trlvec(2,trl) == 0 % if it is an intermission
        
        % -------------------- present intermission ---------------------- %
        
        interm_cnt = interm_cnt+1;
        
        if trlvec(1,trl) > 0 % 1 is dark grey
            cuename = grey_dark;
        elseif trlvec(1,trl) < 0 % -1 is bright grey
            cuename = grey_bright;
        end
        
        cd(paths.supports); cgloadbmp(1,cuename,1280,1024);
        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgfont('Arial',200)
        cgpencol(.2,.2,.2)
        cgtext('+',0,0);
        ic=ic+1; SOT_i(ic,1)=toc;
        SOT_cg_i(f1c) = cgflip(0,0,0);             % load fixation cross
        if eyeflag==1
            Eyelink('Message',[num2str(trl) 'i1'])  % send trialnumber and stimulustype to eyetracker
        end
        if eyeflag==1
            Eyelink('Message',[num2str(trl) 'i2'])  % send trialnumber and stimulustype to eyetracker
        end
        if eyeflag==1
            Eyelink('Message',[num2str(trl) 'i3'])  % send trialnumber and stimulustype to eyetracker
        end
        if eyeflag==1
            Eyelink('Message',[num2str(trl) 'i4'])  % send trialnumber and stimulustype to eyetracker
        end
        wait(config.timing.intermission + config.timing.intermission_jitter(interm_cnt))
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = 'intermission';
        
        % record
        results.SOT.raw.intermission = SOT_i;
        results.SOT.cgflip.intermission = SOT_cg_i;
        
        clear cuename
        
        % ----------------------------------------------------------------- %
        
    else
        
        % ------------------------ present image -------------------------- %
        i = i+1;
        fprintf('\nTrial %d\n', i)
        
        clearkeys;
        
        %%%%%% load cue %%%%%%
        if trlvec(1,trl) == 11     % bright blue
            cuename = blue_bright;
        elseif trlvec(1,trl) == 12 % dark blue
            cuename = blue_dark;
        elseif trlvec(1,trl) == 21 % bright orange
            cuename = oran_bright;
        elseif trlvec(1,trl) == 22 % dark orange
            cuename = oran_dark;
        end
        cd(paths.supports); cgloadbmp(1,cuename,1280,1024);
        
        %%%%%% draw cue screen %%%%%%
        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cc = cc+1; SOT_c(cc,1)=toc;
        SOT_cg_c(cc) = cgflip(0,0,0);
        if eyeflag==1
            Eyelink('Message',[num2str(trl) 'c'])  % send trialnumber and stimulustype to eyetracker
        end
        wait(config.timing.cue)
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = ['cue_' cuename ];
        
        %%%%%% load background %%%%%%
        cd(paths.supports);
        if trlvec(1,trl) == 11     % bright blue
            cgloadbmp(1,blue_bright,1280,1024);
        elseif trlvec(1,trl) == 12 % dark blue
            cgloadbmp(1,blue_dark,1280,1024);
        elseif trlvec(1,trl) == 21 % bright orange
            cgloadbmp(1,oran_bright,1280,1024);
        elseif trlvec(1,trl) == 22 % dark orange
            cgloadbmp(1,oran_dark,1280,1024);
        end
        
        %%%%%% load stimulus %%%%%%
        cd(paths.stim)
        stimname = config.stim.stimlist{i,3};
        cgloadbmp(2,stimname,400,300);
        categories(1:2,i) = trlvec(1:2,trl);
        
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
        
        
        % --------------------- YES/NO response screen -------------------- %
        
        %%%%% load IN/OUT Question screen %%%%%%
        cd(paths.supports); cgloadbmp(1,cuename,1280,1024);
%         cgloadbmp(2,'response_orig.bmp');
        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgfont('Arial',30)
        cgpencol(0,0,0)
        cgtext('Waren Fenster im Bild zu sehen?',0,50);
        cgtext('  JA           NEIN',0,10)
%         cgdrawsprite(2,0,0);
        tmp = toc;
        responset1 = cgflip(0,0,0);
        if eyeflag==1
            Eyelink('Message',[num2str(trl) 'r'])  % send trialnumber and stimulustype to eyetracker
        end
        waitkeydown(config.timing.response,[WinKey,noWinKey]);
        
        % RT
        [keypress1, t1, n1] = getkeydown;
        if n1 == 0; %if they didn't press anything
            response1 = NaN; % mark their response as nan
            keypress1 = NaN;
            rt1 = NaN; % mark their reaction time as nan
            wait(0);
        else % otherwise
            %%%%%% load answer screen %%%%%%
            cd(paths.supports); cgloadbmp(1,cuename,1280,1024);
%             cgloadbmp(2,'response_received_orig.bmp');
            cgsetsprite(0);
            cgdrawsprite(1,0,0); 
%             cgdrawsprite(2,0,0); 
            cgfont('Arial',30)
            cgpencol(0,0,0)
            cgtext('Waren Fenster im Bild zu sehen?',0,50);
            cgtext('  JA           NEIN',0,10)
            cgtext('V',0,-30)
            cgflip(0,0,0);
            response1 = keypress1(1); % their response is whatever key they pressed.
            rt1 = t1(1)/1000 - responset1; rt1 = rt1*1000; % and their reaction time is the time they pressed the button-the time the stimulus apprered
            rc1=rc1+1; SOT_r(rc1,1) = tmp+rt1; SOT_cg_r(rc1,1) = tmp+rt1; clear tmp
            wait(2000-rt1)
        end
        
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = 'response';
        
        % accuracy
        if (response1 == WinKey && trlvec(2,trl) == 1)
            accuracy = 1;
            fprintf('\nCorrect \n')
        elseif (response1 == noWinKey && trlvec(2,trl) == 2)
            accuracy = 1;
            fprintf('\nCorrect \n')
        elseif (response1 == WinKey  && trlvec(2,trl) == 2 && ~isnan(response1)) ...
                || (response1 == noWinKey && trlvec(2,trl) == 1 && ~isnan(response1))
            accuracy = 0;
            fprintf('\nIncorrect \n')
        elseif isnan(response1)
            accuracy = 0;
            fprintf('\nMissed \n')
        else
            accuracy = NaN;
        end
        
        % ----------------------------------------------------------------- %
        
        % ------------------------ present fix x -------------------------- %
        
        %     if config.stim.intermissions(trl) > 0 % 1 is dark grey
        %         cd(paths.supports); cgloadbmp(1,grey_dark,1280,1024);
        %     elseif config.stim.intermissions(trl) < 0 % -1 is bright grey
        %         cd(paths.supports); cgloadbmp(1,grey_bright,1280,1024);
        %     end
        cd(paths.supports); cgloadbmp(1,cuename,1280,1024);
        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgfont('Arial',200)
        cgpencol(.5,.5,.5)
        cgtext('+',0,0);
        f1c=f1c+1; SOT_f1(f1c,1)=toc;
        SOT_cg_f1(f1c) = cgflip(0,0,0);             % load fixation cross
        if eyeflag==1
            Eyelink('Message',[num2str(trl) 'f'])  % send trialnumber and stimulustype to eyetracker
        end
        wait(500)
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = '+';
        
        % ----------------------------------------------------------------- %
        
        % record
        try
            results.keypress(i,1) = keypress1; % response
        catch
            try
                results.keypress(i,1) = keypress1(1);
            catch
                results.keypress(i,1) = NaN;
            end
        end
        results.rt_resp(i,1)  = rt1;
        results.accuracy(i,1) = accuracy;
        results.trl{i,1}      = stimname;
        results.trl{i,2}      = categories(1,i);
        results.trl{i,3}      = categories(2,i);
        
        
        % calculate cumulative SOT
        results.SOT.raw.fix  = SOT_f1;
        results.SOT.raw.stim = SOT_s;
        results.SOT.raw.resp = SOT_r;
        results.SOT.raw.cue  = SOT_c;
        results.SOT.cgflip.stim = SOT_cg_s;
        results.SOT.cgflip.resp = SOT_cg_r;
        results.SOT.cgflip.fix  = SOT_cg_f1;
        results.SOT.cgflip.cue  = SOT_cg_c;
        
        % intermediate save
        dat.bluelight.task.results = results;
        save([paths.behav num2str(ID) '_BL_backup.mat'],'dat')
        
        clear stimname cuename
        
    end % close trial-intermission conditional
    
    
end % close trial loop

%% wrap up

% calculate cumulative SOT
results.SOT.raw.fix  = SOT_f1;
results.SOT.raw.stim = SOT_s;
results.SOT.raw.resp = SOT_r;
results.SOT.raw.cue  = SOT_c;
results.SOT.raw.intermission = SOT_i;

results.SOT.cgflip.stim = SOT_cg_s;
results.SOT.cgflip.resp = SOT_cg_r;
results.SOT.cgflip.fix  = SOT_cg_f1;
results.SOT.cgflip.cue  = SOT_cg_c;
results.SOT.cgflip.intermission = SOT_cg_i;

% save results
dat.bluelight.task.results = results;
dat.bluelight.task.labels  = {'contingencies-> 11=brightBlue, 12=darkBlue, 21=brightOrange, 22=darkOrange / -1=brightGrey, 1=darkGrey / 1=windows, 2=no windows'; ...
    'results.trl: filenames, contingencies'};
dat.bluelight.task.config.keymap = map;
dat.bluelight.task.config = config;
save([paths.behav num2str(ID) '_BL.mat'],'dat')

% terminate protocol
cd(paths.supports)
cgloadbmp(1,'instruction_ende.bmp');
cgsetsprite(0);
cgdrawsprite(1,0,0);
cgflip(0,0,0); wait(3000);
cgshut;


end
