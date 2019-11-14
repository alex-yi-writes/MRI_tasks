%% Face oddball (M/F)
%  In emergency contact Alex:
%              +49 162 713 20 49
%              yeo-jin.yi.15@ucl.ac.uk

%%  work log

%   30_12_2018      created the script

%% START

function [dat] = runOddball_cg(ID,AgeGroup,EyeFlag,MRIFlag,eyetracker_filename,Schedules)

global eyeflag fname_eyetracker

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
        %
        %         % STEP 1
        %         % Open a graphics window on the main screen
        %         % using the PsychToolbox's Screen function.
        %         Screen('Preference', 'SkipSyncTests', 1)
        %         screenNumber=max(Screen('Screens'));
        %         window=Screen('OpenWindow', screenNumber);
        %
        %         % STEP 2
        %         % Provide Eyelink with details about the graphics environment
        %         % and perform some initializations. The information is returned
        %         % in a structure that also contains useful defaults
        %         % and control codes (e.g. tracker state bit and Eyelink key values).
        %         el=EyelinkInitDefaults(window);
        %
        %         % Disable key output to Matlab window:
        %         ListenChar(2);
        %
        %         % STEP 3
        %         % Initialization of the connection with the Eyelink Gazetracker.
        %         % exit program if this fails.
        %         if ~EyelinkInit(dummymode, 1)
        %             fprintf('Eyelink Init aborted.\n');
        %             %         cleanup;  % cleanup function
        %             return;
        %         end
        %
        %         [v vs]=Eyelink('GetTrackerVersion');
        %         fprintf('Running experiment on a ''%s'' tracker.\n', vs );
        %
        %         % make sure that we get gaze data from the Eyelink
        %         Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');
        
        % open file to record data to
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
paths.parent   = '\\fs-md\users\yiy\Documents\DOPE\Oddball\';
paths.stim     = [paths.parent 'stim\'];
paths.pupil    = [paths.parent 'data\pupil\'];
paths.behav    = [paths.parent 'data\behav\'];
paths.cogent   = 'C:\Studies\Grid_experiment2\Cogent2000v1.32';
addpath(genpath(paths.cogent))

% create the data structure
cd(paths.behav)
if exist([num2str(ID) '_OB.mat']) == 2
else
    dat = [];
    dat.ID       = ID;
    dat.AgeGroup = AgeGroup;
    dat.oddball = [];
    dat.oddball.ETinfo  = fname_eyetracker;
    dat.oddball.StimSchedule = Schedules;
end
if Schedules == 1
    load([paths.parent 'scripts\OB_GA_designstruct_ver1.mat'])
elseif Schedules == 2
    load([paths.parent 'scripts\OB_GA_designstruct_ver2.mat'])
end

% scanner preparation
if MRIFlag == 1
    dat.oddball.MRIinfo   = [];
    dat.oddball.MRIinfo.scanPort = 1;
    dat.oddball.MRIinfo.dummy    = 5;        % no. of dummy vols
    dat.oddball.MRIinfo.nslice   = 51;       % no. of slices
    dat.oddball.MRIinfo.TE       = 32;       % time for each slice in msec
    dat.oddball.MRIinfo.TR       = 3.6; % MRIinfo.nslice * MRIinfo.TE/1000; % time for volume in msec
    %     dat.oddball.MRIinfo.nvols    = ceil(MRIinfo.dummy + 5 + ...  % no. of vols per run (5 dummy + 400 + 5 post)
    %         (options.scanblocklength * options.meanTrialLength) / MRIinfo.TR);
    %     dat.oddball.MRIinfo.total_slices = MRIinfo.nvols*MRIinfo.nslice;    % per run!
    %     [MRIinfo] = runScannerSetup_MRIpilot; % redundant
end

%% experimental setups

% timings
fixjit = repmat(500:50:700,1,24);
config.timing.fixation = [(design_struct.eventlist(:,4)+design_struct.eventlist(:,5)-2)']*1000;
config.timing.stim     = 2500;
config.timing.response = 2000; % maximum waiting time for the button press

% --- for breezeing through --- %
% config.timing.fixation = zeros(1,240);
% config.timing.stim     = 200;
% config.timing.response = 200;
% ----------------------------- %

% trials
config.numtrials.total   = 260;
config.numtrials.oddball = 80;
config.numtrials.standard= 160;
config.numtrials.null    = 20;

% stimulus vector
% trlvec   = randperm(config.numtrials.total);
std      = ones(1,config.numtrials.standard);
ob       = ones(1,config.numtrials.oddball).*2;
nulls    = zeros(1,config.numtrials.null);
sex_std  = horzcat(ones(1,config.numtrials.standard/2),...
    ones(1,config.numtrials.standard/2).*2); sex_std = sex_std(randperm(config.numtrials.standard));
sex_ob   = 1:config.numtrials.oddball;
sex_null = zeros(1,config.numtrials.null);
% stimtype = horzcat(std,ob);
% filetype = horzcat(f_std,f_ob);
% stimvec  = vertcat(trlvec,stimtype,filetype);
trlvec   = 1:config.numtrials.total;
stimtype = design_struct.eventlist(:,3)'; % pre-generated presentation schedule (from Tor Wager's script)
filetype = zeros(1,config.numtrials.total);
filetype(find(stimtype==2)) = sex_std;
filetype(find(stimtype==1)) = sex_ob;
filetype(find(stimtype==0)) = sex_null;
stimvec = vertcat(trlvec,stimtype,filetype);

% roster
cd(paths.stim)
tmp = dir('o*.bmp'); fname = [];
[fname.ob{1:config.numtrials.oddball}] = deal(tmp(1:config.numtrials.oddball).name); clear tmp
tmp = dir('s*.bmp');
[fname.std{1:2}] = deal(tmp.name); clear tmp

% insert break
tmp = horzcat(stimvec(:,1:130),[777;NaN;NaN],stimvec(:,131:end));
stimvec = tmp;

config.stim.stimvec  = stimvec;
config.stim.filelist = fname;
stimvec_loop = stimvec(1,:);
stimtype_loop = stimvec(2,:);
filetype_loop = stimvec(3,:);

dat.oddball.config = config;

%% run

% configs
config_display(1, 4, [0 0 0], [1 1 1], 'Arial', 20, 10, 0); % 1024x768
% config_display(screenMode, screenRes, bg ,fontcol, fontName, fontSize, number_of_buffers,0);   % open graphics window
config_sound;
config_keyboard;

monitor = get(0,'screensize');
bg_w = monitor(3); bg_h = monitor(4);

% scanner input
if MRIFlag == 1
    config_serial(dat.oddball.MRIinfo.scanPort); % links to scanner for waitslice
end

%%%%%%%%%%%%
start_cogent
%%%%%%%%%%%%

cgopen(4,0,0,2)
% cgopen(4,0,0,1)

map     = getkeymap; % define keyboard IDs
MalKey  = 28;
FemKey  = 29;
buff    = 2;

Picture = 2;% name buffers
Background = 3;
Questions = 4;
Answers = 5;
Instruction = 6;
Breaks = 7;
Perfgr = 8;

cd(paths.parent);
cgloadbmp(1,'oddball_inst.bmp');
cgsetsprite(0);
cgdrawsprite(1,0,0);
cgflip(0,0,0); waitkeydown(inf,[MalKey,FemKey]);

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
    cd(paths.parent)
    cgloadbmp(1,'operator_hail.bmp');
    cgsetsprite(0);
    cgdrawsprite(1,0,0);
    cgflip(0,0,0);
    waitkeydown(inf,map.Space)
    
    % start the dummy scan
    cd(paths.parent)
    cgloadbmp(1,'ready.bmp');
    cgsetsprite(0);
    cgdrawsprite(1,0,0);
    t0_standby=cgflip(0,0,0); 
    scannerinput = 33;
    
    %%%%%%%%%%%%%%
    tic % ding ding ding
    %%%%%%%%%%%%%%

    [k, dat.oddball.MRIinfo.dummyscanfirst] = waitkeydown(inf,scannerinput); % the first dummy scan
    for dum = 1:(dat.oddball.MRIinfo.dummy - 2)
        waitkeydown(inf,scannerinput);
    end
    [k, dat.oddball.MRIinfo.dummyscanlast] = waitkeydown(inf,scannerinput); % the last dummy scan
    
    cd(paths.parent); cgloadbmp(1,'checker_b.bmp',1152,864);
    cgsetsprite(0);
    cgdrawsprite(1,0,0);
    cgfont('Arial',300);
    cgpencol(.5,.5,.5)
    cgtext('+',0,0);  % Prepare buffer 1 (fixation cross)
    t0_fix0 = cgflip(0,0,0); t0_fix0_raw = toc;  % load fixation cross
    wait(dat.oddball.MRIinfo.TR*1000)
    
else
    
    %%%%%%%%%%%%%%
    tic % ding ding ding
    %%%%%%%%%%%%%%
    
    cd(paths.parent); cgloadbmp(1,'checker_b.bmp',1152,864);
    cgsetsprite(0);
    cgdrawsprite(1,0,0);
    cgfont('Arial',300)
    cgpencol(.5,.5,.5)
    cgtext('+',0,0);
    cgflip(0,0,0);
    wait(1000)
end

results = []; results.SOT = []; 
results.SOT.cgflip.t0_fix0 = t0_fix0;
results.SOT.cgflip.t0_standby = t0_standby;
results.SOT.raw.t0_fix0 = t0_fix0_raw;
eventmarker = 0;
i = 0; f1c = 0; sc = 0; rc = 0; nullc = 0;
SOT_f1=[]; SOT_s=[]; SOT_r=[]; SOT_null = [];
SOT_cg_f1=[]; SOT_cg_s=[]; SOT_cg_r=[]; SOT_cg_n=[];
stimcount = 0;
for trl = 1:length(stimvec_loop)
    
    if (trl <= 130) || (trl >= 132)
        stimcount = stimcount+1;        
                
        % ------------------------ present image -------------------------- %
        i = i+1;
        fprintf('\nTrial %d\n', stimcount)
        
        %%%%%% load background %%%%%%
        cd(paths.parent); 
        cgloadbmp(1,'checker_b.bmp',1152,864);

        %%%%%% load stimulus %%%%%%
        cd(paths.stim)
        if stimtype_loop(trl) == 2 || stimtype_loop(trl) == 1
            
            if stimtype_loop(trl) == 2 % if the stim type is standard
                stimname = char(fname.std(1,filetype_loop(trl)));
            elseif stimtype_loop(trl) == 1 % if the stim type is oddball
                stimname = char(fname.ob(1,filetype_loop(trl)));
            end
            stimname
            cgloadbmp(2,stimname,390,530);
            categories(1,i) = stimtype_loop(trl);
            
            %%%%%% draw stimulus screen %%%%%%
            sc=sc+1; SOT_s(sc,1)=toc;
            cgsetsprite(0);
            cgdrawsprite(1,0,0);
            cgdrawsprite(2,0,0);
            SOT_cg_s(sc) = cgflip(0,0,0);
            if eyeflag==1
                Eyelink('Message',[num2str(stimcount) 's'])  % send trialnumber and stimulustype to eyetracker
            end
            wait(config.timing.stim);
            eventmarker = eventmarker+1;
            results.presentation{eventmarker,1} = stimname;
            
            %%%%%% load Question screen %%%%%%
            tmp = toc; % clearkeys; readkeys;
            cd(paths.parent); 
            cgloadbmp(1,'checker_b.bmp',1152,864);
            cgloadbmp(2,'male_female.bmp');
            cgsetsprite(0);
            cgdrawsprite(1,0,0);
            cgdrawsprite(2,0,0); tmp = toc;
            responset = cgflip(0,0,0);
%             if eyeflag==1
%                 Eyelink('Message',[num2str(stimcount) 'r'])  % send trialnumber and stimulustype to eyetracker
%             end
            waitkeydown(config.timing.response,[MalKey,FemKey]);
            
            % RT
            [keypress, t, n] = getkeydown;
            if n == 0; %if they didn't press anything
                response = NaN; % mark their response as nan
                keypress = NaN;
                rt = NaN; % mark their reaction time as nan
                wait(0);
            else % otherwise
                %%%%%% load answer screen %%%%%%
                rc=rc+1;
                cgloadbmp(1,'checker_b.bmp',1152,864);
                cgloadbmp(2,'male_female_received.bmp');
                cgsetsprite(0);
                cgdrawsprite(1,0,0); cgdrawsprite(2,0,0); SOT_cg_r(rc) = cgflip(0,0,0);
                response = keypress(1); % their response is whatever key they pressed.
                rt = t(1)/1000 - responset; rt = rt*1000; % and their reaction time is the time they pressed the button-the time the stimulus apprered
                SOT_r(rc,1) = tmp+rt; clear tmp
                wait(2000-rt)
            end
            
            eventmarker = eventmarker+1;
            results.presentation{eventmarker,1} = 'response';           
            
            % accuracy
            if (response == MalKey && strcmp(stimname(6),'m'))
                accuracy = 1;
                fprintf('\n***** Correct *****\n')
            elseif (response == FemKey && strcmp(stimname(6),'f'))
                accuracy = 1;
                fprintf('\n***** Correct *****\n')
            elseif (response == FemKey  && strcmp(stimname(6),'m') && ~isnan(response)) ...
                    || (response == MalKey && strcmp(stimname(6),'f') && ~isnan(response))
                accuracy = 0;
                fprintf('\n***** Incorrect ***** \n')
            elseif isnan(response)
                accuracy = 0;
                fprintf('\n***** Missed ***** \n')
            end
            
            
        elseif stimtype_loop(trl) == 0 % if the stim type is null
            stimname = 'null';
            cd(paths.parent); cgloadbmp(1,'checker_b.bmp',1152,864);
            cgsetsprite(0);
            cgdrawsprite(1,0,0);
            cgfont('Arial',300)
            cgpencol(.5,.5,.5)
            cgtext('+',0,0);
            nullc=nullc+1; SOT_null(nullc,1)=toc;
            SOT_cg_n(nullc) = cgflip(0,0,0);             % load fixation cross
            if eyeflag==1
                Eyelink('Message',[num2str(stimcount) 'n'])  % send trialnumber and stimulustype to eyetracker
            end
            wait(2500)
            
            response = NaN; % mark their response as nan
            keypress = NaN;
            rt = NaN; % mark their reaction time as nan
            accuracy = NaN;
        end
        
        % ----------------------------------------------------------------- %
        
        
        
        % ------------------------ present fix x -------------------------- %
        
%         cd(paths.parent); loadpict('checker.bmp',1,0,0,bg_w,bg_h);
%         setforecolour(.5,.5,.5);
%         settextstyle('Arial', 300);
%         preparestring('+', 4);  % Prepare buffer 1 (fixation cross)
        
        cd(paths.parent); cgloadbmp(1,'checker_b.bmp',1152,864);
        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgfont('Arial',300)
        cgpencol(.5,.5,.5)
        cgtext('+',0,0);
        f1c=f1c+1; SOT_f1(f1c,1)=toc;
        SOT_cg_f1(f1c) = cgflip(0,0,0);             % load fixation cross
        if eyeflag==1
            Eyelink('Message',[num2str(stimcount) 'f1'])  % send trialnumber and stimulustype to eyetracker
        end
        wait(config.timing.fixation(stimcount))
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = '+';
        
        clearpict(1);
        
        % ----------------------------------------------------------------- %
        
        
        % record
        try
            results.keypress(i,1) = keypress;
        catch
            results.keypress(i,1) = 999;
        end
        results.rt(i,1)       = rt;
        results.accuracy(i,1) = accuracy;
        results.trl{i,1}      = stimname;
        results.trl{i,2}      = stimtype_loop(trl);
        results.trl{i,3}      = filetype_loop(trl);
        
        % calculate cumulative SOT
%         results.SOT.cumulative.fix  = cumsum(SOT_f1);
%         results.SOT.cumulative.stim = cumsum(SOT_s);
%         results.SOT.cumulative.resp = cumsum(SOT_r);
        results.SOT.raw.fix  = SOT_f1;
        results.SOT.raw.stim = SOT_s;
        results.SOT.raw.resp = SOT_r;
        results.SOT.raw.null = SOT_null;
        results.SOT.cgflip.stim = SOT_cg_s;
        results.SOT.cgflip.resp = SOT_cg_r;
        results.SOT.cgflip.fix  = SOT_cg_f1;
        results.SOT.cgflip.null = SOT_cg_n;
        
        % intermediate save
        dat.oddball.results = results;
        save([paths.behav num2str(ID) '_OB_backup.mat'],'dat')
        
        clear stimname
        
    elseif trl == 131
        
        cd(paths.parent)
        cgloadbmp(1, 'instruction_break.bmp');       
        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgflip(0,0,0);
        waitkeydown(inf,[MalKey,FemKey]);
        
    end
end

%% wrap up
% calculate final cumulative SOT
% results.SOT.cumulative.fix  = cumsum(SOT_f1);
% results.SOT.cumulative.stim = cumsum(SOT_s);
% results.SOT.cumulative.resp = cumsum(SOT_r);
results.SOT.raw.fix  = SOT_f1;
results.SOT.raw.stim = SOT_s;
results.SOT.raw.resp = SOT_r;
results.SOT.raw.nulL = SOT_null;
results.SOT.cgflip.stim = SOT_cg_s;
results.SOT.cgflip.resp = SOT_cg_r;
results.SOT.cgflip.fix  = SOT_cg_f1;
results.SOT.cgflip.null = SOT_cg_n;

% save results
dat.oddball.results = results;
dat.oddball.labels  = {'stimtype->standard(2),oddball(1) / filetype->standard(female_1,male_2),oddball(no indication)'; ...
    'results.trl: filenames, stimtype, filetype'};
dat.oddball.config.keymap = map;
save([paths.behav num2str(ID) '_OB.mat'],'dat')

% terminate protocol
% clearpict(1)
% settextstyle('Arial', 35);
% preparestring('Das Ende der Aufgabe',1,0,0)
% drawpict(1); wait(3000);
% stop_cogent;

% terminate protocol
cd(paths.parent)
cgloadbmp(1,'instruction_ende.bmp');
cgsetsprite(0);
cgdrawsprite(1,0,0); waitkeydown(inf);
cgflip(0,0,0); wait(3000);
cgshut;



end



