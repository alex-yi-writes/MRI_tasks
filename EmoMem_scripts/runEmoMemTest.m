%% EMOTIONAL MEMORY TEST : Scene classification task with emotional/neutral images
%  In emergency contact Alex:
%              +49 162 713 20 49
%              yeo-jin.yi.15@ucl.ac.uk

%%  work log

%   23_01_2019 created the script

%% START
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
function [dat] = runEmoMemTest(ID,AgeGroup,EyeFlag,MRIFlag,eyetracker_filename)

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
paths.stim_old = [paths.parent 'stim\old\'];
paths.stim_new = [paths.parent 'stim\new\'];
paths.pupil    = [paths.parent 'data\pupil\'];
paths.behav    = [paths.parent 'data\behav\'];
paths.cogent   = '\\fs-md\users\yiy\Documents\MATLAB\Cogent2000v1.33';
addpath(genpath(paths.cogent))

% create the data structure
cd(paths.behav)
if exist([num2str(ID) '_EM.mat']) == 2
    load([paths.behav num2str(ID) '_EM.mat'])
    dat.emomem.test = [];
else
    error('Previous task results do not exist!!!')
end


% scanner preparation
if MRIFlag == 1
    dat.emomem.test.MRIinfo   = [];
    dat.emomem.test.MRIinfo.scanPort = 1;
    dat.emomem.test.MRIinfo.dummy    = 5;        % no. of dummy vols
    dat.emomem.test.MRIinfo.nslice   = 51;       % no. of slices
    dat.emomem.test.MRIinfo.TE       = 32;       % time for each slice in msec
    dat.emomem.test.MRIinfo.TR       = 3.6; % MRIinfo.nslice * MRIinfo.TE/1000; % time for volume in msec
    dat.emomem.test.MRIinfo.nvols    = ceil(MRIinfo.dummy + 5 + ...  % no. of vols per run (5 dummy + 400 + 5 post)
        (options.scanblocklength * options.meanTrialLength) / MRIinfo.TR);
    dat.emomem.test.MRIinfo.total_slices = MRIinfo.nvols*MRIinfo.nslice;    % per run!
    %     [MRIinfo] = runScannerSetup_MRIpilot; % redundant
end

%% experimental setups

config = [];

% condition markers
in_emo = 11; in_neu = 12;
out_emo = 21; out_neu = 22;

% timing
config.timing.outcome  = 2500;
fixjit1 = repmat(1000:250:2000,1,48); fixjit2 = repmat(1000:250:2000,1,48);
config.timing.fixation1= fixjit1(randperm(length(fixjit1)));
config.timing.fixation2= fixjit2(randperm(length(fixjit2)));
config.timing.cue      = 2500;
config.timing.resp     = 2000;

% --- for breezeing through --- %
% config.timing.fixation1= zeros(1,240);
% config.timing.fixation2= zeros(1,240);
% config.timing.cue      = 1;
% config.timing.outcome  = 1;
% ----------------------------- %

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


% recall phase stimuli list
config.stim.new = vertcat( horzcat(fname.new_in_emo,fname.new_in_neu,fname.new_out_emo,fname.new_out_neu),...
    num2cell( ones(1,config.numtrials.new_total)*2 ),...
    num2cell( horzcat( ones(1,config.numtrials.new_in_emo)*in_emo,ones(1,config.numtrials.new_in_neu)*in_neu,...
    ones(1,config.numtrials.new_out_emo)*out_emo,ones(1,config.numtrials.new_out_neu)*out_neu ) ));
config.stim.old = vertcat( horzcat(fname.old_in_emo,fname.old_in_neu,fname.old_out_emo,fname.old_out_neu),...
    num2cell( ones(1,config.numtrials.old_total) ),...
    num2cell( horzcat( ones(1,config.numtrials.old_in_emo)*in_emo,ones(1,config.numtrials.old_in_neu)*in_neu,...
    ones(1,config.numtrials.old_out_emo)*out_emo,ones(1,config.numtrials.old_out_neu)*out_neu ) ));
config.stim.all = horzcat(config.stim.new,config.stim.old);


% scramble
stimvec = config.stim.all;
stimvec = stimvec(:,randperm(size(stimvec,2)));

config.stim.stimvec = stimvec;
dat.emomem.test.config = config;


%% run

% configs
config_display(1, 6, [0 0 0], [1 1 1], 'Arial', 35, 4); % 1280*1024
config_sound;
config_keyboard;

monitor = get(0,'screensize');
bg_w = monitor(3); bg_h = monitor(4);

%%%%%%%%%%%%
start_cogent
%%%%%%%%%%%%

map     = getkeymap; % define keyboard IDs
OldKey  = map.B;
NewKey  = map.N;
Surekey = map.B;
Guesskey= map.N;
buff    = 4;

% print instructions
cd(paths.supports)
loadpict('memorytest_instruction1.png',1);
drawpict(1); waitkeydown(inf);

if eyeflag==1
    Eyelinkerror = Eyelink('StartRecording'); % start recording to the file
    if Eyelinkerror ~= 0
        return
        error ('Eyetracker failed')
    end
    Eyelink('Message',num2str([0 100]))
end

clear trl c1 c2
c1 = 0; c2 = 0; results = [];
i = 0; f1c = 0; sc = 0; rc1 = 0; rc2 = 0; vc = 0;
resp = []; accu = []; confi = [];
for trl = 1:config.numtrials.total
    
    tic % ding ding ding
    
    clearpict(1);
    clearpict(2);
    clearpict(3);
    clearpict(4);
    
    cd(paths.supports); loadpict('checker.bmp',1,0,0,bg_w,bg_h);loadpict('checker.bmp',2,0,0,bg_w,bg_h);
    loadpict('checker.bmp',3,0,0,bg_w,bg_h); loadpict('checker.bmp',4,0,0,bg_w,bg_h);
    
    setforecolour(.5,.5,.5);
    settextstyle('Arial', 300);
    preparestring('+', 1);  % Prepare buffer 1 (fixation cross)
    loadpict(stimvec{1,trl},2); % prepare buffer 2 (stim)
    loadpict('newold.bmp',3); % prepare buffer 3 (response)
    loadpict('memorytest_confidence.bmp',4);
    
    
    f1c=f1c+1; SOT_f1(f1c,1)=toc; tic  % reset stopwatch
    drawpict(1) % fixation x
    if eyeflag==1
        Eyelink('Message',[num2str(trl) 'f1'])  % send trialnumber and stimulustype to eyetracker
    end
    wait(config.timing.fixation1(trl));
    
    sc=sc+1; SOT_s(sc,1)=toc;
    drawpict(2) % stim
    if eyeflag==1
        Eyelink('Message',[num2str(trl) 's'])  % send trialnumber and stimulustype to eyetracker
    end
    wait(config.timing.cue);
    
    stimonset = drawpict(3); % response screen
    if eyeflag==1
        Eyelink('Message',[num2str(trl) 'r'])  % send trialnumber and stimulustype to eyetracker
    end
    %     waitkeydown(1) % for breezing
    waitkeydown(config.timing.resp,[OldKey,NewKey]);
    [key_ON, t_ON, n_ON] = getkeydown;
    
    drawpict(4)
    %     waitkeydown(1) % for breezing
    waitkeydown(config.timing.resp,[Surekey,Guesskey])
    [key_SG, t_SG, n_SG] = getkeydown;
    
    % record confidence
    if key_SG == Surekey
        confi = 1;
    elseif key_SG == Guesskey
        confi = 0;
    end
        
    % accuracy
    if key_ON == OldKey
        resp = 1;
        if stimvec{2,trl} == 1
            accu = 1;
        else
            accu = 0;
        end
    elseif key_ON == NewKey
        resp = 2;
        if stimvec{2,trl} == 2
            accu = 1;
        else
            accu = 0;
        end
    end
    
    % RT
    if n_ON == 0 %if they didn't press anything
        rt = NaN; % mark their reaction time as nan
    else % otherwise
        rt = t_ON(1) - stimonset; % and their reaction time is the time they pressed the button-the time the stimulus apprered
    end
    
    % record
    if n_ON == 0 || n_SG == 0
        results.keypress(trl,1:2)   = [NaN NaN];
        results.oldnew_resp(trl,1)  = NaN;
        results.accuracy(trl,1)     = NaN;
        results.confidence(trl,1)   = NaN;
        results.rt(trl,1)           = NaN;
    else
        results.keypress(trl,1:2)   = [key_ON key_SG];
        results.oldnew_resp(trl,1)  = resp;
        results.accuracy(trl,1)     = accu;
        results.confidence(trl,1)   = confi;
        results.rt(trl,1)           = rt;
    end
    results.trl(trl,:)          = stimvec(:,trl)';
    results.all = horzcat(results.oldnew_resp,results.accuracy,results.confidence,results.rt);
    
    % calculate cumulative SOT
    results.SOT.cumulative.fix  = cumsum(SOT_f1);
    results.SOT.cumulative.stim = cumsum(SOT_s);
    results.SOT.raw.fix  = SOT_f1;
    results.SOT.raw.stim = SOT_s;
    
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

%% wrap up

% calculate final cumulative SOT
results.SOT.cumulative.fix  = cumsum(SOT_f1);
results.SOT.cumulative.stim = cumsum(SOT_s);
results.SOT.raw.fix  = SOT_f1;
results.SOT.raw.stim = SOT_s;

% save data
dat.memorytest.labels = {'stimlist_all:' 'filenames' 'Indoor(1)/Outdoor(2)' 'Old(1)/New(2)' ' '; ...
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
% accuracies
dat.emomem.test.results.quickstat.MeanAccuracy = nanmean(results.accuracy);
dat.emomem.test.results.quickstat.MeanAccuracy_Old = ...
    nanmean(results.accuracy(cell2mat(stimlist(:,4))==1));
dat.emomem.test.results.quickstat.MeanAccuracy_New = ...
    nanmean(results.accuracy(cell2mat(stimlist(:,4))==2));
dat.emomem.test.results.quickstat.MeanAccuracy_Rew = ...
    nanmean(results.accuracy(cell2mat(stimlist(:,2))==rewcond));
dat.emomem.test.results.quickstat.MeanAccuracy_Pun = ...
    nanmean(results.accuracy(cell2mat(stimlist(:,2))==puncond));
% confidences
dat.emomem.test.results.quickstat.MeanConfidence_Old = ...
    nanmean(results.confidence(cell2mat(stimlist(:,4))==1));
dat.emomem.test.results.quickstat.MeanConfidence_New = ...
    nanmean(results.confidence(cell2mat(stimlist(:,4))==2));
dat.emomem.test.results.quickstat.MeanConfidence_Rew = ...
    nanmean(results.confidence(cell2mat(stimlist(:,2))==rewcond));
dat.emomem.test.results.quickstat.MeanConfidence_Pun = ...
    nanmean(results.confidence(cell2mat(stimlist(:,2))==puncond));
dat.emomem.test.results.quickstat.MedianConfidence_Old = ...
    nanmedian(results.confidence(cell2mat(stimlist(:,4))==1));
dat.emomem.test.results.quickstat.MedianConfidence_New = ...
    nanmedian(results.confidence(cell2mat(stimlist(:,4))==2));
dat.emomem.test.results.quickstat.MedianConfidence_Rew = ...
    nanmedian(results.confidence(cell2mat(stimlist(:,2))==rewcond));
dat.emomem.test.results.quickstat.MedianConfidence_Pun = ...
    nanmedian(results.confidence(cell2mat(stimlist(:,2))==puncond));
% FA rates
dat.emomem.test.results.quickstat.FA = ...
    nanmean(double(results.accuracy==0 & cell2mat(stimlist(:,4))==1));
dat.emomem.test.results.quickstat.FA_Rew = ...
    nanmean(double(results.accuracy==0 & ...
    cell2mat(stimlist(:,4))==1 & cell2mat(stimlist(:,2))==rewcond));
dat.emomem.test.results.quickstat.FA_Pun = ...
    nanmean(double(results.accuracy==0 & ...
    cell2mat(stimlist(:,4))==1 & cell2mat(stimlist(:,2))==puncond));

save([paths.behav num2str(ID) '_EM.mat'],'dat')

% terminate protocol
clearpict(1); cd(paths.supports)
loadpict('ending.png',1,0,0)
drawpict(1); wait(3000);
stop_cogent;


fprintf('\n%%%%%%%%%%%%%%%%%%%%%%\n COMPLETED  \n%%%%%%%%%%%%%%%%%%%%%%\n')
fprintf('\nMean Accuracy: %3.3f\n', nanmean(results.accuracy))




end
