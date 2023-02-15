%% REPLAY Retrieval 1 
%  In emergency contact Alex:
%              +49 162 713 20 49
%              yeo-jin.yi.15@ucl.ac.uk

%%  work log

%   24_09_2020      created the script

%% START

function [dat] = runREPLAYretieval_1_ptb(ID,AgeGroup,BlockNum,behavfilename,path_ECHO)
%% basic setups

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         MRIFlag = 1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set paths
paths.parent   = path_ECHO;
paths.supports = [paths.parent 'stim\supports\'];
paths.stim     = [paths.parent 'stim\enc\'];
paths.control  = [paths.parent 'stim\control\'];
paths.behav    = [paths.parent 'data\behav\'];
paths.logs     = [paths.parent 'data\logs\'];
paths.cogent   = [paths.parent 'Cogent2000v1.33\'];
paths.coggph   = [paths.parent 'CogGph\'];
addpath(genpath(paths.cogent))
addpath(genpath(paths.coggph))

% create the data structure
cd(paths.behav)
if exist(behavfilename) == 2
    load([paths.behav behavfilename])
else
    dat = [];
    dat.ID       = ID;
    dat.AgeGroup = AgeGroup;
    dat.retrieval1 = [];
end

% load presentation list
load([paths.parent 'stim\stiminfo\RET1_7T_' num2str(BlockNum) '.mat'])
load([paths.parent 'stim\stiminfo\RET1_7T_' num2str(BlockNum) '_answerKeys.mat'])

% scanner preparation
if MRIFlag == 1
    dat.retrieval1.mri          = [];
    dat.retrieval1.mri.scanPort = 1;
    dat.retrieval1.mri.dummy    = 5;        % no. of dummy vols
    dat.retrieval1.mri.nslice   = 51;       % no. of slices
    dat.retrieval1.mri.TE       = 32;       % time for each slice in msec
    dat.retrieval1.mri.TR       = 2.34;     % MRIinfo.nslice * MRIinfo.TE/1000; % time for volume in msec
end

%% experimental setups

% condition markers
cue1 = 1; cue2 = 2; control = 3; nulls = 0;

% timings
config.timing.fixation = eval(['cell2mat(RET17T_' num2str(BlockNum) '(:,5)).*1000']);
config.timing.cue      = 8000;
config.timing.response = 3000; % maximum waiting time for the button press
config.timing.feedback = 3000;
config.timing.null     = 3000;

% --- for breezeing through --- %
% % config.timing.fixation = ones(1,400);
% % config.timing.learning = 1;
% % config.timing.cue      = 1;
% % config.timing.response = 1; % maximum waiting time for the button press
% % config.timing.feedback = 1;
% % config.timing.null     = 1;
% ----------------------------- %

% trials
config.numtrials.total   = length(eval(['RET17T_' num2str(BlockNum)]));
config.numtrials.cue1    = numel(find(cell2mat(eval(['RET17T_' num2str(BlockNum) '(:,4)']))==cue1));
config.numtrials.cue2    = numel(find(cell2mat(eval(['RET17T_' num2str(BlockNum) '(:,4)']))==cue2));
config.numtrials.control = numel(find(cell2mat(eval(['RET17T_' num2str(BlockNum) '(:,4)']))==control));
config.numtrials.null    = numel(find(cell2mat(eval(['RET17T_' num2str(BlockNum) '(:,4)']))==nulls));


config.stim.stimvec  = cell2mat(eval(['RET17T_' num2str(BlockNum) '(:,4)']));
config.stim.stimlist = eval(['RET17T_' num2str(BlockNum) ]);
config.stim.filelist = eval(['RET17T_' num2str(BlockNum) '(:,1:3)']);

% make it easier to see within the script
stimvec = config.stim.stimvec;
rooms   = config.stim.filelist(:,1);
cues    = config.stim.filelist(:,2);
options = config.stim.filelist(:,3);

fprintf('\n\nTASK SETUP DONE')
fprintf('\nINITIALISING REPLAY 7T RETRIEVAL 1 TASK\n\n')

%% run

clear KbCheck;

% check for Opengl compatibility, abort otherwise:
AssertOpenGL;

% Make sure keyboard mapping is the same on all supported operating systems
% Apple MacOS/X, MS-Windows and GNU/Linux:
KbName('UnifyKeyNames');

RightKey   = KbName('3#'); %%%%% check again
MiddleKey  = KbName('2@');
LeftKey    = KbName('1!');
MRItrigger = KbName('t');

%%%%%%%%%%
% button instructions: left to right, 1 - 2 - 3 - 4 
%%%%%%%%%%%

% Get screenNumber of stimulation display. We choose the display with
% the maximum index, which is usually the right one, e.g., the external
% display on a Laptop:
screens=Screen('Screens'); % this should be one, which is the main screen that you're looking at
screenNumber=max(screens);

% Hide the mouse cursor:
HideCursor;

% background colour should be black
% Returns as default the mean black value of screen:
black=BlackIndex(screenNumber);

% Open a double buffered fullscreen window on the stimulation screen
% 'screenNumber' and choose/draw a black background. 'w' is the handle
% used to direct all drawing commands to that window - the "Name" of
% the window. 'wRect' is a rectangle defining the size of the window.
% See "help PsychRects" for help on such rectangles and useful helper
% functions:
%%%%% experimental screens must be 1280x1024!!! %%%%%
oldRes=SetResolution(0,1280,1024);
[w, wRect]=Screen('OpenWindow',screenNumber, black);
[mx, my] = RectCenter(wRect);
W=wRect(3); H=wRect(4);

% Set text size (Most Screen functions must be called after
% opening an onscreen window, as they only take window handles 'w' as
% input:
Screen('TextSize', w, 20);
Screen('TextFont', w, 'Arial');

% Do dummy calls to GetSecs, WaitSecs, KbCheck to make sure
% they are loaded and ready when we need them - without delays
% in the wrong moment:
KbCheck;
WaitSecs(0.1);
GetSecs;

if BlockNum == 1
    
% ----- print instructions ----- %

cd(paths.supports)

% wait a bit between trials
WaitSecs(0.500);

% initialize KbCheck and variables to make sure they're
% properly initialized/allocted by Matlab - this to avoid time
% delays in the critical reaction time measurement part of the
% script:
[KeyIsDown, endrt, KeyCode]=KbCheck;

% ---- instruction 1 ---- %
% read stimulus image into matlab matrix 'imdata':
clear imdata
imdata=imread(char('replay_instruction1_ret.bmp'));
% make texture image out of image matrix 'imdata'
tex=Screen('MakeTexture', w, imdata);
% Draw texture image to backbuffer. It will be automatically
% centered in the middle of the display if you don't specify a
% different destination:
Screen('DrawTexture', w, tex);
% Show stimulus on screen at next possible display refresh cycle,
% and record stimulus onset time in 'startrt':
[VBLTimestamp t0_inst FlipTimestamp]=Screen('Flip', w);
tmp=VBLTimestamp;
clear rsp
[rsp] = getKeys([LeftKey,MiddleKey,RightKey],Inf);
WaitSecs(0.2);

% ---- instruction 2 ---- %
clear imdata
imdata=imread(char('replay_instruction2_ret.bmp'));
tex=Screen('MakeTexture', w, imdata);
Screen('DrawTexture', w, tex);
[VBLTimestamp stimOnsetTime FlipTimestamp]=Screen('Flip', w);
[rsp] = getKeys([LeftKey,MiddleKey,RightKey],Inf);
WaitSecs(0.2);

% ---- instruction 3 ---- %
clear imdata
imdata=imread(char('replay_instruction3_ret.bmp'));
tex=Screen('MakeTexture', w, imdata);
Screen('DrawTexture', w, tex);
[VBLTimestamp stimOnsetTime FlipTimestamp]=Screen('Flip', w);
[rsp] = getKeys([LeftKey,MiddleKey,RightKey],Inf);
WaitSecs(0.2);


else
   
cd(paths.supports)

% ---- instruction reminder ---- %
clear imdata
imdata=imread(char('replay_instruction1_old.bmp'));
tex=Screen('MakeTexture', w, imdata);
Screen('DrawTexture', w, tex);
[VBLTimestamp t0_inst FlipTimestamp]=Screen('Flip', w);
[rsp] = getKeys([LeftKey,MiddleKey,RightKey],Inf);
WaitSecs(0.2);

end


% let the scanner start the task, allow n dummy volumes
if MRIFlag == 1
    
    % hail the operator
    cd(paths.supports)
    clear imdata
    imdata=imread(char('operator_hail.bmp'));
    tex=Screen('MakeTexture', w, imdata);
    Screen('DrawTexture', w, tex);
    [VBLTimestamp stimOnsetTime FlipTimestamp]=Screen('Flip', w);
    [rsp] = getKeys([KbName('space')],Inf);
    WaitSecs(0.2);
    
    % start the dummy scan
    cd(paths.supports)
    clear imdata
    imdata=imread(char('ready.bmp'));
    tex=Screen('MakeTexture', w, imdata);
    Screen('DrawTexture', w, tex);
    [VBLTimestamp t0_standby FlipTimestamp]=Screen('Flip', w);
    
    %%%%%%%%%%%%%%
    tic % ding ding ding
    %%%%%%%%%%%%%%
    
    [rsp] = getKeys([KbName('t')],Inf); 
    trig1=toc; %%%
    dat.retrieval1.mri.dummyscanfirst= rsp;
    disp('1')
    WaitSecs(0.1);
    for dum = 1:(dat.retrieval1.mri.dummy - 2)
        [rsp] = getKeys([KbName('t')],Inf);
        disp(num2str(dum+1))
        WaitSecs(0.1);
    end
    [rsp] = getKeys([KbName('t')],Inf);
    trig5=toc; %%%
    disp('5')
    dat.retrieval1.mri.dummyscanlast = rsp; % the last dummy scan
    WaitSecs(0.1);
    
    cd(paths.supports);
    clear imdata
    imdata=imread(char('fix.bmp'));
    tex=Screen('MakeTexture', w, imdata);
    Screen('DrawTexture', w, tex);
    [VBLTimestamp t0_fix0 FlipTimestamp]=Screen('Flip', w); t0_fix0_raw = toc;
    WaitSecs(dat.retrieval1.mri.TR);

    
else
    
    %%%%%%%%%%%%%%
    tic % ding ding ding
    %%%%%%%%%%%%%%
    
    cd(paths.supports);
    clear imdata
    imdata=imread(char('fix.bmp'));
    tex=Screen('MakeTexture', w, imdata);
    Screen('DrawTexture', w, tex);
    [VBLTimestamp t0_fix0 FlipTimestamp]=Screen('Flip', w); t0_fix0_raw = toc;
    WaitSecs(1);
    
end

results = []; results.SOT = [];
if MRIFlag == 1
results.SOT.cgflip.t0_fix0 = t0_fix0;
results.SOT.cgflip.t0_standby = t0_standby;
results.SOT.cgflip.t0_inst = t0_inst;
results.SOT.raw.trigfirst=trig1;
results.SOT.raw.triglast=trig5;
end
results.SOT.raw.t0_fix0 = t0_fix0_raw;
eventmarker = 0;
i = 0; fixc = 0; cuec = 0; respc = 0; ctrlc = 0; nullc = 0; fbc = 0; cnfc = 0; respc2=0;

SOT_f1=[]; SOT_r = []; SOT_control = []; SOT_cue = []; SOT_null = []; SOT_fb = []; SOT_cnf = [];
for trl = 1:length(stimvec)
    
    i = i+1;
    fprintf('\n============================================\n')
    fprintf('Trial %d\n', trl)
    
    % ------------------------ present fix x -------------------------- %
    
    cd(paths.supports)
    clear imdata
    imdata=imread(char('fix.bmp'));
    tex=Screen('MakeTexture', w, imdata);
    Screen('DrawTexture', w, tex);
    [VBLTimestamp stimOnsetTime FlipTimestamp]=Screen('Flip', w);
    fixc=fixc+1; SOT_f1(fixc,1)=toc;
    WaitSecs(config.timing.fixation(trl)/1000)
    eventmarker = eventmarker+1;
    results.presentation{eventmarker,1} = '+';
        
    % ----------------------------------------------------------------- %

    if stimvec(trl) == cue1 || stimvec(trl) == cue2
        
        % ------------------------- cue screen ------------------------ %
        
        categories(i) = stimvec(trl);
        cd(paths.supports);
        
        clear Qtext tex_Qtext tex_cue cuescreen
        
        % prepare screens
        Qtext=imread(char('question_cue.bmp'));
        
        % Calculate image position (center of the screen)
        imageSize = size(Qtext);
        pos = [(W-imageSize(2))/2 ...
            (H-imageSize(1))/2+110 ...
            (W+imageSize(2))/2 ...
            (H+imageSize(1))/2+110];        
        tex_Qtext=Screen('MakeTexture', w, Qtext);
        Screen('DrawTexture', w, tex_Qtext, [], pos); clear pos imageSize
        
        cd(paths.stim);
        cuename = cues{trl};
        cuescreen=imread(cuename);
        % Calculate image position (center of the screen)
        imageSize = size(cuescreen);
        pos = [(W-imageSize(2))/2 ...
            (H-imageSize(1))/2-70 ...
            (W+imageSize(2))/2 ...
            (H+imageSize(1))/2-70];        
        tex_cue=Screen('MakeTexture', w, cuescreen);
        Screen('DrawTexture', w, tex_cue, [], pos); clear pos imageSize
        
        [VBLTimestamp stimOnsetTime FlipTimestamp]=Screen('Flip', w);
        cuec = cuec+1; SOT_cue(cuec,1)=toc;
        
        WaitSecs(config.timing.cue/1000);
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = cuename;
        
        % ------------------------------------------------------------- %
        
        
        
        % ---------------------- response screen ---------------------- %
        
        %%%%%% load objects %%%%%%
        option1 = options{trl,1}{1,1}; 
        option2 = options{trl,1}{1,2}; 
        option3 = options{trl,1}{1,3};
        
        cd(paths.stim);
        
        % === option 1 === %
        opt1=imread(char(option1));
        % Calculate image position (center of the screen)
        imageSize = size(opt1);
        pos = [(W-imageSize(2))/2-300 ...
            (H-imageSize(1))/2-40 ...
            (W+imageSize(2))/2-300 ...
            (H+imageSize(1))/2-40];        
        tex_opt1=Screen('MakeTexture', w, opt1);
        Screen('DrawTexture', w, tex_opt1,[], pos); clear pos imageSize
        
        % === option 2 === %
        opt2=imread(char(option2));
        % Calculate image position (center of the screen)
        imageSize = size(opt2);
        pos = [(W-imageSize(2))/2 ...
            (H-imageSize(1))/2-40 ...
            (W+imageSize(2))/2 ...
            (H+imageSize(1))/2-40];        
        tex_opt2=Screen('MakeTexture', w, opt2);
        Screen('DrawTexture', w, tex_opt2,[], pos); clear pos imageSize
        
        % === option 3 === %
        opt3=imread(char(option3));
        % Calculate image position (center of the screen)
        imageSize = size(opt2);
        pos = [(W-imageSize(2))/2+300 ...
            (H-imageSize(1))/2-40 ...
            (W+imageSize(2))/2+300 ...
            (H+imageSize(1))/2-40];
        tex_opt3=Screen('MakeTexture', w, opt3);
        Screen('DrawTexture', w, tex_opt3,[], pos); clear pos imageSize
        
        
        [VBLTimestamp responset1 FlipTimestamp]=Screen('Flip', w); tmp = toc;
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = options{trl,1};
        [rsp] = getKeys([RightKey,MiddleKey,LeftKey],config.timing.response/1000);
        
        % RT
        if ~isstruct(rsp) %if they didn't press anything
            response1 = NaN; % mark their response as nan
            keypress1 = NaN;
            rt1 = NaN; % mark their reaction time as nan
            wait(0);
        else % otherwise
            response1 = KbName(rsp.keyName); % their response is whatever key they pressed.
            keypress1 = rsp.keyName;
            fprintf(['\nkey pressed: %d \n'],KbName(rsp.keyName))
            rt1 = rsp.RT; % and their reaction time is the time they pressed the button-the time the stimulus apprered
            respc=respc+1; SOT_r(respc,1) = tmp+rt1;
            clear tmp
            
            
            %%%%%% load answer screen %%%%%%
            %%%%%% load objects %%%%%%
            
            option1 = options{trl,1}{1,1};
            option2 = options{trl,1}{1,2};
            option3 = options{trl,1}{1,3};
            
            cd(paths.stim);
            
            % === option 1 === %
            opt1=imread(char(option1));
            % Calculate image position (center of the screen)
            imageSize = size(opt1);
            pos = [(W-imageSize(2))/2-300 ...
                (H-imageSize(1))/2-40 ...
                (W+imageSize(2))/2-300 ...
                (H+imageSize(1))/2-40];
            tex_opt1=Screen('MakeTexture', w, opt1);
            Screen('DrawTexture', w, tex_opt1,[], pos); clear pos imageSize
            
            % === option 2 === %
            opt2=imread(char(option2));
            % Calculate image position (center of the screen)
            imageSize = size(opt2);
            pos = [(W-imageSize(2))/2 ...
                (H-imageSize(1))/2-40 ...
                (W+imageSize(2))/2 ...
                (H+imageSize(1))/2-40];
            tex_opt2=Screen('MakeTexture', w, opt2);
            Screen('DrawTexture', w, tex_opt2,[], pos); clear pos imageSize
            
            % === option 3 === %
            opt3=imread(char(option3));
            % Calculate image position (center of the screen)
            imageSize = size(opt2);
            pos = [(W-imageSize(2))/2+300 ...
                (H-imageSize(1))/2-40 ...
                (W+imageSize(2))/2+300 ...
                (H+imageSize(1))/2-40];
            tex_opt3=Screen('MakeTexture', w, opt3);
            Screen('DrawTexture', w, tex_opt3,[], pos); clear pos imageSize
            
            cd(paths.supports);            
            if response1 == RightKey
                cursorpos = 300;
            elseif response1 == MiddleKey
                cursorpos = 0;
            elseif  response1 == LeftKey
                cursorpos = -300;
            end
            if ~isnan(response1)
            cursor=imread(char('cursor.bmp'));
            % Calculate image position (center of the screen)
            imageSize = size(cursor);
            pos = [(W-imageSize(2))/2+cursorpos ...
                (H-imageSize(1))/2+150 ...
                (W+imageSize(2))/2+cursorpos ...
                (H+imageSize(1))/2+150];
            tex_cursor=Screen('MakeTexture', w, cursor);
            Screen('DrawTexture', w, tex_cursor,[], pos); clear pos imageSize
            [VBLTimestamp responset1 FlipTimestamp]=Screen('Flip', w);
            end
            WaitSecs(config.timing.response/1000-rt1);
        end
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = 'response';
        
        % accuracy
        
        % ----------- correct responses ----------- %
        if (response1 == RightKey && Answers(trl)==3)
            accuracy = 1;
            fprintf('\nCorrect \n')
        elseif (response1 == MiddleKey && Answers(trl)==2)
            accuracy = 1;
            fprintf('\nCorrect \n')
        elseif (response1 == LeftKey && Answers(trl)==1)
            accuracy = 1;
            fprintf('\nCorrect \n')
            
        % -------- internal lure responses -------- %
        elseif (response1 == LeftKey && InternalLures(trl)==1 && ~isnan(response1)) ...
                || (response1 == MiddleKey && InternalLures(trl)==2 && ~isnan(response1)) ...
                || (response1 == RightKey && InternalLures(trl)==3 && ~isnan(response1))
            accuracy = -1;
            fprintf('\nIncorrect, internal lure \n')
            
        % -------- external lure responses -------- %
        elseif (response1 == LeftKey && ExternalLures(trl)==1 && ~isnan(response1)) ...
                || (response1 == MiddleKey && ExternalLures(trl)==2 && ~isnan(response1)) ...
                || (response1 == RightKey && ExternalLures(trl)==3 && ~isnan(response1))
            accuracy = 0;
            fprintf('\nIncorrect, external lure \n')
            
        % ----------- missed responses! ----------- %
        elseif isnan(response1)
            accuracy = NaN;
            fprintf('\nMissed \n')
            
        else
            accuracy = NaN;
            fprintf('\nwhat''s going on?! \n')

        end
        
        % ----------------------------------------------------------------- %
        
        
        
        % --------------------- confidence rating ------------------------- %
        
        %%%%%% load stimulus %%%%%%
        cd(paths.supports);
        confscreen=imread('confidence.bmp');
        imageSize = size(confscreen);
            pos = [(W-imageSize(2))/2 ...
                (H-imageSize(1))/2-10 ...
                (W+imageSize(2))/2 ...
                (H+imageSize(1))/2-10];
        tex_confscreen=Screen('MakeTexture', w, confscreen);
        Screen('DrawTexture', w, tex_confscreen,[], pos); clear pos imageSize
        
        %%%%%% draw stimulus screen %%%%%%
        cnfc = cnfc+1; SOT_cnf(cnfc,1)=toc;
        [VBLTimestamp responset2 FlipTimestamp]=Screen('Flip', w);
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = 'confidence';
        [rsp] = getKeys([RightKey,LeftKey],config.timing.response/1000);
        

        % RT
        if ~isstruct(rsp) %if they didn't press anything
            response2 = NaN; % mark their response as nan
            keypress2 = NaN;
            rt2 = NaN; % mark their reaction time as nan
        else % otherwise
            response2 = KbName(rsp.keyName); % their response is whatever key they pressed.
            keypress2 = rsp.keyName;
            fprintf(['\nkey pressed: %d \n'], KbName(rsp.keyName))
            rt2 = rsp.RT; % and their reaction time is the time they pressed the button-the time the stimulus apprered
            respc2=respc2+1; %SOT_r(respc2,1) = tmp+rt2;
            clear tmp
        end
        
        % record response
        if response2 == LeftKey
            confidence = 1;
            fprintf('\Sure \n')
            
        elseif response2 == RightKey
            confidence = 0;
            fprintf('\nNot sure \n')
            
        else
            confidence = NaN;
            fprintf('\nwhat''s going on?! \n')

        end
        WaitSecs(config.timing.response/1000-rt1);
        
        % ----------------------------------------------------------------- %
        
        
        
        % ----------------------- feedback screen ------------------------- %
        
        %%%%%% load stimulus %%%%%%
        roomname = rooms{trl};
        cd(paths.stim);
        roomscreen=imread(roomname);
        imageSize = size(roomscreen);
            pos = [(W-imageSize(2))/2 ...
                (H-imageSize(1))/2-70 ...
                (W+imageSize(2))/2 ...
                (H+imageSize(1))/2-70];
        tex_roomscreen=Screen('MakeTexture', w, roomscreen);
        Screen('DrawTexture', w, tex_roomscreen,[], pos); clear pos imageSize
                
        %%%%%% draw stimulus screen %%%%%%
        fbc=fbc+1; SOT_fb(fbc,1)=toc;
        [VBLTimestamp responset2 FlipTimestamp]=Screen('Flip', w);
        WaitSecs(config.timing.feedback/1000);
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = [roomname '_fb'];
        
        % ----------------------------------------------------------------- %
        
        % record
        
        try
            results.keypress(i,1) = keypress1; % response
            results.keypress_conf(i,1) = keypress2; % confidence
        catch
            try
                results.keypress(i,1) = keypress1(1); % response
                results.keypress_conf(i,1) = keypress2(1); % confidence
            catch
                results.keypress(i,1) = NaN; % response
                results.keypress_conf(i,1) = NaN; % confidence
            end
        end
        
        results.rt_resp(i,1)  = rt1;
        results.rt_conf(i,1)  = rt2;
        results.accuracy(i,1) = accuracy;
        results.confidence(i,1) = confidence;
        results.trl{i,1}      = roomname;
        results.trl{i,2}      = categories(i);
        
        
        
    elseif stimvec(trl) == control
        
        %%%%%% load stimulus %%%%%%
        cd(paths.supports);
        ctrlscreen=imread('question_control.bmp');
        imageSize = size(ctrlscreen);
        pos = [(W-imageSize(2))/2 ...
            (H-imageSize(1))/2+70 ...
            (W+imageSize(2))/2 ...
            (H+imageSize(1))/2+70];        
        tex_ctrlscreen=Screen('MakeTexture', w, ctrlscreen);
        Screen('DrawTexture', w, tex_ctrlscreen,[], pos); clear pos imageSize
       
        cd(paths.control); clear roomname
        roomname = rooms{trl};
        roomscreen=imread(roomname);
        imageSize = size(roomscreen);
        pos = [(W-imageSize(2))/2 ...
            (H-imageSize(1))/2-70 ...
            (W+imageSize(2))/2 ...
            (H+imageSize(1))/2-70];        
        tex_roomscreen=Screen('MakeTexture', w, roomscreen);
        Screen('DrawTexture', w, tex_roomscreen,[], pos); clear pos imageSize
        
        %%%%%% draw stimulus screen %%%%%%
        ctrlc = ctrlc+1; SOT_control(ctrlc,1)=toc;
        [VBLTimestamp responset2 FlipTimestamp]=Screen('Flip', w);
        WaitSecs(config.timing.cue/1000)
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = roomname;
        
        categories(i) = stimvec(trl);
        response1 = NaN; % mark their response as nan
        keypress1 = NaN;
        rt1 = NaN; % mark their reaction time as nan
        response2 = NaN;
        keypress2 = NaN;
        rt2 = NaN;
        
        % record
        results.keypress(i,1) = NaN; % response
        results.rt_resp(i,1)  = NaN;
        results.accuracy(i,1) = NaN;
        results.confidence(i,1) = NaN;
        results.trl{i,1}      = 'control';
        results.trl{i,2}      = 0;
        results.keypress_conf(i,1) = NaN;
        results.rt_conf(i,1)  = NaN;
        
    elseif  stimvec(trl) == nulls
        
        cd(paths.supports)
        cgsetsprite(0);
        cgloadbmp(1,'bg.bmp');
        cgloadbmp(2,'fix.bmp');
        cgdrawsprite(1,0,0);
        cgdrawsprite(2,0,40);
        nullc=nullc+1; SOT_null(nullc,1)=toc;
        cgflip(0,0,0);
        wait(config.timing.null)
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = 'Null';
        
        % record
        results.keypress(i,1) = NaN; % response
        results.rt_resp(i,1)  = NaN;
        results.accuracy(i,1) = NaN;
        results.confidence(i,1) = NaN;
        results.trl{i,1}      = 'null';
        results.trl{i,2}      = 0;
        results.keypress_conf(i,1) = NaN;
        results.rt_conf(i,1)  = NaN;
        
    end
    
    % calculate cumulative SOT
    results.SOT.raw.fix  = SOT_f1;
    results.SOT.raw.resp = SOT_r;
    results.SOT.raw.cue  = SOT_cue;
    results.SOT.raw.control= SOT_control;
    results.SOT.raw.null = SOT_null;
    results.SOT.raw.feedback = SOT_fb;
    results.SOT.raw.confidence = SOT_cnf;
    
    % intermediate save
    dat.retrieval1.results = results;
    save([paths.behav num2str(ID) '_ret1' num2str(BlockNum) '_backup.mat'],'dat')
    
    fprintf('\n============================================\n')  
    
end


% wrap up
% terminate protocol
cd(paths.supports)
if BlockNum == 1
    clear imdata
    imdata=imread(char('instruction_ende1.bmp'));
    tex=Screen('MakeTexture', w, imdata);
    Screen('DrawTexture', w, tex);
elseif BlockNum == 2
    clear imdata
    imdata=imread(char('instruction_ende2.bmp'));
    tex=Screen('MakeTexture', w, imdata);
    Screen('DrawTexture', w, tex);
elseif BlockNum == 3
    clear imdata
    imdata=imread(char('instruction_ende3.bmp'));
    tex=Screen('MakeTexture', w, imdata);
    Screen('DrawTexture', w, tex);
elseif BlockNum == 4
    clear imdata
    imdata=imread(char('instruction_ende4.bmp'));
    tex=Screen('MakeTexture', w, imdata);
    Screen('DrawTexture', w, tex);
end
[VBLTimestamp t_TaskEnd FlipTimestamp]=Screen('Flip', w);
t_TaskEnd_raw = toc;
WaitSecs(3);
sca; ShowCursor;
    
% calculate final cumulative SOT
results.SOT.raw.fix  = SOT_f1;
results.SOT.raw.resp = SOT_r;
results.SOT.raw.cue  = SOT_cue;
results.SOT.raw.control= SOT_control;
results.SOT.raw.null = SOT_null;
results.SOT.raw.feedback = SOT_fb;
results.SOT.raw.confidence = SOT_cnf;
results.SOT.raw.end  = t_TaskEnd_raw;


% save results
dat.retrieval1.results = results;
dat.retrieval1.labels  = {'categories-> 1=cue1, 2=cue2, 3=control, 0=null'; ...
    'results.trl: accuracy-> 1=correct, -1=internal, 0=external, NaN=missed)'};
dat.retrieval1.config  = config;
dat.retrieval1.config.keymap = KbName('KeyNames');
save([paths.behav behavfilename],'dat')

disp('*******************************************')
disp('******* end of the retrieval task 1 *******')


end
