%% ENGRAMS experiment: recognition

function [dat] = engrams_recognition(id,block,condition,phase,mri,pupil,taskpath)

sca % close any open Psychtoolbox windows

% temp settings
path_ECHO      = taskpath;
ID             = id;
Block          = block;

if strcmp(condition,'Emotional')==1
    conditions = 'emo';
elseif strcmp(condition,'Neutral')==1
    conditions = 'neu';
end

if strcmp(phase,'Original')==1
    phases = 'orig';
elseif strcmp(phase,'Recombination')==1
    phases = 'recombi';
end

if strcmp(mri,'Yes')==1
    MRIFlag = 1;
elseif strcmp(mri,'No')==1
    MRIFlag = 0;
end

if strcmp(pupil,'Yes')==1
    EyeFlag = 1;
elseif strcmp(pupil,'No')==1
    EyeFlag = 0;
end


behavfilename  = [num2str(ID) '_rcg_' phases '_' num2str(Block) '.mat'];

% set paths
paths.parent   = path_ECHO;
paths.stimuli  = [paths.parent conditions '/'];
paths.stimlist = [paths.parent conditions '/Excel/'];
paths.obj      = [paths.parent conditions '/images/objects/'];
paths.scenes   = [paths.parent conditions '/images/scenes/'];
paths.elements = [paths.parent 'elements/'];
paths.bg       = paths.parent;
paths.behav    = [paths.parent 'data/behav/'];
paths.logs     = [paths.parent 'data/logs/'];

addpath(genpath([paths.parent 'scripts']))

if exist(strcat(paths.logs, num2str(ID), '_rcg_', phases, '_', num2str(Block), '_log_', string(date), '.txt')) == 2
K=string(datestr(now));
try
diary(strcat(paths.logs, num2str(ID), '_rcg_', phases, '_', num2str(Block), '_log_', string(date),'_',K{1}(end-1:end),'.txt'))
catch
diary(strcat(paths.logs, num2str(ID), '_rcg_', phases, '_', num2str(Block), '_log_', string(date),'_',K(end-1:end),'.txt'))
end
else
diary(strcat(paths.logs, num2str(ID), '_rcg_', phases, '_', num2str(Block), '_log_', string(date), '.txt'))
end

disp(strcat('ENGRAMS RECOGNITION ', string(datetime)))

% create the data structure
cd(paths.behav)
if exist(behavfilename) == 2
    load([paths.behav behavfilename])
else
    dat = [];
    dat.ID  = ID;
    dat.rcg = [];
    dat.rcg.condition = conditions;
    dat.rcg.phase     = phases;
    dat.rcg.measured_on = datetime;
end

% scanner preparation
if MRIFlag == 1
    dat.rcg.mri          = [];
    dat.rcg.mri.scanPort = 1;
    dat.rcg.mri.dummy    = 5;        % no. of dummy vols
    dat.rcg.mri.nslice   = 51;       % no. of slices
    dat.rcg.mri.TE       = 32;       % time for each slice in msec
    dat.rcg.mri.TR       = 2.34;     % MRIinfo.nslice * MRIinfo.TE/1000; % time for volume in msec
end

disp('environment prepared')
disp(['subject ID ' num2str(ID) ', block number ' num2str(Block) ', condition ' conditions])

%% prepare stimuli

config = [];

% ======== prepare stimuli ======== %

rng('shuffle'); % prepare for true randomisation

stimlist = [];
stimlist = readtable([paths.stimlist 'recognition_' phases '_' conditions '.xlsx']);

config.stimuli.numtrials = size(stimlist,1);

v = 1:config.stimuli.numtrials;
randomIndex = randperm(length(v));
randomizedVector = v(randomIndex);

config.stimuli.stimlist.scene       = stimlist(randomizedVector,1);
config.stimuli.stimlist.object      = stimlist(randomizedVector,2);
config.stimuli.stimlist.options     = stimlist(randomizedVector,3:5);
config.stimuli.stimlist.answerkeys  = table2array(stimlist(randomizedVector,6:8));
config.stimuli.stimlist.all         = stimlist(randomizedVector,:);

% ================================= %


% ======== timing information ======== %

rng('shuffle'); % shuffle again
config.timing.fixation = randi([10, 40], 60, 1) * 100;%eval(['cell2mat(RET17T_' num2str(BlockNum) '(:,5)).*1000']);
config.timing.intermission = [randi([10, 15], 60, 1) * 100, randi([5, 10], 60, 1) * 100];%eval(['cell2mat(RET17T_' num2str(BlockNum) '(:,5)).*1000']);
config.timing.scene    = 2000;
config.timing.object   = 2000;
config.timing.together = 2000;
config.timing.confidenceQ = 2000;

% ================================= %

disp('task setup done')
disp('initialising...')

%% run


clear KbCheck;

% check for Opengl compatibility, abort otherwise:
AssertOpenGL;
Screen('Preference', 'SkipSyncTests', 1);

% Make sure keyboard mapping is the same on all supported operating systems
% Apple MacOS/X, MS-Windows and GNU/Linux:
KbName('UnifyKeyNames');

%%%%%%%%%%
% button instructions: left to right, a(left thumb) - b(left index) - c(right index) - d(right thumb) 
%%%%%%%%%%%

FarRightKey      = KbName('d'); %%%%% check again
MiddleRightKey   = KbName('c');
MiddleLeftKey    = KbName('b');
FarLeftKey       = KbName('a');
MRItrigger       = KbName('s');


% Get screenNumber of stimulation display. We choose the display with
% the maximum index, which is usually the right one, e.g., the external
% display on a Laptop:
screens=Screen('Screens'); % this should be one, which is the main screen that you're looking at
screenNumber=max(screens);

HideCursor; % Hide the mouse cursor:

% background colour should be black
% Returns as default the mean black value of screen:
black=BlackIndex(screenNumber);

% Open a double buffered fullscreen window on the stimulation screen
% 'screenNumber' and choose/draw a black background. 'w' is the handle
% used to direct all drawing commands to that window - the "Name" of
% the window. 'wRect' is a rectangle defining the size of the window.
% See "help PsychRects" for help on such rectangles and useful helper
% functions:
% oldRes=SetResolution(0,2560,1440);
oldRes=SetResolution(0,1920,1080); % at 7T, only this resolution works
[w, wRect]=Screen('OpenWindow',screenNumber, black);
[mx, my] = RectCenter(wRect);
W=wRect(3); H=wRect(4);

% Store the current font size and style
oldSize = Screen('TextSize', w);
oldStyle = Screen('TextStyle', w);

Screen('TextFont',w,'Arial');
Screen('TextSize',w,17);

% Get the size of the screen
[screenXpixels, screenYpixels] = Screen('WindowSize', w);

% Load the background
bg = imread([paths.bg 'background.jpg']);

% Create a texture from the image
bgtexture = Screen('MakeTexture', w, bg);

% Get the size of the image
[bgHeight, bgWidth, ~] = size(bg);

% Calculate the number of tiles needed to cover the screen
tilesX = ceil(screenXpixels / bgWidth);
tilesY = ceil(screenYpixels / bgHeight);

% Tile the bg image across the screen
for i = 0:(tilesX-1)
    for j = 0:(tilesY-1)
        destRect = [i * bgWidth, j * bgHeight, (i+1) * bgWidth, (j+1) * bgHeight];
        Screen('DrawTexture', w, bgtexture, [], destRect);
    end
end

% Text presentation doesn't work at the connectome scanner, so we opt for
% presenting an image

% Load the instruction text
instruction = imread([paths.elements 'instr_' phases '_recognition.png']);

% Create a texture from the image
tex=Screen('MakeTexture', w, instruction);
Screen('DrawTexture', w, tex);

% Show stimulus on screen at next possible display refresh cycle,
% and record stimulus onset time in 'startrt':
[VBLTimestamp, t0_inst, FlipTimestamp]=Screen('Flip', w);
tmp=VBLTimestamp;
clear rsp
[rsp] = getKeys([FarLeftKey,MiddleLeftKey,MiddleRightKey,FarRightKey],Inf);

% Let the scanner start the task, allow n dummy volumes
if MRIFlag == 1
    
    % hail the operator
    
    % load operator information screen
    clear tex
    opinfo = imread([paths.elements 'operator_hail.bmp']);
    % Create a texture from the image
    tex=Screen('MakeTexture', w, opinfo);
    Screen('DrawTexture', w, tex);
    
    [VBLTimestamp, ~, FlipTimestamp]=Screen('Flip', w);
    
    [rsp] = getKeys(KbName('space'),Inf);
    
    WaitSecs(0.2);
    
    % start the dummy scan
    % Tile the bg image across the screen
    for i = 0:(tilesX-1)
        for j = 0:(tilesY-1)
            destRect = [i * bgWidth, j * bgHeight, (i+1) * bgWidth, (j+1) * bgHeight];
            Screen('DrawTexture', w, bgtexture, [], destRect);
        end
    end
    
    clear tex
    standby = imread([paths.elements 'standby.bmp']);
    % Create a texture from the image
    tex=Screen('MakeTexture', w, standby);
    Screen('DrawTexture', w, tex);

    [VBLTimestamp, t0_standby, FlipTimestamp]=Screen('Flip', w);
    
    %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%
         tic % ding ding ding
    %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%
   
    
    % Initialize counters and flags
    numTriggers = 0;
    trig1 = NaN;
    trig5 = NaN;
    
    % Loop until 5 triggers are detected
    while numTriggers < dat.rcg.mri.dummy
        % Check the state of the keyboard
        [keyIsDown, ~, keyCode] = KbCheck;
        
        % If a key is pressed
        if keyIsDown
            if keyCode(MRItrigger)
                % Increment the trigger counter
                numTriggers = numTriggers + 1;
                
                % Record the timing for the 1st and 5th triggers
                if numTriggers == 1
                    trig1 = toc;
                elseif numTriggers == 5
                    trig5 = toc;
                    break; % Exit the loop after the 5th trigger
                end
            end
            
            % Wait until all keys are released to avoid detecting the same press again
            while KbCheck; end
        end
        
        % Optional: insert a brief pause to reduce CPU load
        WaitSecs(0.01);
    end
    
%     [rsp] = getKeys(MRItrigger,Inf); 
%     trig1=toc; %%%
%     dat.rcg.mri.dummyscanfirst= rsp;
%     disp('5')
%     WaitSecs(0.1);
%     for dum = 1:(dat.rcg.mri.dummy - 2)
%         [rsp] = getKeys(MRItrigger,Inf);
%         disp(num2str(4-dum+1))
%         WaitSecs(0.1);
%     end
%     [rsp] = getKeys(MRItrigger,Inf);
%     trig5=toc; %%%
%     disp('1')
%     dat.rcg.mri.dummyscanlast = rsp; % the last dummy scan
    
    
    % draw the first fixation
    % Tile the bg image across the screen
    for i = 0:(tilesX-1)
        for j = 0:(tilesY-1)
            destRect = [i * bgWidth, j * bgHeight, (i+1) * bgWidth, (j+1) * bgHeight];
            Screen('DrawTexture', w, bgtexture, [], destRect);
        end
    end
    
    % Define fixation cross parameters
    crossSize = 50; % Size of the cross in pixels
    crossColour = [200, 200, 200]; % Grey ish, easy on the eyes
    lineWidth = 10; % Adjust the line width as needed
    % Define the coordinates for the cross lines
    crossLines = [-crossSize, crossSize, 0, 0; 0, 0, -crossSize, crossSize];
    % Draw the fixation cross
    Screen('DrawLines', w, crossLines, lineWidth, crossColour, [mx, my]);
        
    [VBLTimestamp, t0_fix0, FlipTimestamp]=Screen('Flip', w); t0_fix0_raw = toc;
    
    WaitSecs(dat.rcg.mri.TR);

    
else
    
    %%%%%%%%%%%%%%
    tic % ding ding ding
    %%%%%%%%%%%%%%
    
    % draw the first fixation
    % Tile the bg image across the screen
    for i = 0:(tilesX-1)
        for j = 0:(tilesY-1)
            destRect = [i * bgWidth, j * bgHeight, (i+1) * bgWidth, (j+1) * bgHeight];
            Screen('DrawTexture', w, bgtexture, [], destRect);
        end
    end
    
    % Define fixation cross parameters
    crossSize = 50; % Size of the cross in pixels
    crossColour = [200, 200, 200]; % Grey ish, easy on the eyes
    lineWidth = 10; % Adjust the line width as needed
    % Define the coordinates for the cross lines
    crossLines = [-crossSize, crossSize, 0, 0; 0, 0, -crossSize, crossSize];
    % Draw the fixation cross
    Screen('DrawLines', w, crossLines, lineWidth, crossColour, [mx, my]);
    
    [VBLTimestamp, t0_fix0, FlipTimestamp]=Screen('Flip', w); t0_fix0_raw = toc;
    
    WaitSecs(1);
    
end


% start the trials

results = []; %results.SOT = [];
% if MRIFlag == 1
% results.SOT.ptb.t0_fix0 = t0_fix0;
% results.SOT.ptb.t0_standby = t0_standby;
% results.SOT.ptb.t0_inst = t0_inst;
% results.SOT.raw.trigfirst=trig1;
% results.SOT.raw.triglast=trig5;
% end
% results.SOT.raw.t0_fix0 = t0_fix0_raw;
eventmarker = 0;
i = 0; fixc = 0; scenec = 0; objc = 0; togc = 0; rcgQc = 0; confQc = 0; resp1c = 0; resp2c = 0;
SOT_f1=[]; SOT_together = []; SOT_obj = []; SOT_scene = []; SOT_confQ =[]; SOT_rcgQ =[]; SOT_resprcgQ = []; SOT_respconfQ = [];
rcgQ=[]; confQ=[]; valQ_objs=[]; RT_rcgQ=[]; RT_obj=[];
OptionChosen_rcgQ=[]; accuracies=[]; accuracies_bin=[]; confidences=[];
dat.rcg.results = [];

for trl = 1%:length(config.stimuli.numtrials)
    
    i = i+1;
    fprintf('\n============================================\n')
    fprintf('Trial %d\n', trl)
    
    % ------------------------ present fix x -------------------------- %
    
    for i = 0:(tilesX-1)
        for j = 0:(tilesY-1)
            destRect = [i * bgWidth, j * bgHeight, (i+1) * bgWidth, (j+1) * bgHeight];
            Screen('DrawTexture', w, bgtexture, [], destRect);
        end
    end
    
    % Define fixation cross parameters
    crossSize = 50; % Size of the cross in pixels
    crossColour = [200, 200, 200]; % Grey ish, easy on the eyes
    lineWidth = 10; % Adjust the line width as needed
    % Define the coordinates for the cross lines
    crossLines = [-crossSize, crossSize, 0, 0; 0, 0, -crossSize, crossSize];
    % Draw the fixation cross
    Screen('DrawLines', w, crossLines, lineWidth, crossColour, [mx, my]);

    [VBLTimestamp, stimOnsetTime, FlipTimestamp]=Screen('Flip', w);
    fixc=fixc+1; SOT_f1(fixc,1)=toc;
    
    WaitSecs(config.timing.fixation(trl)/1000)
    eventmarker = eventmarker+1;
    results.presentation{eventmarker,1} = '+';
        
    % ----------------------------------------------------------------- %
    
    
    
    % ------------------------ present scene alone -------------------------- %
        
    % Load Scene
    
    for i = 0:(tilesX-1)
        for j = 0:(tilesY-1)
            destRect = [i * bgWidth, j * bgHeight, (i+1) * bgWidth, (j+1) * bgHeight];
            Screen('DrawTexture', w, bgtexture, [], destRect);
        end
    end
    
    % Define fixation cross parameters
    crossSize = 50; % Size of the cross in pixels
    crossColour = [200, 200, 200]; % Grey ish, easy on the eyes
    lineWidth = 10; % Adjust the line width as needed
    % Define the coordinates for the cross lines
    crossLines = [-crossSize, crossSize, 0, 0; 0, 0, -crossSize, crossSize];
    % Draw the fixation cross
    Screen('DrawLines', w, crossLines, lineWidth, crossColour, [mx, my]);

    SceneStim = imread([paths.stimuli config.stimuli.stimlist.scene{trl,1}{1}]);
    SceneTexture = Screen('MakeTexture', w, SceneStim);
    
    % Define size and position of the image
    [s1, s2, ~] = size(SceneStim);
    aspectRatio = s2 / s1;
    imageHeight = screenYpixels * 0.4; 
    imageWidth = imageHeight * aspectRatio;
    
    % Adjust the vertical offset as neededaa
    verticalOffset = screenYpixels * 0.25; % Adjust this value as needed
    dstRect = CenterRectOnPointd([0 0 imageWidth imageHeight], mx, verticalOffset);

    
    % Draw the image
    Screen('DrawTexture', w, SceneTexture, [], dstRect);
    
    [VBLTimestamp, stimOnsetTime, FlipTimestamp]=Screen('Flip', w);
    scenec=scenec+1; SOT_scene(scenec,1)=toc;
    WaitSecs(config.timing.scene/1000);
    eventmarker = eventmarker+1;
    results.presentation{eventmarker,1} = config.stimuli.stimlist.scene{trl,1}{1,1};
    
    % ----------------------------------------------------------------- %
    
    
    % ------------------------ present intermission x -------------------------- %
    
    for i = 0:(tilesX-1)
        for j = 0:(tilesY-1)
            destRect = [i * bgWidth, j * bgHeight, (i+1) * bgWidth, (j+1) * bgHeight];
            Screen('DrawTexture', w, bgtexture, [], destRect);
        end
    end
    
    % Define fixation cross parameters
    crossSize = 50; % Size of the cross in pixels
    crossColour = [200, 200, 200]; % Grey ish, easy on the eyes
    lineWidth = 10; % Adjust the line width as needed
    % Define the coordinates for the cross lines
    crossLines = [-crossSize, crossSize, 0, 0; 0, 0, -crossSize, crossSize];
    % Draw the fixation cross
    Screen('DrawLines', w, crossLines, lineWidth, crossColour, [mx, my]);
    
    [VBLTimestamp, stimOnsetTime, FlipTimestamp]=Screen('Flip', w);
    
    WaitSecs(config.timing.intermission(trl,1)/1000)
        
    % ----------------------------------------------------------------- %
    
    
    
    % -------------------- present options (obj) ----------------------- %
    
    for i = 0:(tilesX-1)
        for j = 0:(tilesY-1)
            destRect = [i * bgWidth, j * bgHeight, (i+1) * bgWidth, (j+1) * bgHeight];
            Screen('DrawTexture', w, bgtexture, [], destRect);
        end
    end
    
    % Define fixation cross parameters
    crossSize = 50; % Size of the cross in pixels
    crossColour = [200, 200, 200]; % Grey ish, easy on the eyes
    lineWidth = 10; % Adjust the line width as needed
    % Define the coordinates for the cross lines
    crossLines = [-crossSize, crossSize, 0, 0; 0, 0, -crossSize, crossSize];
    % Draw the fixation cross
    Screen('DrawLines', w, crossLines, lineWidth, crossColour, [mx, my]);
   
    
    % Draw the question
    img_rcgQ       = imread([paths.elements 'question_recognition.png']);
    
    img_objL    = imread([paths.stimuli config.stimuli.stimlist.options{trl,1}{1}]);
    img_objM    = imread([paths.stimuli config.stimuli.stimlist.options{trl,2}{1}]);
    img_objR    = imread([paths.stimuli config.stimuli.stimlist.options{trl,3}{1}]);
    
    rcgQtex     = Screen('MakeTexture', w, img_rcgQ);
    
    objLtex     = Screen('MakeTexture', w, img_objL);
    objMtex     = Screen('MakeTexture', w, img_objM);
    objRtex     = Screen('MakeTexture', w, img_objR);
    
    % Define the starting x position for the question iamge
    questionDestRectUp = [screenXpixels / 2 - size(img_rcgQ, 2) / 2, screenYpixels * 0.325 - size(img_rcgQ, 1) / 2, screenXpixels / 2 + size(img_rcgQ, 2) / 2, screenYpixels * 0.325 + size(img_rcgQ, 1) / 2];
    Screen('DrawTexture', w, rcgQtex, [], questionDestRectUp);
    
    % Define spacing between the option images
    optionSpacing = 100; % Adjust this value to increase or decrease the space between images
    
    % Calculate the total width of the option images including spacing
    optionTextures = [objLtex, objMtex, objRtex];
    totalImageWidth = sum(arrayfun(@(x) size(Screen('GetImage', x), 2), optionTextures));
    totalWidthWithSpacing = totalImageWidth + optionSpacing * (length(optionTextures) - 1);
    
    % Calculate the starting X position for the first option image
    startX = ((screenXpixels - totalWidthWithSpacing) / 2) + 50;
    
    % Draw the option images with spacing
    for i = 1:length(optionTextures)
        [imgH, imgW, ~] = size(Screen('GetImage', optionTextures(i)));
        optionDestRect = [startX, screenYpixels * 0.7 - (imgH*0.8) / 2, startX + (imgW*0.8), screenYpixels * 0.7 + (imgH*0.8) / 2];
        Screen('DrawTexture', w, optionTextures(i), [], optionDestRect);
        startX = startX + imgW + optionSpacing; % Move to the next position with spacing
    end
    
    [VBLTimestamp, stimOnsetTime, FlipTimestamp]=Screen('Flip', w); tmp = toc;
    rcgQc=rcgQc+1; SOT_rcgQ(rcgQc,1)=toc;
    eventmarker = eventmarker+1;
    results.presentation{eventmarker,1} = 'recognitionQ';
    
    
    [rsp] = getKeys([FarLeftKey,MiddleLeftKey,MiddleRightKey,FarRightKey],Inf); tmp = toc;
    
    
    % RT
    if ~isstruct(rsp) %if they didn't press anything
        response1 = NaN; % mark their response as nan
        keypress1 = NaN;
        rt1 = NaN; % mark their reaction time as nan
        wait(0);
        rcgQ(trl,1) = NaN;
        RT_rcgQ(trl,1) = NaN;
        
    else % otherwise
        response1 = KbName(rsp.keyName); % their response is whatever key they pressed.
        keypress1 = rsp.keyName;
        fprintf(['\nkey pressed: %d \n'],KbName(rsp.keyName))
        
        if response1==FarLeftKey
            rcgQ(trl,1) = 1;
        elseif response1==MiddleLeftKey
            rcgQ(trl,1) = 2;
        elseif response1==MiddleRightKey
            rcgQ(trl,1) = 2;
        elseif response1==FarRightKey
            rcgQ(trl,1) = 3;
        else
            warning('unknown key pressed')
            KbName(rsp.keyName);
        end            
        
        rt1 = rsp.RT; % and their reaction time is the time they pressed the button-the time the stimulus apprered
        RT_rcgQ(trl,1) = rt1;
        resp1c=resp1c+1; SOT_resprcgQ(resp1c,1) = tmp+rt1;
        clear tmp
    end
    
    
    % calculate accuracy (hit=1, internallure=-1, externallure=0)
    
    OptionChosen_rcgQ(trl,1) = config.stimuli.stimlist.answerkeys(trl,rcgQ(trl,1));
    if OptionChosen_rcgQ(trl,1)==1
        accuracies(trl,1)=1;
        accuracies_bin(trl,1)=1;
    elseif OptionChosen_rcgQ(trl,1)==-1
        accuracies(trl,1)=-1;
        accuracies_bin(trl,1)=0;
    elseif OptionChosen_rcgQ(trl,1)==0
        accuracies(trl,1)=0;
        accuracies_bin(trl,1)=0;
    else
        accuracies(trl,1)=NaN;
        accuracies_bin(trl,1)=NaN;
    end
        
    % ----------------------------------------------------------------- %

    
     % ------------------------ present intermission x -------------------------- %
    
    for i = 0:(tilesX-1)
        for j = 0:(tilesY-1)
            destRect = [i * bgWidth, j * bgHeight, (i+1) * bgWidth, (j+1) * bgHeight];
            Screen('DrawTexture', w, bgtexture, [], destRect);
        end
    end
    
    % Define fixation cross parameters
    crossSize = 50; % Size of the cross in pixels
    crossColour = [200, 200, 200]; % Grey ish, easy on the eyes
    lineWidth = 10; % Adjust the line width as needed
    % Define the coordinates for the cross lines
    crossLines = [-crossSize, crossSize, 0, 0; 0, 0, -crossSize, crossSize];
    % Draw the fixation cross
    Screen('DrawLines', w, crossLines, lineWidth, crossColour, [mx, my]);
    
    [VBLTimestamp, stimOnsetTime, FlipTimestamp]=Screen('Flip', w);
    
    WaitSecs(0.25)
        
    % ----------------------------------------------------------------- %
    
    
    
    % ------------------ present confidence rating -------------------- %

    for i = 0:(tilesX-1)
        for j = 0:(tilesY-1)
            destRect = [i * bgWidth, j * bgHeight, (i+1) * bgWidth, (j+1) * bgHeight];
            Screen('DrawTexture', w, bgtexture, [], destRect);
        end
    end
    
    % Define fixation cross parameters
    crossSize = 50; % Size of the cross in pixels
    crossColour = [200, 200, 200]; % Grey ish, easy on the eyes
    lineWidth = 10; % Adjust the line width as needed
    % Define the coordinates for the cross lines
    crossLines = [-crossSize, crossSize, 0, 0; 0, 0, -crossSize, crossSize];
    % Draw the fixation cross
    Screen('DrawLines', w, crossLines, lineWidth, crossColour, [mx, my]);
   
    
    % Draw the question
    img_confQ       = imread([paths.elements 'rating_confidence.png']);
    
    img_garnicht    = imread([paths.elements 'gar_nicht.png']);
    img_etwas       = imread([paths.elements 'etwas.png']);
    img_ziemlich    = imread([paths.elements 'ziemlich.png']);
    img_sehr        = imread([paths.elements 'sehr.png']);
    
    rcgQtex         = Screen('MakeTexture', w, img_confQ);
    garnichttex     = Screen('MakeTexture', w, img_garnicht);
    etwastex        = Screen('MakeTexture', w, img_etwas);
    ziemlichtex     = Screen('MakeTexture', w, img_ziemlich);
    sehrtex         = Screen('MakeTexture', w, img_sehr);
    
    % Define the starting x position for the question iamge
    questionDestRectUp = [screenXpixels / 2 - size(img_confQ, 2) / 2, screenYpixels * 0.325 - size(img_confQ, 1) / 2, screenXpixels / 2 + size(img_confQ, 2) / 2, screenYpixels * 0.325 + size(img_confQ, 1) / 2];
    Screen('DrawTexture', w, rcgQtex, [], questionDestRectUp);
    
    % Define spacing between the option images
    optionSpacing = 50; % Adjust this value to increase or decrease the space between images
    
    % Calculate the total width of the option images including spacing
    optionTextures = [garnichttex, etwastex, ziemlichtex, sehrtex];
    totalImageWidth = sum(arrayfun(@(x) size(Screen('GetImage', x), 2), optionTextures));
    totalWidthWithSpacing = totalImageWidth + optionSpacing * (length(optionTextures) - 1);
    
    % Calculate the starting X position for the first option image
    startX = (screenXpixels - totalWidthWithSpacing) / 2;
    
    % Draw the option images with spacing
    for i = 1:length(optionTextures)
        [imgH, imgW, ~] = size(Screen('GetImage', optionTextures(i)));
        optionDestRect = [startX, screenYpixels * 0.65 - imgH / 2, startX + imgW, screenYpixels * 0.65 + imgH / 2];
        Screen('DrawTexture', w, optionTextures(i), [], optionDestRect);
        startX = startX + imgW + optionSpacing; % Move to the next position with spacing
    end
    
    [VBLTimestamp, stimOnsetTime, FlipTimestamp]=Screen('Flip', w); tmp = toc;
    confQc=confQc+1; SOT_confQ(confQc,1)=toc;
    eventmarker = eventmarker+1;
    results.presentation{eventmarker,1} = 'confidenceRating';
    
    
    [rsp] = getKeys([FarLeftKey,MiddleLeftKey,MiddleRightKey,FarRightKey],Inf); tmp = toc;
    
    
    % RT
    if ~isstruct(rsp) %if they didn't press anything
        response2 = NaN; % mark their response as nan
        keypress2 = NaN;
        rt2 = NaN; % mark their reaction time as nan
        wait(0);
        confQ(trl,1) = NaN;
        RT_confQ(trl,1) = NaN;
        
    else % otherwise
        response2 = KbName(rsp.keyName); % their response is whatever key they pressed.
        keypress2 = rsp.keyName;
        fprintf(['\nkey pressed: %d \n'],KbName(rsp.keyName))
        
        if response2==FarLeftKey
            confQ(trl,1) = 0;
        elseif response2==MiddleLeftKey
            confQ(trl,1) = 1;
        elseif response2==MiddleRightKey
            confQ(trl,1) = 2;
        elseif response2==FarRightKey
            confQ(trl,1) = 3;
        else
            warning('unknown key pressed')
            KbName(rsp.keyName);
        end            
        
        rt2 = rsp.RT; % and their reaction time is the time they pressed the button-the time the stimulus apprered
        RT_confQ(trl,1) = rt2;
        resp2c=resp2c+1; SOT_respconfQ(resp2c,1) = tmp+rt2;
        clear tmp
    end
    
    
    % record confidence ratings
    confidences(trl,1)=confQ(trl,1);
    
    % ----------------------------------------------------------------- %
    
    
    
    % ------------------- present intermission x ---------------------- %
   
    for i = 0:(tilesX-1)
        for j = 0:(tilesY-1)
            destRect = [i * bgWidth, j * bgHeight, (i+1) * bgWidth, (j+1) * bgHeight];
            Screen('DrawTexture', w, bgtexture, [], destRect);
        end
    end
    
    % Define fixation cross parameters
    crossSize = 50; % Size of the cross in pixels
    crossColour = [200, 200, 200]; % Grey ish, easy on the eyes
    lineWidth = 10; % Adjust the line width as needed
    % Define the coordinates for the cross lines
    crossLines = [-crossSize, crossSize, 0, 0; 0, 0, -crossSize, crossSize];
    % Draw the fixation cross
    Screen('DrawLines', w, crossLines, lineWidth, crossColour, [mx, my]);

    [VBLTimestamp, stimOnsetTime, FlipTimestamp]=Screen('Flip', w);
    
    WaitSecs(config.timing.intermission(trl,2)/1000)
        
    % ----------------------------------------------------------------- %
    
    
    % ------------------------ present together -------------------------- %

    for i = 0:(tilesX-1)
        for j = 0:(tilesY-1)
            destRect = [i * bgWidth, j * bgHeight, (i+1) * bgWidth, (j+1) * bgHeight];
            Screen('DrawTexture', w, bgtexture, [], destRect);
        end
    end
    
    % Define fixation cross parameters
    crossSize = 50; % Size of the cross in pixels
    crossColour = [200, 200, 200]; % Grey ish, easy on the eyes
    lineWidth = 10; % Adjust the line width as needed
    % Define the coordinates for the cross lines
    crossLines = [-crossSize, crossSize, 0, 0; 0, 0, -crossSize, crossSize];
    % Draw the fixation cross
    Screen('DrawLines', w, crossLines, lineWidth, crossColour, [mx, my]);
    
    clear SceneStim
    SceneStim = imread([paths.stimuli config.stimuli.stimlist.scene{trl,1}{1,1}]);
    SceneTexture = Screen('MakeTexture', w, SceneStim);
    
    % Define size and position of the image
    [s1, s2, ~] = size(SceneStim);
    aspectRatio = s2 / s1;
    imageHeight = screenYpixels * 0.4; 
    imageWidth = imageHeight * aspectRatio;
    
    % Adjust the vertical offset as needed
    verticalOffset = screenYpixels * 0.25; % Adjust this value as needed
    dstRect = CenterRectOnPointd([0 0 imageWidth imageHeight], mx, verticalOffset);

    
    % Draw the image
    Screen('DrawTexture', w, SceneTexture, [], dstRect);
    
    clear ObjectStim
    ObjectStim = imread([paths.stimuli config.stimuli.stimlist.object{trl,1}{1,1}]);
    ObjTexture = Screen('MakeTexture', w, ObjectStim);
    
    % Define size and position of the image
    [s1, s2, ~] = size(ObjectStim);
    aspectRatio = s2 / s1;
    imageHeight = screenYpixels * 0.35; 
    imageWidth = imageHeight * aspectRatio;
    
    % Adjust the vertical offset as needed
    verticalOffset = screenYpixels * 0.725;
    horizontalOffset = screenXpixels * 0.625;
    
    dstRect = CenterRectOnPointd([0 0 imageWidth imageHeight], horizontalOffset, verticalOffset);
    
    % Draw the image
    Screen('DrawTexture', w, ObjTexture, [], dstRect);
    
    
    [VBLTimestamp, stimOnsetTime, FlipTimestamp]=Screen('Flip', w);
    togc=togc+1; SOT_together(togc,1)=toc;
    WaitSecs(config.timing.together/1000);
    eventmarker = eventmarker+1;
    results.presentation{eventmarker,1} = [config.stimuli.stimlist.scene{trl,1}{1,1} '_' config.stimuli.stimlist.object{trl,1}{1,1}];
    
    % ----------------------------------------------------------------- %

    % mid-save the things
    
    if MRIFlag==1
    dat.rcg.results.SOT.t0_standby_PTB = t0_standby;
    dat.rcg.results.SOT.t0_instruction_PTB = t0_inst;
    dat.rcg.results.SOT.trig1 = trig1;
    dat.rcg.results.SOT.trig5 = trig5;
    dat.rcg.results.SOT.t0_fix0 = t0_fix0_raw;
    dat.rcg.results.SOT.fixationX = SOT_f1;
    dat.rcg.results.SOT.scene = SOT_scene;
    dat.rcg.results.SOT.together = SOT_together;
    dat.rcg.results.SOT.RecognitionQ = SOT_rcgQ;
    dat.rcg.results.SOT.resp_RecognitionQ = SOT_resprcgQ;
    dat.rcg.results.SOT.ConfidenceQ = SOT_confQ;
    dat.rcg.results.SOT.resp_confidenceQ = SOT_respconfQ;
    
    dat.rcg.results.RT.RecognitionQ = RT_rcgQ;
    dat.rcg.results.RT.ConfidenceQ = RT_confQ;
    
    dat.rcg.results.accuracy = accuracies;
    dat.rcg.results.confidence = confidences;
    
    dat.rcg.results.presentationOrder = results.presentation;
    
    dat.rcg.results.codingInfo = 'confidence: 0=garnicht,1=etwas,2=ziemlich,3=sehr / recognition: 1=hit, -1=internalLure, 0=externalLure';
    
    dat.rcg.config.keymap = KbName('KeyNames');
    dat.rcg.config.stimlist = config.stimuli.stimlist;
    
    save([paths.behav 'tmp_' behavfilename],'dat')
    
    else
    dat.rcg.results.RT.RecognitionQ = RT_rcgQ;
    dat.rcg.results.RT.ConfidenceQ = RT_confQ;
    
    dat.rcg.results.accuracy = accuracies;
    dat.rcg.results.confidence = confidences;
    
    dat.rcg.results.presentationOrder = results.presentation;
    
    dat.rcg.results.codingInfo = 'confidence: 0=garnicht,1=etwas,2=ziemlich,3=sehr / recognition: 1=hit, -1=internalLure, 0=externalLure';
    
    dat.rcg.config.keymap = KbName('KeyNames');
    dat.rcg.config.stimlist = config.stimuli.stimlist;
    
    save([paths.behav 'tmp_' behavfilename],'dat')
    
    end

end

clear ending
ending=imread([paths.elements 'enc_end.png']);
tex=Screen('MakeTexture', w, ending);
Screen('DrawTexture', w, tex);

% Draw the text inside the grey box
% DrawFormattedText(w, finaltext, 'center', 'center', [200 200 200], [], [], [], 1.5, [], centeredRect);
[VBLTimestamp, t0_inst, FlipTimestamp]=Screen('Flip', w); t_TaskEnd_raw = toc;
WaitSecs(3);
sca; ShowCursor;

if MRIFlag==1
    dat.rcg.results.SOT.t0_standby_PTB = t0_standby;
    dat.rcg.results.SOT.t0_instruction_PTB = t0_inst;
    dat.rcg.results.SOT.trig1 = trig1;
    dat.rcg.results.SOT.trig5 = trig5;
    dat.rcg.results.SOT.t0_fix0 = t0_fix0_raw;
    dat.rcg.results.SOT.fixationX = SOT_f1;
    dat.rcg.results.SOT.scene = SOT_scene;
    dat.rcg.results.SOT.together = SOT_together;
    dat.rcg.results.SOT.RecognitionQ = SOT_rcgQ;
    dat.rcg.results.SOT.resp_RecognitionQ = SOT_resprcgQ;
    dat.rcg.results.SOT.ConfidenceQ = SOT_confQ;
    dat.rcg.results.SOT.resp_confidenceQ = SOT_respconfQ;
    
    dat.rcg.results.RT.RecognitionQ = RT_rcgQ;
    dat.rcg.results.RT.ConfidenceQ = RT_confQ;
    
    dat.rcg.results.accuracy = accuracies;
    dat.rcg.results.confidence = confidences;
    
    dat.rcg.results.presentationOrder = results.presentation;
    
    dat.rcg.results.codingInfo = 'confidence: 0=garnicht,1=etwas,2=ziemlich,3=sehr / recognition: 1=hit, -1=internalLure, 0=externalLure';
    
    dat.rcg.config.keymap = KbName('KeyNames');
    dat.rcg.config.stimlist = config.stimuli.stimlist;
    
    save([paths.behav behavfilename],'dat')
    
else
    dat.rcg.results.RT.RecognitionQ = RT_rcgQ;
    dat.rcg.results.RT.ConfidenceQ = RT_confQ;
    
    dat.rcg.results.accuracy = accuracies;
    dat.rcg.results.confidence = confidences;
    
    dat.rcg.results.presentationOrder = results.presentation;
    
    dat.rcg.results.codingInfo = 'confidence: 0=garnicht,1=etwas,2=ziemlich,3=sehr / recognition: 1=hit, -1=internalLure, 0=externalLure';
    
    dat.rcg.config.keymap = KbName('KeyNames');
    dat.rcg.config.stimlist = config.stimuli.stimlist;
    
    save([paths.behav behavfilename],'dat')
    
end



% % Close textures
% Screen('Close', texture);

% Close the screen
Screen('CloseAll');
disp('*********************************************')
disp('******* the end of the encoding phase *******')

disp(datetime)

diary off;

end
