%% REPLAY Retrieval 1 
%  In emergency contact Alex:
%              +49 162 713 20 49
%              yeo-jin.yi.15@ucl.ac.uk

%%  work log

%   24_09_2020      created the script

%% START

function [dat] = runREPLAYretieval_1_cg(ID,AgeGroup,BlockNum,behavfilename,path_ECHO)
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

% configs
config_display(1, 5, [0 0 0], [1 1 1], 'Arial', 20, 10, 0); % 1024x768
% config_display(screenMode, screenRes, bg ,fontcol, fontName, fontSize, number_of_buffers,0);   % open graphics window
config_sound;
config_keyboard;

% % configs
% config_display(1, 5, [0 0 0], [1 1 1], 'Arial', 20, 10, 0); % 1024x768
% % config_display(screenMode, screenRes, bg ,fontcol, fontName, fontSize, number_of_buffers,0);   % open graphics window
% config_sound;
% config_keyboard;


% scanner input
if MRIFlag == 1
    config_serial(dat.retrieval1.mri.scanPort); % links to scanner for waitslice
end

monitor = get(0,'screensize');
bg_w = monitor(3); bg_h = monitor(4);

%%%%%%%%%%%%
% start_cogent
%%%%%%%%%%%%

cgopen(5,0,0,1)
start_cogent

map        = getkeymap; % define keyboard IDs
RightKey   = 28; %%%%% check again
MiddleKey  = 29;
LeftKey    = 30;


%%%%%%%%%%
% button instructions: left to right, 1 - 2 - 3 - 4 
%%%%%%%%%%%



if BlockNum == 1
% print instructions
cd(paths.supports)
cgsetsprite(0);
cgloadbmp(1,'bg.bmp');
cgloadbmp(2,'replay_instruction1_ret.bmp');
cgdrawsprite(1,0,0);
cgdrawsprite(2,0,0);
cgflip(0,0,0); waitkeydown(inf,[RightKey,MiddleKey,LeftKey]); clearkeys

cgsetsprite(0);
cgloadbmp(1,'bg.bmp');
cgloadbmp(2,'replay_instruction2_ret.bmp');
cgdrawsprite(1,0,0);
cgdrawsprite(2,0,0);
cgflip(0,0,0); waitkeydown(inf,[RightKey,MiddleKey,LeftKey]); clearkeys

cgsetsprite(0);
cgloadbmp(1,'bg.bmp');
cgloadbmp(2,'replay_instruction2_confi7T.bmp');
cgdrawsprite(1,0,0);
cgdrawsprite(2,0,0);
cgflip(0,0,0); waitkeydown(inf,[RightKey,MiddleKey,LeftKey]); clearkeys

cgsetsprite(0);
cgloadbmp(1,'bg.bmp');
cgloadbmp(2,'replay_instruction3_ret.bmp');
cgdrawsprite(1,0,0);
cgdrawsprite(2,0,40);
t0_inst=cgflip(0,0,0); waitkeydown(inf,[RightKey,MiddleKey,LeftKey]); clearkeys

else
% print instructions
cd(paths.supports)
% cgsetsprite(0);
cgloadbmp(1,'bg.bmp');
cgloadbmp(2,'replay_instruction1_old.bmp');
cgdrawsprite(1,0,0);
cgdrawsprite(2,0,0);
t0_inst=cgflip(0,0,0); waitkeydown(inf,[RightKey,MiddleKey,LeftKey]); clearkeys

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
    
    [k, dat.retrieval1.mri.dummyscanfirst] = waitkeydown(inf,scannerinput); % the first dummy scan
    trig1=toc;
    for dum = 1:(dat.retrieval1.mri.dummy - 2)
        waitkeydown(inf,scannerinput);
    end
    [k, dat.retrieval1.mri.dummyscanlast] = waitkeydown(inf,scannerinput); % the last dummy scan
    trig5=toc;
    
%     %%%%%%%%%%%%%%
%     tic % ding ding ding
%     %%%%%%%%%%%%%%    
    
    cd(paths.supports);
    cgloadbmp(1,'bg.bmp');
    cgloadbmp(2,'fix.bmp');
    cgdrawsprite(1,0,0);
    cgdrawsprite(2,0,40);
    t0_fix0 = cgflip(0,0,0);   t0_fix0_raw = toc;    % load fixation cross
    wait(dat.retrieval1.mri.TR*1000)
    clearpict(1);
    
else
    
    %%%%%%%%%%%%%%
    tic % ding ding ding
    %%%%%%%%%%%%%%
    
    cd(paths.supports);
    cgloadbmp(1,'bg.bmp');
    cgloadbmp(2,'fix.bmp');
    cgdrawsprite(1,0,0);
    cgdrawsprite(2,0,40);
    t0_fix0 = cgflip(0,0,0);   t0_fix0_raw = toc;    % load fixation cross
    wait(1000)
    
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
    
    clearkeys;    

    % ------------------------ present fix x -------------------------- %
    
    cd(paths.supports)
    cgsetsprite(0);
    cgloadbmp(1,'bg.bmp');
    cgloadbmp(2,'fix.bmp');
    cgdrawsprite(1,0,0);
    cgdrawsprite(2,0,40);
    fixc=fixc+1; SOT_f1(fixc,1)=toc;
    cgflip(0,0,0);
    wait(config.timing.fixation(trl))
    eventmarker = eventmarker+1;
    results.presentation{eventmarker,1} = '+';
        
    % ----------------------------------------------------------------- %

    if stimvec(trl) == cue1 || stimvec(trl) == cue2
        
        % ------------------------- cue screen ------------------------ %
        
        categories(i) = stimvec(trl);
        cd(paths.supports);
        cgsetsprite(0);
        cgloadbmp(1,'bg.bmp',2000,2000);
        cgloadbmp(2,'question_cue.bmp');
        cd(paths.stim);
        cuename = cues{trl};
        cgloadbmp(3,cuename);
        cgdrawsprite(1,0,0);
        cgdrawsprite(2,0,-110);
        cgdrawsprite(3,0,70);
        cuec = cuec+1; SOT_cue(cuec,1)=toc;
        cgflip(0,0,0);
        wait(config.timing.cue);
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = cuename;
        
        % ------------------------------------------------------------- %
        
        
        
        % ---------------------- response screen ---------------------- %
        
        %%%%%% load objects %%%%%%
        option1 = options{trl,1}{1,1}; 
        option2 = options{trl,1}{1,2}; 
        option3 = options{trl,1}{1,3};
        cd(paths.supports);
        cgloadbmp(1,'bg.bmp');
        cd(paths.stim);
        cgloadbmp(2,option1);
        cgloadbmp(3,option2);
        cgloadbmp(4,option3);

        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgdrawsprite(2,300,40);
        cgdrawsprite(3,0,40);
        cgdrawsprite(4,-300,40);
        tmp = toc;
        responset1 = cgflip(0,0,0);
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = options{trl,1};        
        waitkeydown(config.timing.response,[RightKey,MiddleKey,LeftKey]);
        
        % RT
        [keypress1, t1, n1] = getkeydown;
        if n1 == 0 %if they didn't press anything
            response1 = NaN; % mark their response as nan
            keypress1 = NaN;
            rt1 = NaN; % mark their reaction time as nan
            wait(0);
        else % otherwise
            response1 = keypress1(1); % their response is whatever key they pressed.
            fprintf(['\nkey pressed: %d \n'],keypress1(1))
            rt1 = t1(1)/1000 - responset1; rt1 = rt1*1000; % and their reaction time is the time they pressed the button-the time the stimulus apprered
            respc=respc+1; SOT_r(respc,1) = tmp+rt1;
            clear tmp
            
            %%%%%% load answer screen %%%%%%
            %%%%%% load objects %%%%%%
            option1 = options{trl,1}{1,1};
            option2 = options{trl,1}{1,2};
            option3 = options{trl,1}{1,3};
            cd(paths.supports);
            cgloadbmp(1,'bg.bmp');
            cd(paths.stim);
            cgloadbmp(2,option1);
            cgloadbmp(3,option2);
            cgloadbmp(4,option3)
            cd(paths.supports);
            cgloadbmp(5,'cursor.bmp')
            cgsetsprite(0);
            cgdrawsprite(1,0,0);
            cgdrawsprite(2,300,40);
            cgdrawsprite(3,0,40);
            cgdrawsprite(4,-300,40);
            cgsetsprite(0);
            if response1 == RightKey
                cursorpos = -300;
            elseif response1 == MiddleKey
                cursorpos = 0;
            elseif  response1 == LeftKey
                cursorpos = 300;
            end
            if ~isnan(response1)
            cgdrawsprite(5,cursorpos,-150);
            end
            cgflip(0,0,0);
            
            wait(2000-rt1)
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
        cgloadbmp(1,'bg.bmp');
        cgloadbmp(2,'confidence.bmp');
        
        
        %%%%%% draw stimulus screen %%%%%%
        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgdrawsprite(2,0,10);
        cnfc = cnfc+1; SOT_cnf(cnfc,1)=toc;
        responset2 = cgflip(0,0,0);
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = 'confidence';
        waitkeydown(config.timing.response,[RightKey,LeftKey]);

        % RT
        [keypress2, t2, n2] = getkeydown;
        if n2 == 0 %if they didn't press anything
            response2 = NaN; % mark their response as nan
            keypress2 = NaN;
            rt2 = NaN; % mark their reaction time as nan
        else % otherwise
            response2 = keypress2(1); % their response is whatever key they pressed.
            fprintf(['\nkey pressed: %d \n'],keypress2(1))
            rt2 = t2(1)/1000 - responset2; rt2 = rt2*1000; % and their reaction time is the time they pressed the button-the time the stimulus apprered
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
        
        % ----------------------------------------------------------------- %
        
        
        
        % ----------------------- feedback screen ------------------------- %
        
        %%%%%% load stimulus %%%%%%
        roomname = rooms{trl};
        cd(paths.supports);
        cgloadbmp(1,'bg.bmp');
        cd(paths.stim);
        cgloadbmp(2,roomname,800,600);
        
        %%%%%% draw stimulus screen %%%%%%
        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgdrawsprite(2,0,70);
        fbc=fbc+1; SOT_fb(fbc,1)=toc;
        cgflip(0,0,0);
        wait(config.timing.feedback);
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
        cgloadbmp(1,'bg.bmp');
        cgloadbmp(2,'question_control.bmp');
        cd(paths.control); clear roomname
        roomname = rooms{trl};
        cgloadbmp(3,roomname);
        
        
        %%%%%% draw stimulus screen %%%%%%
        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgdrawsprite(2,0,-70);
        cgdrawsprite(3,0,70);
        ctrlc = ctrlc+1; SOT_control(ctrlc,1)=toc;
        cgflip(0,0,0);
        wait(config.timing.cue)
        
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


%% wrap up
% terminate protocol
cd(paths.supports)
if BlockNum == 1
cgloadbmp(1,'instruction_ende1.bmp');
elseif BlockNum == 2
cgloadbmp(1,'instruction_ende2.bmp');
elseif BlockNum == 3
cgloadbmp(1,'instruction_ende3.bmp');
elseif BlockNum == 4
cgloadbmp(1,'instruction_ende4.bmp');
end

cgsetsprite(0);
cgdrawsprite(1,0,0); 
t_TaskEnd = cgflip(0,0,0);   t_TaskEnd_raw = toc; 
wait(3000);
cgshut; stop_cogent;

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
dat.retrieval1.config.keymap = map;
save([paths.behav behavfilename],'dat')

cgshut
stop_cogent

disp('*******************************************')
disp('******* end of the retrieval task 1 *******')


end
