%% REPLAY ENCODING 
%  In emergency contact Alex:
%              +49 162 713 20 49
%              yeo-jin.yi.15@ucl.ac.uk

%%  work log

%   23_09_2020      created the script

%% START

function [dat] = runREPLAYencoding_cg(ID,AgeGroup,BlockNum,behavfilename,path_ECHO)


%% basic setups

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
load([paths.parent 'stim\stiminfo\ENC7T_' num2str(BlockNum) '.mat'])
load([paths.parent 'stim\stiminfo\ENC7T_' num2str(BlockNum) '_answerKeys.mat'])


%% experimental setups

% condition markers
cue1 = 1; cue2 = 2; control = 0;

% timings
dotjit = repmat(500:100:1000,1,8);
fixjit = repmat(500:100:2000,1,8);
config.timing.fixation = fixjit;
config.timing.dot      = dotjit;
config.timing.learning = 6000;
config.timing.cue      = 3000;
config.timing.response = 2000; % maximum waiting time for the button press
config.timing.feedback = 3000;
config.timing.null     = 3000;


% --- for breezeing through --- %
% config.timing.fixation = fixjit;
% config.timing.learning = 1;
% config.timing.cue      = 1;
% config.timing.response = 1; % maximum waiting time for the button press
% config.timing.feedback = 1;
% ----------------------------- %

% trials
config.numtrials.total   = length(eval(['ENC7T_' num2str(BlockNum)]));
config.numtrials.cue1    = numel(find(cell2mat(eval(['ENC7T_' num2str(BlockNum) '(:,4)']))==cue1));
config.numtrials.cue2    = numel(find(cell2mat(eval(['ENC7T_' num2str(BlockNum) '(:,4)']))==cue2));
config.numtrials.control = numel(find(cell2mat(eval(['ENC7T_' num2str(BlockNum) '(:,4)']))==control));

config.stim.stimvec  = cell2mat(eval(['ENC7T_' num2str(BlockNum) '(:,4)']));
config.stim.stimlist = eval(['ENC7T_' num2str(BlockNum) ]);
config.stim.filelist = eval(['ENC7T_' num2str(BlockNum) '(:,1:3)']);

% make it easier to see within the script
stimvec = config.stim.stimvec;
rooms   = config.stim.filelist(:,1);
cues    = config.stim.filelist(:,2);
options = config.stim.filelist(:,3);


fprintf('\n\nTASK SETUP DONE')
fprintf('\nINITIALISING REPLAY 7T ENCODING TASK\n\n')


%% run

% configs
config_display(0, 6, [0 0 0], [1 1 1], 'Arial', 20, 10, 0); % 1024x768
% config_display(screenMode, screenRes, bg ,fontcol, fontName, fontSize, number_of_buffers,0);   % open graphics window
config_sound;
config_keyboard;

monitor = get(0,'screensize');
bg_w = monitor(3); bg_h = monitor(4);

%%%%%%%%%%%%
start_cogent
%%%%%%%%%%%%

cgopen(6,0,0,0)

map     = getkeymap; % define keyboard IDs
RightKey   = 28; %%%%% check again
MiddleKey  = 29;
LeftKey    = 30;


% print instructions
cd(paths.supports)
% cgsetsprite(0);
cgloadbmp(1,'bg.bmp');
cgloadbmp(2,'replay_instruction1_enc.bmp');
cgdrawsprite(1,0,0);
cgdrawsprite(2,0,40);
cgflip(0,0,0); waitkeydown(inf,[RightKey,MiddleKey,LeftKey]); clearkeys
cgloadbmp(1,'bg.bmp');
cgloadbmp(2,'replay_instruction2_enc.bmp');
cgdrawsprite(1,0,0);
cgdrawsprite(2,0,40);
cgflip(0,0,0); waitkeydown(inf,[RightKey,MiddleKey,LeftKey]); clearkeys
cgloadbmp(1,'bg.bmp');
cgloadbmp(2,'replay_instruction3_enc.bmp');
cgdrawsprite(1,0,0);
cgdrawsprite(2,0,40);
cgflip(0,0,0); waitkeydown(inf,[RightKey,MiddleKey,LeftKey]); clearkeys
cgloadbmp(1,'bg.bmp');
cgloadbmp(2,'replay_instruction4_enc.bmp');
cgdrawsprite(1,0,0);
cgdrawsprite(2,0,40);
t0_inst=cgflip(0,0,0); waitkeydown(inf,[RightKey,MiddleKey,LeftKey]); clearkeys

    
%%%%%%%%%%%%%%
tic % ding ding ding
%%%%%%%%%%%%%%

cd(paths.supports);
cgloadbmp(1,'bg.bmp');
cgloadbmp(2,'fix.bmp');
cgdrawsprite(1,0,0);
cgdrawsprite(2,0,40);

cgflip(0,0,0);
wait(1000)

results = []; results.SOT = [];
results.SOT.cgflip.t0_inst = t0_inst;
eventmarker = 0;
i = 0; fix1c = 0; stimc = 0; respc1 = 0; nc = 0; fb = 0;
SOT_f1=[]; SOT_learn=[]; SOT_r=[]; SOT_control = []; SOT_cue =[];
for trl = 1:length(stimvec)
        
    i = i+1;
    fprintf('\n============================================\n')
    fprintf('Trial %d\n', trl)
    
    clearkeys;
    cd(paths.supports);
    
    
    % ------------------------ present fix x -------------------------- %
    
    cd(paths.supports)
    % cgsetsprite(0);
    cgloadbmp(1,'bg.bmp');
    cgloadbmp(2,'fix.bmp');
    cgdrawsprite(1,0,0);
    cgdrawsprite(2,0,40);
    fix1c=fix1c+1; SOT_f1(fix1c,1)=toc;
    cgflip(0,0,0);
    wait(config.timing.fixation(trl))
    eventmarker = eventmarker+1;
    results.presentation{eventmarker,1} = '+';
        
    % ----------------------------------------------------------------- %
    
    if stimvec(trl) ~= 0
        
        % --------------------- learning screen ----------------------- %

        %%%%%% load stimulus %%%%%%
        clear roomname
        roomname = rooms{trl};
        cd(paths.supports);
        cgloadbmp(1,'bg.bmp');
        cd(paths.stim);
        cgloadbmp(2,roomname,800,600);
        categories(1,i) = stimvec(trl);
        
        %%%%%% draw stimulus screen %%%%%%
%         cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgdrawsprite(2,0,70);
        stimc=stimc+1; SOT_learn(stimc,1)=toc;
        cgflip(0,0,0);
%         SOT_cg_s(sc) = cgflip(0,0,0);
        wait(config.timing.learning);
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = roomname;
        
        % ------------------------------------------------------------- %
        
        
        
        % ------------------------ dot screen ------------------------- %
        
        cd(paths.supports)
        % cgsetsprite(0);
        cgloadbmp(1,'bg.bmp');
        cgloadbmp(2,'dot.bmp');
        cgdrawsprite(1,0,0);
        cgdrawsprite(2,0,40);
        cgflip(0,0,0);
        wait(config.timing.dot(trl))
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = '.';
        
        % ------------------------------------------------------------- %
        
        
        
        % ------------------------- cue screen ------------------------ %
        
        %%%%%% load cue %%%%%%
        cd(paths.supports);
        cgloadbmp(1,'bg.bmp',2000,2000);
        cgloadbmp(2,'question_cue.bmp');
        cd(paths.stim);
        cuename = cues{trl};
        cgloadbmp(3,cuename);
%         cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgdrawsprite(2,0,-110);
        cgdrawsprite(3,0,70);
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

%         cgsetsprite(0);
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
            respc1=respc1+1; SOT_r(respc1,1) = tmp+rt1;
            clear tmp
            
            %%%%%% load answer screen %%%%%%
            option1 = options{trl,1}{1,1};
            option2 = options{trl,1}{1,2};
            option3 = options{trl,1}{1,3};
            cd(paths.supports);
            cgloadbmp(1,'bg.bmp');
            cgloadbmp(5,'cursor.bmp')
            cd(paths.stim);
            cgloadbmp(2,option1);
            cgloadbmp(3,option2);
            cgloadbmp(4,option3);
            %         cgsetsprite(0);
            if response1 == RightKey
                cursorpos = -300;
            elseif response1 == MiddleKey
                cursorpos = 0;
            elseif  response1 == LeftKey
                cursorpos = 300;
            end
            cgdrawsprite(1,0,0);
            cgdrawsprite(2,300,40);
            cgdrawsprite(3,0,40);
            cgdrawsprite(4,-300,40);
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
        
        
        
        
        
        % ----------------------- feedback screen ------------------------- %
        
        %%%%%% load stimulus %%%%%%
        cd(paths.supports);
        cgloadbmp(1,'bg.bmp');
        cd(paths.stim);
        cgloadbmp(2,roomname,800,600);
        
        %%%%%% draw stimulus screen %%%%%%
        cgsetsprite(0);
        cgdrawsprite(1,0,0);
        cgdrawsprite(2,0,70);
        fb=fb+1; SOT_fb(fb,1)=toc;
        cgflip(0,0,0);
        wait(config.timing.feedback);
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = [roomname '_fb'];
        
        
        % ----------------------------------------------------------------- %
        
        
    elseif stimvec(trl) == 0 % if the stim type is control
        
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
        nc = nc+1; SOT_control(nc,1)=toc;
        cgflip(0,0,0);
        wait(config.timing.learning)
        
        eventmarker = eventmarker+1;
        results.presentation{eventmarker,1} = roomname;
        
        categories(1,i) = stimvec(trl);
        response1 = NaN; % mark their response as nan
        keypress1 = NaN;
        rt1 = NaN; % mark their reaction time as nan
        
    end
        
    
    if stimvec(trl) == 0 % if the stim type is control
        results.keypress(i,1) = NaN; % response
        results.rt_resp(i,1)  = NaN;
        results.accuracy(i,1) = NaN;
        results.trl{i,1}      = 'control';
        results.trl{i,2}      = 0;
    else
        % record
        
        try
            results.keypress(i,1) = keypress1; % response
        catch
            try
                results.keypress(i,1) = keypress1(1); % response
            catch
                results.keypress(i,1) = NaN; % response
            end
        end
        
        results.rt_resp(i,1)  = rt1;
        results.accuracy(i,1) = accuracy;
        results.trl{i,1}      = roomname;
        results.trl{i,2}      = categories(1,i);
        
        % calculate cumulative SOT
        results.SOT.raw.fix  = SOT_f1;
        results.SOT.raw.stim = SOT_learn;
        results.SOT.raw.resp = SOT_r;
        results.SOT.raw.cue  = SOT_cue;
        results.SOT.raw.control= SOT_control;

    end
    
    % intermediate save
    dat.encoding.results = results;
    save([paths.behav num2str(ID) '_enc' num2str(BlockNum) '_backup.mat'],'dat')
    
    fprintf('\n============================================\n')    
end

%% wrap up
% terminate protocol
cd(paths.supports)
if BlockNum == 1
cgloadbmp(1,'instruction_ende1.bmp');
elseif BlockNum == 2
cgloadbmp(1,'instruction_ende_enc.bmp');
end
cgsetsprite(0);
cgdrawsprite(1,0,0); 
t_TaskEnd = cgflip(0,0,0);   t_TaskEnd_raw = toc; 
wait(3000);
cgshut; stop_cogent;

% calculate final cumulative SOT
results.SOT.raw.fix  = SOT_f1;
results.SOT.raw.stim = SOT_learn;
results.SOT.raw.resp = SOT_r;
results.SOT.raw.cue  = SOT_cue;
results.SOT.raw.control= SOT_control;
results.SOT.raw.end  = t_TaskEnd_raw;


% save results
dat.encoding.results = results;
dat.encoding.labels  = {'categories-> 1=cue1, 2=cue2, 0=control'; ...
    'results.trl: accuracy-> 1=correct, -1=internal, 0=external, NaN=missed)'};
dat.encoding.config  = config;
dat.encoding.config.keymap = map;
save([paths.behav behavfilename],'dat')

disp('****************************************')
disp('******* end of the encoding task *******')


end