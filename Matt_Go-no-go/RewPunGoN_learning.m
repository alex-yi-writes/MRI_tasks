%%%%%%%%%%%%%%%
%%  Variables  %
%%%%%%%%%%%%%%%%

cd('\\fs-md\users\yiy\Desktop\Matt_Go-no-go\')


%%%% Before beginning experiment must switch on line 25 of this script and
%%%% line 18-20 of the 2 task display (learning and task)

%%subject-dependent
subj_name       =input('Subject ID: ', 's');
subj_number     =input('Subject number: ');
if length(subj_number)==0
    error ('Subject number missing ')
end
randomization=input('press 0 if stimuli have not been randomized, press 1 if they have been randomized ');
scr = input('press 0 if no scr, 1 if scr is measured')
scanner = input('press 0 if no scanner and 1 if scanner')
% session = input ('enter 0 for training and 1 to 4 for task ')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

config_display(1, 3, [0.5 0.5 0.5], [1 1 1], 'Helvetica', 30, 7, 0);
config_keyboard;
% config_data('Pictures2.dat');
config_log( ['Cogent' num2str(datevec(now),'-%02.0f') subj_name '.log'] );
if scr==1
    startportb(888);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% buffers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
InstructScreen=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparing matrices and stimuli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch randomization
    case 0
        Stimuli={'S1.bmp','S2.bmp','S3.bmp','S4.bmp'};
        StimOrd=randperm(4)
        save(['StimOrd_' num2str(subj_number) '_' subj_name],'StimOrd')
        CSs=[Stimuli(StimOrd(1)),Stimuli(StimOrd(2)),Stimuli(StimOrd(3)),Stimuli(StimOrd(4))];
        RandCSs=CSs
        save(['RandCSs_' num2str(subj_number) '_' subj_name],'RandCSs')
    case 1
        load (['RandCSs_' num2str(subj_number) '_' subj_name])
end

%%% Trial type matrix
TrialTypeLearning=repmat([1 2 3 4],[1,60]);  % 1 go reward; 2 go punishment; 3 no-go reward; 4 no-go punishment
RandTrialTypeLearning=TrialTypeLearning(randperm(size(TrialTypeLearning,2)));

% TrialType=repmat([1 2 3 4],[1,20]);
% RandTrialType=TrialType(randperm(size(TrialType,2)));

TrialTarget=repmat([1 2],[1,60]);   % 1 real trial; 2 sham trial
RandTrialTarget=TrialTarget(randperm(size(TrialTarget,2)));

TrialNumber=[randperm(length(TrialTarget))' randperm(length(TrialTarget))' randperm(length(TrialTarget))' randperm(length(TrialTarget))'];

TrialCounter=[0 0 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TargetTime=250 + rand(1,length(TrialTypeLearning))*3250;
TargetDisplayTime=1500;
ITI=750 + rand(1,length(TrialTypeLearning))*500;

if scanner==1
    port=1;  %for scanner
    config_serial(1);
    config_keyboard_monitoring('led');
    %%keys
    LeftKey         =80;              %m
    RightKey        =81;
else
    LeftKey         =97;              %m
    RightKey        =98;
end

xlocation1     =-100;
xlocation2     =+100;
ylocation      =0;

NumDummies=6;
SlicesVol=35;

start_cogent;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Target practice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

preparestring('Probedurchlauf', InstructScreen,0,275);
preparestring('Sobald Sie den Kreis sehen',1,0,110);
preparestring('geben Sie mit Knopfdruck an auf welcher Seite',1,0,70);
preparestring('wenn Sie es rechtzeitig korrekt ausf�hren erhalten Sie eine OK Nachricht',1,0,30);
preparestring('wenn Sie versagen wird Ihnen die Ursache mitgeteilt',1,0,-10);
preparestring('...Dr�cken Sie Leertaste zum Fortfahren', InstructScreen,250,-275);
drawpict(InstructScreen);
waitkeydown(inf,71);
%     wait(3000);
clearpict(InstructScreen);

for count=1:10
    ITI_trial=ITI(count);
    [data] = screen_display (0,0,1,xlocation1,xlocation2,TargetDisplayTime,ITI_trial,LeftKey,RightKey);
    PracticeData(count,:)=data;
end
save([num2str(subj_number) '_' subj_name '_PracticeData'],'PracticeData');
% eval( ['save Cogent_results_' subj_name '_' num2str(subj_number) '_learning_' num2str(datevec(now),'-%02.0f')] );


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Cues and Task Learning
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
preparestring('Aufgabenlernen', InstructScreen,0,275);
preparestring('Zun�chst werden Sie das erste Bild sehen',1,0,150);
preparestring('Auf die Bilder wird der Kreis folgen:',1,0,70);
preparestring('manchmal m�ssen Sie die Position des Kreises angeben,',1,0,30);
preparestring('manchmal m�ssen Sie nichts tun.',1,0,-10);
preparestring('In jedem Durchgang werden Sie sehen ob Sie 50 Cent gewonnen haben oder nicht.',1,0,-50);
preparestring('Ihre Aufgabe wird es sein die beste Antwort nach jedem Bild herauszufinden',1,0,-90);
preparestring('um Ihren Gesamtgewinn zu maximieren.',1,0,-130);
preparestring('...Dr�cken Sie Leertaste zum Fortfahren',InstructScreen,250,-275);
drawpict(InstructScreen);
waitkeydown(inf,71);
clearpict(InstructScreen);

preparestring('Bereiten Sie sich vor',InstructScreen,0,0);
drawpict(InstructScreen);
wait(1500);
clearpict(InstructScreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Learning the task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TargetTime=250 + rand(1,length(TrialTypeLearning))*3250 %in primary behavioural experiment 1750;
TargetDisplayTime=1500;
ITI=750 + rand(1,length(TrialTypeLearning))*750;

for count=1:length(RandTrialTypeLearning)
    TrialCue=RandTrialTypeLearning(count);
    Time2Target=TargetTime(count);
    ITI_trial=ITI(count);
    if rem(count,60)==0
        preparestring('Sie d�rfen sich nun ausruhen',InstructScreen,0,0);
        preparestring('...Dr�cken Sie Leertaste zum Fortfahren', InstructScreen,250,-275);
        drawpict(InstructScreen);
        waitkeydown(inf,71);
        clearpict(InstructScreen);
    end
    [data] = task_display_learning (TrialCue,Time2Target,1,xlocation1,xlocation2,TargetDisplayTime,ITI_trial,RandCSs,scr,LeftKey,RightKey);
    TaskDataLearning(count,:)=[count data];
end
save([num2str(subj_number) '_' subj_name '_TaskDataLearning'],'TaskDataLearning');
% eval( ['save Cogent_results_' subj_name '_' num2str(subj_number) '_learning_' num2str(datevec(now),'-%02.0f')] );

preparestring(['Sie haben in dieser Runde  ' num2str(sum(TaskDataLearning(:,17))) 'euro gewonnen'],InstructScreen)
drawpict(InstructScreen);
waitkeydown(inf,71);
clearpict(InstructScreen);

% preparestring('+',InstructScreen)
% drawpict(InstructScreen);
% waitkeydown(inf,71);
% clearpict(InstructScreen);

stop_cogent;
