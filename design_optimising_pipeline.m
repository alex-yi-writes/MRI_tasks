%% create and reconfigure design matrix

close all;
clear;
warning off

addpath(genpath('/Users/yeojin/Desktop/B_scripts/BD_tasks/task_modelling')) % path to your optimisation functions
addpath(genpath('/Users/yeojin/Documents/MATLAB/CanlabCore-master')) % path to the GA toolbox

%% main task

%% CB1&2 DAY1: private(3)-90, public(4)-90


for ii = 1:2 % number of unique matrices you'd like to create
    
    cnt = 0; errs = 0; % reset the while loop conditional
    while cnt == errs % try until it doesn't error
        try
            close all;
            
            cd('/Users/yeojin/Desktop/B_scripts/BD_tasks/task_scripts/MRPET/models/transmat/') % where would you like to save your models and figures?
            
            [meanrecipvif, vifs, design_struct] = generate_jittered_er_design_MainTask; % run the specific GA optimisation function
            
            % save diagnostics figure
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(['GA_MainTask_ver' num2str(ii) '_Diagnostics_D1'],'-dpng','-r0')
            close figure 1
            
            % save predicted HRF function figure
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(['GA_MainTask_ver' num2str(ii) '_PredictedActivity_D1'],'-dpng','-r0')
            close figure 2
            
            % tailor the design matrix to your task
            tmp = design_struct.eventlist(:,3);
            tmp2(find(tmp==1),1)= 1;
            tmp2(find(tmp==2),1)= 2;
%             tmp2(find(tmp==3),1)= 1;
%             tmp2(find(tmp==4),1)= 2; % 3:1 reward ratio
            design_struct.eventlist(:,3) = tmp2; % assign the new regime
            
            % save the design matrix and efficiency infos
            save(['/Users/yeojin/Desktop/B_scripts/BD_tasks/task_scripts/MRPET/models/transmat/Main_GA_designstruct_D1_ver' num2str(ii) '.mat' ],...
                'design_struct','meanrecipvif','vifs')
            keep ii cnt errs
            
        catch WhatWasTheError
            errs = errs+1; 
        end % close the try loop
        cnt = cnt+1;
    end % close the while loop
    
end % close the for loop


%% CB1&2 DAY2: nature(1)-45, urban(2)-135

close all
clear

for ii2 = 1:2
    
    cnt = 0; errs = 0;
    while cnt == errs
        try
            close all;
            cd('/Users/yeojin/Desktop/B_scripts/BD_tasks/task_scripts/MRPET/models/transmat/')
            
            [meanrecipvif, vifs, design_struct] = generate_jittered_er_design_MainTask;
            
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(['GA_MainTask_ver' num2str(ii2) '_Diagnostics_D2'],'-dpng','-r0')
            close figure 1
            
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(['GA_MainTask_ver' num2str(ii2) '_PredictedActivity_D2'],'-dpng','-r0')
            close figure 2
            
            tmp = design_struct.eventlist(:,3);
            tmp2(find(tmp==1),1)= 1;
            tmp2(find(tmp==2),1)= 2;
%             tmp2(find(tmp==3),1)= 2;
%             tmp2(find(tmp==4),1)= 2;
            design_struct.eventlist(:,3) = tmp2;
            
            save(['/Users/yeojin/Desktop/B_scripts/BD_tasks/task_scripts/MRPET/models/transmat/Main_GA_designstruct_D2_ver' num2str(ii2) '.mat' ],'design_struct','meanrecipvif','vifs')
            keep ii cnt errs
            
        catch WhatWasTheError
            errs = errs+1;
        end
        cnt = cnt+1;
    end
    
end

%% memory test

%% CB1&2 DAY1: nature(2)-135, urban(1)-45

close all
clear

for ii = 1:2
    cnt = 0; errs = 0;
    while cnt == errs
        try
            close all;
            
            cd('/Users/yeojin/Desktop/B_scripts/BD_tasks/task_scripts/MRPET/models/transmat/')
            
            [meanrecipvif, vifs, design_struct] = generate_jittered_er_design_MemoryTest;
            
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(['GA_ImmMem_ver' num2str(ii) '_Diagnostics_D1'],'-dpng','-r0')
            close figure 1
            
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(['GA_ImmMem_ver' num2str(ii) '_PredictedActivity_D1'],'-dpng','-r0')
            close figure 2
            
            tmp = design_struct.eventlist(:,3);
            tmp2(find(tmp==1),1)= 11;
            tmp2(find(tmp==2),1)= 21;
            tmp2(find(tmp==3),1)= 12;
            tmp2(find(tmp==4),1)= 22;
%             tmp2(find(tmp==5),1)= 12;
%             tmp2(find(tmp==6),1)= 12;
%             tmp2(find(tmp==7),1)= 12;
%             tmp2(find(tmp==8),1)= 22;
            design_struct.eventlist(:,3) = tmp2;
            
            save(['/Users/yeojin/Desktop/B_scripts/BD_tasks/task_scripts/MRPET/models/transmat/ImmMem_GA_designstruct_D1_ver' num2str(ii) '.mat' ],'design_struct','meanrecipvif','vifs')
            keep ii cnt errs
            
        catch WhatWasTheError
            errs = errs+1;
        end
        cnt = cnt+1;
    end
    
end


%% CB1&2 DAY2: private(3)-45, public(4)-135

for ii2 = 1:2
    cnt = 0; errs = 0;
    while cnt == errs
        try
            close all;
            cd('/Users/yeojin/Desktop/B_scripts/BD_tasks/task_scripts/MRPET/models/transmat/')
            
            [meanrecipvif, vifs, design_struct] = generate_jittered_er_design_MemoryTest;
            
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(['GA_ImmMem_ver' num2str(ii2) '_Diagnostics_D2'],'-dpng','-r0')
            close figure 1
            
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(['GA_ImmMem_ver' num2str(ii2) '_PredictedActivity_D2'],'-dpng','-r0')
            close figure 2
            
            tmp = design_struct.eventlist(:,3);
            tmp2(find(tmp==1),1)= 11;
            tmp2(find(tmp==2),1)= 21;
            tmp2(find(tmp==3),1)= 12;
            tmp2(find(tmp==4),1)= 22;
%             tmp2(find(tmp==5),1)= 12;
%             tmp2(find(tmp==6),1)= 22;
%             tmp2(find(tmp==7),1)= 22;
%             tmp2(find(tmp==8),1)= 22;
            design_struct.eventlist(:,3) = tmp2;
            
            save(['/Users/yeojin/Desktop/B_scripts/BD_tasks/task_scripts/MRPET/models/transmat/ImmMem_GA_designstruct_D2_ver' num2str(ii2) '.mat' ],'design_struct','meanrecipvif','vifs')
            keep ii cnt errs
            
        catch WhatWasTheError
            errs = errs+1;
        end
        cnt = cnt+1;
    end
end

