% STARTER SCRIPT : ENGRAMS
%  In emergency contact Alex:
%              +49 162 713 20 49
%              yeo-jin.yi.15@ucl.ac.uk

%%  work log

%   12_12_2023      created the script

%% experiment preparation
%  always run this cell before starting an experiment.434

clear all; close all; clc;

taskpath       = 'E:/ENGRAMS/';
savepath_behav = [taskpath 'data/behav/'];
logpath        = [taskpath 'data/logs/'];
addpath(genpath([taskpath 'scripts']))

% Create a UIFigure with a larger size
fig = uifigure('Name', 'Data Entry Dialog', 'Position', [100, 100, 300, 500]);

% Create labels and input components with improved spacing
xPosition = 50;
yPosition = 470; % Adjusted yPosition for better spacing
spacing = 70;

idLabel             = uilabel(fig, 'Text', 'ID:', 'Position', [xPosition, yPosition, 100, 22]);
idEdit              = uieditfield(fig, 'numeric', 'Position', [xPosition, yPosition - 30, 100, 22]);
yPosition           = yPosition - spacing;

blockLabel          = uilabel(fig, 'Text', 'Block Number:', 'Position', [xPosition, yPosition, 150, 22]);
blockEdit           = uieditfield(fig, 'numeric', 'Position', [xPosition, yPosition - 30, 100, 22]);
yPosition           = yPosition - spacing;

conditionLabel      = uilabel(fig, 'Text', 'Condition:', 'Position', [xPosition, yPosition, 100, 22]);
conditionDropdown   = uidropdown(fig, 'Items', {'Emotional', 'Neutral'}, 'Position', [xPosition, yPosition - 30, 100, 22]);
yPosition           = yPosition - spacing;

phaseLabel          = uilabel(fig, 'Text', 'Phase:', 'Position', [xPosition, yPosition, 100, 22]);
phaseDropdown       = uidropdown(fig, 'Items', {'Original', 'Recombination'}, 'Position', [xPosition, yPosition - 30, 100, 22]);
yPosition           = yPosition - spacing;

TaskTypeLabel       = uilabel(fig, 'Text', 'Task Type:', 'Position', [xPosition, yPosition, 100, 22]);
TaskTypeDropdown    = uidropdown(fig, 'Items', {'Encoding', 'Recognition'}, 'Position', [xPosition, yPosition - 30, 100, 22]);
yPosition           = yPosition - spacing;

mriLabel            = uilabel(fig, 'Text', 'MRI:', 'Position', [xPosition, yPosition, 150, 22]);
mriDropdown         = uidropdown(fig, 'Items', {'Yes', 'No'}, 'Position', [xPosition, yPosition - 30, 100, 22]);
yPosition           = yPosition - spacing;

pupilLabel          = uilabel(fig, 'Text', 'Pupil Recording:', 'Position', [xPosition, yPosition, 150, 22]);
pupilDropdown       = uidropdown(fig, 'Items', {'Yes', 'No'}, 'Position', [xPosition, yPosition - 30, 100, 22]);
yPosition           = yPosition - spacing;


% Create a button to submit the data
submitBtn = uibutton(fig, 'push', 'Text', 'Submit', 'Position', [160, 10, 80, 30]);

% Define a callback function for the submit button
submitBtn.ButtonPushedFcn = @(btn, event) submitData(idEdit, blockEdit, conditionDropdown, phaseDropdown, TaskTypeDropdown, mriDropdown, pupilDropdown, fig);

% Display the UIFigure
fig.Visible = 'on';
waitfor(fig);

%% Start Experiment

if strcmp(confirmation,'Yes')==1
    
    disp('booting up...')
    
    if strcmp(tasktype,'Encoding')==1
        
        if strcmp(phase,'Original')==1
            behavfilename  = [num2str(id) '_enc_orig_' num2str(block) '.mat'];
        elseif strcmp(phase,'Recombination')==1
            behavfilename  = [num2str(id) '_enc_recombi_' num2str(block) '.mat'];
        end
        
        [dat] = engrams_encoding(id,block,condition,phase,mri,pupil,taskpath)
        
    elseif strcmp(tasktype,'Recognition')==1
        if strcmp(phase,'Original')==1
            behavfilename  = [num2str(id) '_rcg_orig_' num2str(block) '.mat'];
        elseif strcmp(phase,'Recombination')==1
            behavfilename  = [num2str(id) '_rcg_recombi_' num2str(block) '.mat'];
        end
        
        [dat] = engrams_recognition(id,block,condition,phase,mri,pupil,taskpath)
        
    end
    
    % save results to make it doubly sure!
    save([savepath_behav behavfilename],'dat')
    
    disp('****************************************')
    
elseif strcmp(confirmation,'No')==1
    % Display the pop-up message
    choice = questdlg('Run the starter script again fresh!', 'Message', 'Proceed', 'Proceed');
    
else
    error('check inputs again')
end

%% functions

% Function to submit the data, display confirmation, and save as variables
function submitData(idEdit, blockEdit, conditionDropdown, phaseDropdown, TaskTypeDropdown, mriDropdown, pupilDropdown, fig)

% Get the values from input fields
id = idEdit.Value;
block = blockEdit.Value;
condition = conditionDropdown.Value;
phase = phaseDropdown.Value;
tasktype = TaskTypeDropdown.Value;
mri = mriDropdown.Value;
pupil = pupilDropdown.Value;

% Save the collected data as variables
assignin('base', 'id', id);
assignin('base', 'block', block);
assignin('base', 'condition', condition);
assignin('base', 'phase', phase);
assignin('base', 'tasktype', tasktype);
assignin('base', 'mri', mri);
assignin('base', 'pupil', pupil);

% Create a summary of the user's input
userInputSummary = sprintf('Do you confirm the following entered data? \n\nID: %4.0f\nBlock Number: %2.0f\nCondition: %s\nPhase: %s\nTask Type: %s\nMRI: %s\nPupil Recording: %s', ...
    id, block, condition, phase, tasktype, mri, pupil);

% Ask for confirmation and show the user's input
confirmationMessage = userInputSummary;
confirmation = uiconfirm(fig, confirmationMessage, 'Confirmation', 'Options', {'Yes', 'No'}, 'DefaultOption', 1, 'CancelOption', 2);

% Save the confirmation response
assignin('base', 'confirmation', confirmation);

% Close the figure after confirmation
close(fig);

if strcmp(confirmation,'Yes')==1
    disp('parameters confirmed')
elseif strcmp(confirmation,'No')==1
    % Display the pop-up message
    choice = questdlg('Run the script again!', 'Message', 'Proceed', 'Proceed');
    
else
end
end

% % Function to handle the confirmation dialog response
% function handleConfirmation(dlg)
% % Check the selected option
% if strcmp(dlg.SelectedObject.Text, 'Yes')
%     confirmations(end+1) = 1; % Store 'Yes' as 1
% else
%     confirmations(end+1) = 0; % Store 'No' as 0
% end
% end
