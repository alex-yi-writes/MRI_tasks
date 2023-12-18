% Function to submit the data and show confirmation
function submitData(idEdit, blockEdit, conditionDropdown, phaseDropdown, mriDropdown, pupilDropdown, fig)
% Get the values from input fields
id = idEdit.Value;
block = blockEdit.Value;
condition = conditionDropdown.Value;
phase = phaseDropdown.Value;
mri = mriDropdown.Value;
pupil = pupilDropdown.Value;

% Create a confirmation dialog
confirmationDlg = uiconfirm(fig, ...
    sprintf('ID: %f\nBlock Number: %f\nCondition: %s\nPhase: %s\nMRI: %s\nPupil Recording: %s', ...
    id, block, condition, phase, mri, pupil), ...
    'Confirm Inputs', ...
    'Options', {'Yes', 'No'}, ...
    'DefaultOption', 'No');

% Define a callback function for the confirmation dialog
confirmationDlg.ButtonPushedFcn = @(dlg, event) handleConfirmation(dlg, fig);

% Function to handle the confirmation dialog response
    function handleConfirmation(dlg, fig)
        if strcmp(dlg.SelectedObject.Text, 'Yes')
            % Close the dialog
            delete(fig);
        end
    end

end

