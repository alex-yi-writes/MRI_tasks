function [rsp] = getKeys(keysWanted,waitTime) % keysWanted is a matrix list of keys you are 
                                           % waiting for e.g, [124 125 kbName('space')]
%   FlushEvents('keydown');
%   success = 0;
%   while success == 0
%    pressed = 0;
%    while pressed == 0
%     [pressed, secs, kbData] = KbCheck;
%    end;
%     for i = 1:length(keysWanted)
%       if kbData(keysWanted(i)) == 1
%        success = 1;
%        keyPressed = keysWanted(i);
%        FlushEvents('keydown');
%        break;
%       end;
%     end;
%     FlushEvents('keydown');
%    end;

% KbName('UnifyKeyNames');

% specify key names of interest in the study
activeKeys = keysWanted;

% set value for maximum time to wait for response (in seconds)
t2wait = waitTime; 
% if the wait for presses is in a loop, 
% then the following two commands should come before the loop starts
% restrict the keys for keyboard input to the keys we want
RestrictKeysForKbCheck(activeKeys);
% suppress echo to the command line for keypresses
ListenChar(2);
% get the time stamp at the start of waiting for key input 
% so we can evaluate timeout and reaction time
% tStart can also be the timestamp for the onset of the stimuli, 
% for example the VBLTimestamp returned by the 'Flip'
tStart = GetSecs;
% repeat until a valid key is pressed or we time out
timedout = false;
% initialise fields for rsp variable 
% that would contain details about the response if given
rsp.RT = NaN; rsp.keyCode = []; rsp.keyName = [];
while ~timedout,
    % check if a key is pressed
    % only keys specified in activeKeys are considered valid
    [ keyIsDown, keyTime, keyCode ] = KbCheck; 
      if(keyIsDown), break; end
      if( (keyTime - tStart) > t2wait), timedout = true; 
      rsp=0;
      end
  end
  % store code for key pressed and reaction time
  if(~timedout)
      rsp.RT            = keyTime - tStart;
      rsp.keyCode       = keyCode;
      rsp.keyName       = KbName(rsp.keyCode);
      rsp.rawKeytime    = keyTime;
      rsp.tStart        = tStart;
  end
% if the wait for presses is in a loop, 
% then the following two commands should come after the loop finishes
% reset the keyboard input checking for all keys
RestrictKeysForKbCheck([]);
% re-enable echo to the command line for key presses
% if code crashes before reaching this point 
% CTRL-C will reenable keyboard input
ListenChar(1)