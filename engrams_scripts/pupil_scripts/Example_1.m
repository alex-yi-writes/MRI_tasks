%% Example_1

% Clear the workspace and the display
sca;    
close all;       
clearvars; 

%loading default values into the Psychtoolbox
PsychDefaultSetup(2); 

%getting number of Screens to show the presentationwindow on the
%second monitor   
screens = Screen('Screens');
screenNumber = max(screens);
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber); 

%start a connection to the Eye-Tracking Software at the local pc on the port 5257
eye_connect('localhost',5257);

%get the name,the version and timestamp of the Eye-Tracking Software 
eye_get_version()
eye_get_timestamp()

%opening a psychwindow of 1200x800
[window, windowRect] = PsychImaging('OpenWindow',screenNumber,white, [30 50 830 650]);
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
ifi = Screen('GetFlipInterval', window);
%setting the text format
Screen('TextFont', window, 'Ariel');
Screen('TextSize', window, 30);
Screen('BlendFunction',window,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
[xCenter, yCenter] = RectCenter(windowRect);

%sending the setup to the Eye-Tracking Software 
eye_set_display_parameter(screenXpixels,screenYpixels,360,0.333);
eye_set_display_offset(30,50);

%starting the the calibration(if eye was detected once) and drawing each calibration 
%point to the screen until the Eye-Tracking Software is calibrated
if(eye_start_calibrate(9) == 0)
    %i is a counter checking the calibration
    %if the eye isn't found the counter increase; if it hits 10 the 
    %calibration stops
    i = 0;     
    while((eye_get_status() == 1) && (~KbCheck))
        if(eye_get_parameter('eye_status_detection') == 1)  
            dotColor = [255 255 255];
            dotSizePix = 20;
            i = i+1;
            if(i > 10)
                eye_stop_calibration();
                text = 'calibration stopped';
                %To show the text upside down the window must be mirrored
                Screen('glTranslate', window, xCenter, yCenter  , 0);
                Screen('glScale', window, 1, -1, 1);
                Screen('glTranslate', window, -xCenter, -yCenter, 0);
                DrawFormattedText(window, text, 'center', 'center', black);
                break
            end
        else
            dotColor = rand(1,3);
            dotSizePix = rand(1)*10  + 10;
        end
            %showing the calibration_point on the screen
            point = eye_get_calibration_point();
            Screen('DrawDots', window,[point(2) point(3)], dotSizePix, dotColor, [], 2);
            Screen('Flip', window);
    end    
end

 dotColor = [0 0 0];
 dotSizePix = 10;

% Loop the animation until a key is pressed
while((eye_get_status() == 2) && (~KbCheck))

    % Get the current position of the eye
    gaze = eye_get_gaze();
    %Drawing the current position of the eye (if gaze is available)
    if(gaze(1) == 0) 
        Screen('DrawDots', window,[gaze(2) gaze(3)] , dotSizePix, dotColor, [], 1);
    end
    
    % Flip to the screen
    Screen('Flip', window);

end

% Clear the screen and disconnect 
sca;
eye_disconnect();
 