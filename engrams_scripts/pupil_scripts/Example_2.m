%% Example_2

%clear the display
sca;  
close all;        
clearvars; 
 
%loading some default values into the Psychtoolbox 
PsychDefaultSetup(2);
 
%getting number of Screens to show show the presentation window on the
%second monitor
screens = Screen('Screens');
screenNumber = max(screens);
white = WhiteIndex(screenNumber); 
black = BlackIndex(screenNumber);

%start a connection to the Eye-Tracking Software at the local pc on the port 5257
eye_connect('localhost',5257);  

%get the version of the Eye-Tracking Software 
eye_get_version()

%enable saving both tracking data and video
eye_set_parameter('eye_save_tracking_and_video','true');
%start saving the tracking and video data
eye_set_parameter('eye_save_tracking','true');
%opening a psychwindow (800*600)
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, white, [30 50 830 650]);
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
ifi = Screen('GetFlipInterval', window); 
%setting the text format
Screen('TextFont', window, 'Ariel');
Screen('TextSize', window, 30);
Screen('BlendFunction',window,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
[xCenter, yCenter] = RectCenter(windowRect);

%setting the Eye-Tracking display
eye_set_displaymode(screenXpixels,screenYpixels);
eye_set_display_offset(30,50);

text = 'To calibrate watch directly at the points';

%To show the text upside down the window must be mirrored
Screen('glTranslate', window, xCenter, yCenter, 0);
Screen('glScale', window, 1, -1, 1);
Screen('glTranslate', window, -xCenter, -yCenter, 0);
DrawFormattedText(window, text, 'center', 'center', black);
%show text on the screen
Screen('Flip', window);

%wait for the person to read the text
pause(3);
eye_set_software_event('Start calibration');
%start the the calibartion and draw each calibration point to the
%psychwindow until the Software is calibrated
if(eye_start_calibrate(9) == 0)
    while((eye_get_status() == 1) && (~KbCheck))
        point = eye_get_calibstate();
        dotColor = rand(1,3);
        dotSizePix = rand(1)*10  + 10;
        Screen('DrawDots', window,[point(1) point(2)], dotSizePix, dotColor, [], 2);
        Screen('Flip', window);
    end    
end


dotSizePix = 10; 
dotColor = [0 0 0]; 

eye_set_software_event('Start point catching');
text = 'Try to catch the point below';
DrawFormattedText(window, text, 'center', 'center', black);     
pos = [screenXpixels/2 screenYpixels*2/3];  
Screen('DrawDots', window,pos, dotSizePix ,dotColor, [], 2);
Screen('Flip', window); 
counter = 0;

%starts stream for faster events 
eye_start_stream(1);

%while the Software is connected random points appear until you look at it.
%if the gaze is in a defined area around the random point and stays there
%15 detected gaze points will  be count until the point disappears.
while((eye_get_status() == 2) && (~KbCheck))
  
   gaze = eye_get_gaze();
     
     if(gaze(1) == 0)
       if((gaze(2)>(pos(1)-30)) && (gaze(2)<(pos(1)+30)) && (gaze(3)>(pos(2)-30)) && (gaze(3)<(pos(2)+30)))
           counter = counter + 1;
       end 
    end

    if counter>15       
        pos = [screenXpixels*(0.1+ 0.8*rand(1)) screenYpixels*(0.1+ 0.8*rand(1))]; 
        
        Screen('DrawDots', window,[pos(1) pos(2)], dotSizePix, dotColor, [], 2);
     
        Screen('Flip',window);
        
        counter = 0;
    end    
end

%at the end of the example the recording and stream are stopped and the last error
%description caught
eye_set_software_event('Stop point catching');
%stop saving
eye_set_parameter('eye_save_tracking','false');

last_error = eye_get_last_error();
%disconnect from eyetracking system
eye_disconnect();

% Clear the screen
sca;