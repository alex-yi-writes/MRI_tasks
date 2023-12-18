clc
clear

sync_box = serialport("COM4", 57600);

write(sync_box,"S","string")
pause(0.1)
confirmation = read(sync_box,1, "string")
pause(0.1)

for recieved = 1:1000
    trigger = read(sync_box, 1, "string")
end

% write(sync_box,"A","string")
% pause(0.1)
% confirmation = read(sync_box,1, "string")
% pause(0.1)
% write(sync_box,"D","string")
% pause(0.1)
% confirmation = read(sync_box,1, "string")
% pause(0.1)
% 
% clear sync_box;