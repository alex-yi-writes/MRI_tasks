function trigger = sync_trigger_bmmr(x)

sync_box = serialport("COM4", 57600);

write(sync_box,"S","string")
pause(0.1)
confirmation = read(sync_box,1, "string")
pause(0.1)

trigger_1= cell(1,x);

%     for filename = 1:x
%         trigger_anzahl{filename} = filename;
%     end


for recieved = 1:x
trigger = read(sync_box, 1, "string")
trigger_1 = trigger
%disp(recieved)
disp(trigger)
end

write(sync_box,"A","string")
pause(0.1)
confirmation = read(sync_box,1, "string")
pause(0.1)

write(sync_box,"D","string")
pause(0.1)
confirmation = read(sync_box,1, "string")
pause(0.1)

clear sync_box;
end