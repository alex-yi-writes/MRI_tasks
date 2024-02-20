function  [data] = screen_display (TrialTask, Time2Target, Sham, xlocation1, xlocation2, TargetDisplayTime, ITI, RightKey, LeftKey) 

if rand<=0.5
    TargetPosition=xlocation1
else
    TargetPosition=xlocation2
end
% 
% if TrialTask>0
%     
% end


clearpict(1);
drawpict(1);

cgpencol(0.5,0.5,0.5);
cgrect(0,0,1000,1000);

cgpencol(1,1,1);
cgpenwid(5);
cgellipse(TargetPosition,0,70,70);

clearkeys;
TimeTarget=cgflip;
logstring(['time Target is ' mat2str(TimeTarget)]);
wait(TargetDisplayTime);

readkeys;
logkeys;
[KeyResp, KeyTime]=lastkeydown;

%this measure valid responses
if (KeyResp==LeftKey | KeyResp==RightKey )  
    RT=KeyTime-TimeTarget*1000;
    Key=KeyResp;
    if RT<2000
        if TargetPosition==xlocation1 & KeyResp==RightKey
            Response=1; %correct response
        elseif TargetPosition==xlocation2 & KeyResp==LeftKey
            Response=1;
        else
            Response=2; %incorrect responses
        end
    elseif RT>700 & ((TargetPosition==xlocation2 & KeyResp==RightKey) | (TargetPosition==xlocation1 & KeyResp==LeftKey))
        Response=3; 
    else
        Response=4;
    end
else
    RT=0;
    Key=0;
    KeyTime=0;
    Response=0;  %missed response
end

clearpict(3);
preparestring('+',3);
drawpict(3);
wait(1000);

clearpict(4);

if Response==1
    preparestring(['OK'],4);
elseif Response==2
    preparestring('falsche Taste',4);
elseif Response==3
    preparestring('falsche Taste und zu langsam',4);
elseif Response==4
    preparestring('zu langsam',4);
else
    preparestring('Sie haben nicht reagiert',4)
end
TimeOutcome=drawpict(4);
wait(1000);

clearpict(3);
preparestring('+',3);
drawpict(3);
wait(ITI);

data=[TimeTarget,TargetPosition,KeyResp,KeyTime,RT,Response,TargetDisplayTime,ITI];