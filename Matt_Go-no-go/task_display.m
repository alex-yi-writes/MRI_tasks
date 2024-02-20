function  [data] = task_display (TrialCue, Time2Target, Sham, xlocation1, xlocation2, TargetDisplayTime, ITI, RandCSs, scr)

n=0;
Money=1;

if rand<=0.5
    TargetPosition=xlocation1
else
    TargetPosition=xlocation2
end

Pmat = [ones(1,7),2,2,2];
Pmat = Pmat(randperm(10));
P = Pmat(1);

clearkeys;
clearpict(2);
preparepict(loadpict(RandCSs{TrialCue}),2);
if scr==1
    outportb(888,255);
    wait(10);
    outportb(888,0);
end
TimeCue=drawpict(2);
wait(1000);

readkeys;
logkeys;
[KeyRespC1, KeyTimeC1, n]=getkeydown;

if n~=0
    ResponseCue1=1
else 
    ResponseCue1=0
end
   
clearkeys;
clearpict(3);
preparestring('+',3);
drawpict(3)
wait(Time2Target);

readkeys;
logkeys;
[KeyRespC2, KeyTimeC2, n]=getkeydown;

if n~=0
    ResponseCue2=1
else 
    ResponseCue2=0
end

if Sham==1 
    clearpict(1);
    drawpict(1);

    cgpencol(0.5,0.5,0.5);
    cgrect(0,0,1000,1000);

    cgpencol(1,1,1);
    cgpenwid(5);
    cgellipse(TargetPosition,0,70,70);

    clearkeys;
    TimeTarget=cgflip;
    logstring(['time Target is ' mat2str(TimeTarget*1000)]);
    wait(TargetDisplayTime);

    readkeys;
    logkeys;
    [KeyResp, KeyTime]=lastkeydown;

    %this measure valid responses
    if (KeyResp==80 | KeyResp==81 )
        RT=KeyTime-TimeTarget*1000;
        Key=KeyResp;
        if RT<700
            if TargetPosition==xlocation1 & KeyResp==80
                Response=1; %correct response
            elseif TargetPosition==xlocation2 & KeyResp==81
                Response=1;
            else
                Response=2; %incorrect response
            end
        else
            Response=3;
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
    
    x=[0,3,-3];
    y=[5,0,0];
    
    if TrialCue==1
        if Response==1 & P==1
            cgscale(60);
            cgpencol(0,1,0)
            cgrect(0,-4,2,8);
            cgpolygon(x,y);
            Won=Money,
        else
            cgscale(60);
            cgpencol(1,1,0)
            cgrect(0,0,8,2);
            Won=0;
        end
    elseif TrialCue==2
        if Response==1 & P==1
            cgscale(60);
            cgpencol(1,1,0);
            cgrect(0,0,8,2);
            Won=0;
        else
            cgscale(60);
            cgpencol(1,0,0);
            cgrect(0,4,2,8);
            cgpolygon(-x,-y);
            Won=-Money;
        end
    elseif TrialCue==3
        if Response==0 & P==1
            cgscale(60);
            cgpencol(0,1,0);
            cgrect(0,-4,2,8);
            cgpolygon(x,y);
            Won=Money;
        else
            cgscale(60);
            cgpencol(1,1,0);
            cgrect(0,0,8,2);
            Won=0;
        end
    elseif TrialCue==4
        if Response==0 & P==1
            cgscale(60);
            cgpencol(1,1,0);
            cgrect(0,0,8,2);
            Won=0;
        else
            cgscale(60);
            cgpencol(1,0,0)
            cgrect(0,4,2,8);
            cgpolygon(-x,-y);
            Won=-Money;
        end
    end
    TimeOutcome=drawpict(4);
    wait(1000);
    cgscale;
    
    clearpict(3);
    preparestring('+',3);
    drawpict(3);
    wait(ITI);
else
    clearpict(3);
    preparestring('+',3);
    drawpict(3);
    wait(500);
    TimeTarget=0;
    KeyResp=0;
    KeyTime=0;
    RT=0;
    Response=0;
    TargetDisplayTime=0;
    Won=0;
    TimeOutcome=0;
end

data=[TrialCue,TimeCue,ResponseCue1,Time2Target,TimeTarget,ResponseCue2,Sham,TargetPosition,KeyResp,KeyTime,RT,Response,TargetDisplayTime,TimeOutcome,ITI,Won];
