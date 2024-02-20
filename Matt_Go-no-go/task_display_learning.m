function  [data] = task_display (TrialCue, Time2Target, Sham, xlocation1, xlocation2, TargetDisplayTime, ITI, RandCSs, scr, KeyLeft, KeyRight)

Money=0.5;

if rand<=0.5
    TargetPosition=xlocation1;
else
    TargetPosition=xlocation2;
end

Pmat = [ones(1,8),2,2];
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
[KeyResp, KeyTime, n]=getkeydown;

if n~=0
    clearpict(3);
    cgscale(60);
    cgfont('Arial',2);
    preparestring('Sie haben die Taste zu früh gedrückt',3);
    drawpict(3);
    wait(1000);
    cgscale;
    
    TimeTarget=0; KeyResp=0; KeyTime=0; RT=0; ResponseCue1=1; ResponseCue2=0; Response=0; TargetDisplayTime=0; TimeOutcome=0; Won=0;
else

    ResponseCue1=0
    
    clearkeys;
    clearpict(3);
    preparestring('+',3);
    drawpict(3);
    wait(Time2Target);

    readkeys;
    logkeys;
    [KeyResp, KeyTime, n]=getkeydown;

    if n~=0
        clearpict(3);
        cgscale(60);
        cgfont('Arial',2);
        preparestring('Sie haben die Taste zu früh gedrückt',3);
        drawpict(3);
        wait(1000);
        cgscale;
        TimeTarget=0; KeyResp=0; KeyTime=0; RT=0; ResponseCue2=1; Response=0; TargetDisplayTime=0; TimeOutcome=0; Won=0;
    else
        ResponseCue2=0;
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
        if (KeyResp==KeyLeft | KeyResp==KeyRight )
            RT=KeyTime-TimeTarget*1000;
            Key=KeyResp;
            if RT<2000
                if TargetPosition==xlocation1 & KeyResp==KeyLeft
                    Response=1; %richtige Antwort
                elseif TargetPosition==xlocation2 & KeyResp==KeyRight
                    Response=1;
                else
                    Response=2; %falsche Antwort
                end
            else
                Response=3; % zu spät
            end
        else
            RT=0;
            Key=0;
            KeyTime=0;
            Response=0;  %keine Antwort
        end

        clearpict(3);
        preparestring('+',3);
        drawpict(3);
        wait(1000);

        clearpict(4);

        x=[0,3,-3];
        y=[5,0,0];

        if TrialCue==1
            if (Response==1 & P==1) | (Response==0 & P==2)
                Won=Money,
                cgscale(60);
                cgpencol(0,1,0)
                cgrect(0,-4,2,8);
                cgpolygon(x,y);
%                 cgpencol(1,1,1);
%                 cgfont('Arial',2)
%                 cgtext(['OK, you won £' num2str(Money)],0,-10);
   
            else
                Won=0;
                cgscale(60);
                cgpencol(1,1,0)
                cgrect(0,0,8,2);
%                 cgpencol(1,1,1);
%                 cgfont('Arial',2);
%                 if Response==1
%                     cgtext(['OK, but you did not win£' num2str(Money) 'anyway'],0,-10);
%                 else
%                     cgtext('UNSUCCESSFUL, your response was too late or incorrect',0,-10);
%                 end
            end
        elseif TrialCue==2
            if (Response==1 & P==1) | (Response==0 & P==2)
                Won=0;
                cgscale(60);
                cgpencol(1,1,0);
                cgrect(0,0,8,2);
%                 cgpencol(1,1,1);
%                 cgfont('Arial',2)
%                 cgtext(['OK, you avoided losing £' num2str(Money)],0,-10);
            else
                Won=-Money;
                cgscale(60);
                cgpencol(1,0,0);
                cgrect(0,4,2,8);
                cgpolygon(-x,-y);
%                 cgpencol(1,1,1);
%                 cgfont('Arial',2)
%                 if Response==1
%                     cgtext(['OK, but you lost £' num2str(Money) 'anyway'],0,-10);
%                 else
%                     cgtext('UNSUCCESSFUL, your response was too late or incorrect',0,-10);
%                 end
            end
        elseif TrialCue==3
            if (Response==0 & P==1)  | (Response>0 & P==2)
                Won=Money;
                cgscale(60);
                cgpencol(0,1,0);
                cgrect(0,-4,2,8);
                cgpolygon(x,y);
%                 cgpencol(1,1,1);
%                 cgfont('Arial',2)
%                 cgtext(['OK, you won £' num2str(Money)],0,-10);
            else
                Won=0;
                cgscale(60);
                cgpencol(1,1,0);
                cgrect(0,0,8,2);
%                 cgpencol(1,1,1);
%                 cgfont('Arial',2)
%                 if Response==0
%                     cgtext(['OK, but you did not win £' num2str(Money) 'anyway'],0,-10);
%                 else
%                     cgtext('UNSUCCESSFUL, you should have not responded',0,-10);
%                 end
            end
        elseif TrialCue==4
            if (Response==0 & P==1)  | (Response>0 & P==2)
                Won=0
                cgscale(60);
                cgpencol(1,1,0);
                cgrect(0,0,8,2);
%                 cgpencol(1,1,1);
%                 cgfont('Arial',2)
%                 cgtext(['OK, you avoided losing £' num2str(Money)],0,-10);
            else
                Won=-Money;
                cgscale(60);
                cgpencol(1,0,0)
                cgrect(0,4,2,8);
                cgpolygon(-x,-y);
%                 cgpencol(1,1,1);
%                 cgfont('Arial',2)
%                 if Response==0
%                     cgtext(['OK, but you lost £' num2str(Money) 'anyway'],0,-10);
%                 else
%                     cgtext('UNSUCCESSFUL, you should have not responded',0,-10);
%                 end
            end
        end
        TimeOutcome=drawpict(4);
        wait(2000);
        cgscale;
    end
end
clearpict(3);
preparestring('+',3);
drawpict(3);
wait(ITI);

% cgfont('Helvetica',3)

data=[TrialCue,TimeCue,ResponseCue1,ResponseCue2,Time2Target,TimeTarget,TargetPosition,KeyResp,KeyTime,RT,Response,P,TargetDisplayTime,TimeOutcome,ITI,Won];
