% Size ERP 
% Screen('Preference', 'SkipSyncTests', 1);
clear variables

% subject = 'TEST';

%%%
daq_id = DaqDeviceIndex;
if ~isempty(daq_id), DaqDConfigPort(daq_id,0,0); end

% set the trigger value
START = 1;
SIZE1 = 11;
SIZE2 = 12;
SIZE3 = 13;
SIZE4 = 14;
SIZE5 = 15;
TRIALEND = 99;
END = 32;d

%%% display params
dp.screenNum = max(Screen('Screens'));

dp.dist = 55;                  % (cm) 
dp.width = 60;                 % 
dp.bkColor = 0.5;
dp.textColor = zeros(1,3);

dp.screenNum = max(Screen('Screens'));
d = Screen('Resolution',dp.screenNum);
dp.resolution = [d.width d.height];
dp.frameRate = d.hz;
dp.aspectRatio = dp.resolution(2) / dp.resolution(1);
dp.ppd = dp.resolution(1)/((2*atan(dp.width/2/dp.dist))*180/pi);

% experiment params
s = rng('shuffle');
param.randSeed = s;

param.trialDur = 0.25;
param.tf = 8;               %Hz - counter-phase flicker
param.f2 = 2*param.tf;
param.sf = 5;
param.period = param.trialDur;
param.numFlipsPerCycle = param.period*param.tf*2;   %half cycle (default)
param.numTextures = param.numFlipsPerCycle;     %default

param.fixSize = 0.3*dp.ppd;         %fixation point size
param.ovalSize = 0.7*dp.ppd; 
param.height = dp.resolution(2)/dp.ppd; % in degree
param.iti = 3*ones(1,6) + [0, 0.1, 0.2, 0.3, 0.4, 0.5];

param.useKbQueueCheck = 1;

cond.size = round(10.^linspace(log10(3.5), log10(12),5),1); %[3, 6, 9, 12, 15];
cond.nReps = 70;
cond.nTrials = length(cond.size) * cond.nReps;
cond.trial = [];
for i = 1: cond.nReps
    cond.trial = [cond.trial randperm(length(cond.size))];
end

%%% create textures
%make cartesian matries
[x,y] = meshgrid( linspace(-1/dp.aspectRatio,1/dp.aspectRatio, dp.resolution(1)), linspace(-1,1,dp.resolution(2)));

%make polar angle matrices
r = sqrt(x.^2 + y.^2);

% str = sprintf('making textures');
% h = waitbar(0,str);
for i = 1: length(cond.size)
    check = sign(sin(2*pi*param.sf*x) .* sin(2*pi*param.sf*y));
    for j = 1:param.numTextures
        boundH = cond.size(i)/param.height;
        check(r>boundH) = 0;
        img = check; %reset your image to the basic polar checkerboard each time
        %flip it from B/W to W/B
        img = (-1)^j*img;

        %make and draw your texture
        texture(:,:,j,i) = (img+1)*0.5;

        %advance waitbar
%         waitbar(j/param.numTextures, h);
    end
end
% close(h)

%%%
ListenChar(0);
HideCursor;
dp = OpenWindow(dp);
Screen('Flip', dp.wPtr);
kb = init_keyboard(param);    

%for fixation square
fixationOut = [(dp.resolution(1)-param.fixSize)/2, (dp.resolution(2)-param.fixSize)/2, (dp.resolution(1)+param.fixSize)/2, (dp.resolution(2)+param.fixSize)/2];
fixationIn = [(dp.resolution(1)-param.fixSize+5)/2, (dp.resolution(2)-param.fixSize+5)/2, (dp.resolution(1)+param.fixSize-5)/2, (dp.resolution(2)+param.fixSize-5)/2];
fixationOval = [(dp.resolution(1)-param.ovalSize)/2, (dp.resolution(2)-param.ovalSize)/2, (dp.resolution(1)+param.ovalSize)/2, (dp.resolution(2)+param.ovalSize)/2];
%% start
if ~isempty(daq_id), err=DaqDOut(daq_id,0,START); end
str = 'Press space key to start';
strBounds = Screen('TextBounds', dp.wPtr, str);
Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-strBounds(4)/2-50, [0 0 0]);
Screen('FillRect', dp.wPtr, [1 1 1], fixationOut);
Screen('FillRect', dp.wPtr, [0 0 0], fixationIn);
Screen('Flip', dp.wPtr);
isResponse = 0;
while ~isResponse
    [keyIsDown, secs, keyCode] = KbCheck(-1);
    if keyCode(kb.spaceKey)
        isResponse = 1;
    elseif keyCode(kb.escKey)
        sca;
    end
end

Screen('FillRect', dp.wPtr, [1 1 1], fixationOut);
Screen('FillRect', dp.wPtr, [0 0 0], fixationIn);
Screen('Flip', dp.wPtr);
startTime = GetSecs;
randITI = randperm(length(param.iti));
while (GetSecs - startTime) < param.iti(randITI(1)), ;, end

%% animate textures
for iTrial = 1: cond.nTrials
    if ~mod(iTrial,71)
        str = 'Take a short break';
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-strBounds(4)/2-50, [0 0 0]);
        Screen('FillRect', dp.wPtr, [1 1 1], fixationOut);
        Screen('FillRect', dp.wPtr, [0 0 0], fixationIn);
        Screen('Flip', dp.wPtr);
        isResponse = 0;
        while ~isResponse
            [keyIsDown, secs, keyCode] = KbCheck(-1);
            if keyCode(kb.spaceKey)
                isResponse = 1;
            elseif keyCode(kb.escKey)
                sca;
            end
        end
        Screen('FillRect', dp.wPtr, [1 1 1], fixationOut);
        Screen('FillRect', dp.wPtr, [0 0 0], fixationIn);
        Screen('Flip', dp.wPtr);
        startTime = GetSecs;
        randITI = randperm(length(param.iti));
        while (GetSecs - startTime) < param.iti(randITI(1)), ;, end
    end
    
    switch cond.trial(iTrial)
        case 1
            TRIG = SIZE1;
        case 2
            TRIG = SIZE2;
        case 3
            TRIG = SIZE3;
        case 4
            TRIG = SIZE4;
        case 5
            TRIG = SIZE5;
    end
    
    % make texture here....
    textureIndex = Screen('MakeTexture', dp.wPtr, texture(:,:,mod(iTrial-1,2)+1,cond.trial(iTrial)));
    % DRAW TEXTURE
    Screen('DrawTexture', dp.wPtr, textureIndex);
    % draw fixation
    Screen('FillOval', dp.wPtr, [.5 .5 .5],fixationOval);
    Screen('FillRect', dp.wPtr, [1 1 1], fixationOut);
    Screen('FillRect', dp.wPtr, [0 0 0], fixationIn);

    if ~isempty(daq_id), err=DaqDOut(daq_id,0,TRIG); end
    Screen('Flip',dp.wPtr);  %flip!
    startTime = GetSecs;
    while (GetSecs - startTime) < param.trialDur
        [pressed, firstPress, firstRelease, lastPress, lastRelease] = KbQueueCheck();
        if firstPress(kb.escKey)
            sca;
            ShowCursor;
            ListenChar(1);
        end
    end
    elapsedTime(iTrial) = GetSecs-startTime;
    
    %clear out memory
    Screen('Close', textureIndex)

    startTime = GetSecs;

    TRIG = TRIALEND;
    Screen('FillRect', dp.wPtr, [1 1 1], fixationOut);
    Screen('FillRect', dp.wPtr, [0 0 0], fixationIn);
    
    if ~isempty(daq_id), err=DaqDOut(daq_id,0,TRIG); end
    flipTime = Screen('Flip',dp.wPtr);  %flip!

    randITI = randperm(length(param.iti));
    while (GetSecs - startTime) < param.iti(randITI(1)), ;, end
    itiTime(iTrial) = GetSecs-startTime;
end

if ~isempty(daq_id), err=DaqDOut(daq_id,0,END); end

str = 'Thank you!!!';
strBounds = Screen('TextBounds', dp.wPtr, str);
Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-strBounds(4)/2-50, [0 0 0]);
Screen('FillRect', dp.wPtr, [1 1 1], fixationOut);
Screen('FillRect', dp.wPtr, [0 0 0], fixationIn);
Screen('Flip', dp.wPtr);
isResponse = 0;
while ~isResponse
    [keyIsDown, secs, keyCode] = KbCheck(-1);
    if keyCode(kb.spaceKey)
        isResponse = 1;
    elseif keyCode(kb.escKey)
        sca;
    end
end

sca;
ShowCursor;
ListenChar(1);