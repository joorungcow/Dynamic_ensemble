function practice_symmetry(params,center,dp,kb,subID,expVer,iRun)

showInfo(params,center,dp,kb,iRun)
make_fixation(dp,center,params.fixSize,params.fixColor)
Screen('Flip', dp.wPtr);
t1 = GetSecs;
while GetSecs-t1 <= 1; end

%%
tSizes = [params.sSize-0.2 params.sSize-0.1 params.sSize+0.1 params.sSize+0.2];

for trial = 1 : 30
    params.randLoc(trial) = floor(rand(1)+.5);
    params.randSym(trial) = floor(rand(1)+.5);
    params.sSizeSdBig(trial) = floor(rand(1)+.5);
    
    if params.debugFullDots
        params.randSym(trial) = 1;
    end

    if trial <= 15
        duraIndx = 1;
    else
        duraIndx = 2;
    end

    tnum=Sample(1:length(tSizes));
    MeanSize = tSizes(tnum);
    [params, MeanSize] = makeStimuli(params, trial, MeanSize, dp, center, duraIndx);
    
    KbQueueCheck;
    isResponse = 0;
    tstart = GetSecs;

    while ~isResponse
        %[keyIsDown, secs, keyCode] = KbCheck(-1);
        [pressed, firstPress]=KbQueueCheck([]);
        if pressed
            tend = GetSecs;
            if firstPress(kb.leftKey) && ~params.randLoc(trial) || firstPress(kb.rightKey) && params.randLoc(trial) 
                if MeanSize > params.sSize
                    response = 1;
                else
                    response = 0;
                end
                isResponse = 1;  
            elseif firstPress(kb.leftKey) && params.randLoc(trial) || firstPress(kb.rightKey) && ~params.randLoc(trial) 
                if MeanSize < params.sSize
                    response = 1;
                else
                    response = 0;
                end
                isResponse = 1;  
            elseif firstPress(kb.escKey)
                %Screen('LoadNormalizedGammaTable', dp.wPtr, repmat(linspace(0,1, 256)',1,3));
                ListenChar(1);
                ShowCursor;
                sca;
                break
            end
        end
    end

    practiceResult.Location(trial) = params.randLoc(trial); % true: test right visual field
    practiceResult.intensity(trial) = MeanSize;
    practiceResult.response(trial) = response;
    practiceResult.responseTime(trial) = tend - tstart;
    practiceResult.Symmetry(trial) = params.randSym(trial); % true: test dots is symmetry
    practiceResult.sSizeSdBig(trial) = params.sSizeSdBig(trial);
    practiceResult.dotSizes(trial,:,:) = params.dotSizes(trial,:,:);
    practiceResult.duration(trial) = params.duration(duraIndx);

    if response == 1
        make_fixation(dp,center,params.fixSize,[0 .8 0]);
    elseif response == 0
        make_fixation(dp,center,params.fixSize,[.8 0 0]);
    end
    
    Screen('Flip', dp.wPtr);
    t1 = GetSecs;
    while GetSecs-t1 <= .8; end
end

%%
if ~exist('results')
    mkdir('results')
end

cd('results')

if ~exist('practice')
    mkdir('practice')
end

cd('practice')

expVerN = num2str(expVer);
iRunN = num2str(iRun);
fn = sprintf('practice_%s_%s_%s_%s.mat',subID, expVerN, iRunN, datestr(now,'yyyymmdd_HHMMSS'));
save(fn,'practiceResult','params');

str = double('잠시만 기다려주세요.');
strBounds = Screen('TextBounds', dp.wPtr, str);
Screen(dp.wPtr,'TextSize', 30);
Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-120);
make_fixation(dp,center,params.fixSize,params.fixColor)
Screen('Flip', dp.wPtr);
t1 = GetSecs;
while GetSecs-t1 <= 2; end

cd('..')
cd('..')
