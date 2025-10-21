function run_symmetry(params,center,dp,kb,subID,expVer,iRun)

showInfo(params,center,dp,kb,iRun)
make_fixation(dp,center,params.fixSize,params.fixColor)
Screen('Flip', dp.wPtr);
t1 = GetSecs;
while GetSecs-t1 <= 1; end
%%
PF = @PAL_Logistic; %assumed psychometric function
grain = 201; %grain of posterior, high numbers make method more precise at the cost of RAM and time to compute.
             %Always check posterior after method completes [using e.g., :
             %image(PAL_Scale0to1(PM.pdf)*64)] to check whether appropriate
             %grain and parameter ranges were used.
stimRange = linspace(params.sSize-0.2,params.sSize+0.2,21); 
priorAlphaRange = linspace(params.sSize-0.2,params.sSize+0.2,grain); %threshold range
priorBetaRange = linspace(log10(.25),log10(25),grain); 
priorGammaRange = 0;  %baseline
priorLambdaRange = .02; %lapsing range
%Initialize PM structure

for i = 1 : params.nCon % 1,2: random 3.4: symmetry
    stair(i) = PAL_AMPM_setupPM('priorAlphaRange',priorAlphaRange,...
        'priorBetaRange',priorBetaRange,...
        'priorGammaRange',priorGammaRange,...
        'priorLambdaRange',priorLambdaRange,...
        'numtrials',params.NumTrials,...
        'PF' , PF,...
        'stimRange',stimRange);
end

%%
%indexing
for i = 1 : 4
    Index(i) = 1;
end

duraIndx = iRun+2;

%stair randomizing
stairCond = [repmat([1;2;3;4],params.NumTrials/4,1)];
stairCond = Shuffle(stairCond);

%variance randomizing
params.sSizeSdBig = [repmat([0;1],params.NumTrials/2,1)];
params.sSizeSdBig = Shuffle(params.sSizeSdBig);

%left, right visual field randomizing
params.randLoc = [repmat([0;1],params.NumTrials/2,1)];
params.randLoc = Shuffle(params.randLoc);

for trial = 1 : params.NumTrials
    if stairCond(trial) == 1 || stairCond(trial) == 2 % random
        MeanSize = stair(stairCond(trial)).xCurrent;
        params.randSym(trial) = 0;
        if params.debugFullDots
            params.randSym(trial) = 1;
        end

    else %symmetry
        MeanSize = stair(stairCond(trial)).xCurrent;
        params.randSym(trial) = 1;
    end
    [params] = makeStimuli(params, trial, MeanSize, dp, center, duraIndx);
    
    KbQueueCheck;
    isResponse = 0;
    tstart = GetSecs;

    while ~isResponse
        %[keyIsDown, secs, keyCode] = KbCheck(-1);
        [pressed, firstPress]=KbQueueCheck([]);
        if pressed
            if firstPress(kb.leftKey) && ~params.randLoc(trial) || firstPress(kb.rightKey) && params.randLoc(trial)
                response = 1; % 1:test is bigger 0:standard is bigger
                tend = GetSecs;
                isResponse = 1;  
            elseif firstPress(kb.leftKey) && params.randLoc(trial) || firstPress(kb.rightKey) && ~params.randLoc(trial)
                response = 0; % 1:test is bigger 0:standard is bigger
                tend = GetSecs;
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
    
    stair(stairCond(trial)) = PAL_AMPM_updatePM(stair(stairCond(trial)), response);
    result(stairCond(trial)).trialIndex(Index(stairCond(trial))) = trial;
    result(stairCond(trial)).Location(Index(stairCond(trial))) = params.randLoc(trial); % true: test right visual field            
    result(stairCond(trial)).responseTime(Index(stairCond(trial))) = tend - tstart;
    result(stairCond(trial)).Symmetry(Index(stairCond(trial))) = params.randSym(trial); % true: test dots is symmetry
    result(stairCond(trial)).sSizeSdBig(Index(stairCond(trial))) = params.sSizeSdBig(trial); 
    Index(stairCond(trial)) = Index(stairCond(trial)) + 1;
            
    make_fixation(dp,center,params.fixSize,[0 0 0]);
    
    Screen('Flip', dp.wPtr);
    t1 = GetSecs;
    while GetSecs-t1 <= .8; end
end

for i = 1 : params.nCon
    result(i).threshold = stair(i).threshold;
    result(i).intensity = stair(i).x;
    result(i).response = stair(i).response;
    result(i).slope = stair(i).slope;
    result(i).stimRange = stair(i).stimRange;
    result(i).duration = params.duration(duraIndx);
end

if ~exist('results')
    mkdir('results')
end

cd('results')
expVerN = num2str(expVer);
iRunN = num2str(iRun);
fn = sprintf('%s_%s_%s_%s.mat',subID, expVerN, iRunN, datestr(now,'yyyymmdd_HHMMSS'));
save(fn,'result','params');

str = double('잠시만 기다려주세요.');
strBounds = Screen('TextBounds', dp.wPtr, str);
Screen(dp.wPtr,'TextSize', 30);
Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-120);
make_fixation(dp,center,params.fixSize,params.fixColor)
Screen('Flip', dp.wPtr);
t1 = GetSecs;
while GetSecs-t1 <= 2; end

cd('..')
