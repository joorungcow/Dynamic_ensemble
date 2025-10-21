function run_symmetry(params,center,dp,kb,subID,expVer,iRun)

showInfo(params,center,dp,kb,iRun)
make_fixation(dp,center,params.fixSize,params.fixColor)
Screen('Flip', dp.wPtr);
WaitSecs(1);
%%
PF = @PAL_Logistic; %assumed psychometric function
grain = 201; %grain of posterior, high numbers make method more precise at the cost of RAM and time to compute.
             %Always check posterior after method completes [using e.g., :
             %image(PAL_Scale0to1(PM.pdf)*64)] to check whether appropriate
             %grain and parameter ranges were used.
stimRange = linspace(0.81,1.21,21); %0.81부터 1.21까지 .02씩 (1~1.4로 하니 threshold가 1.03근방으로 나와서 더 범위를 좁힘)
priorAlphaRange = linspace(0.81,1.21,grain); %threshold range
priorBetaRange = linspace(log10(.25),log10(25),grain); %slop range가 무엇인지?
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
%MeanSize = params.initSize;
SymIndex = 1;
RanIndex = 1;
%stair randomizing
stairCond = [repmat([1;2;3;4],params.NumTrials/4,1)];
stairCond = Shuffle(stairCond);
for trial = 1 : params.NumTrials
    params.randLoc(trial) = floor(rand(1)+.5); %test left, right visual field
    %params.randSym(trial) = floor(rand(1)+.5);
    if stairCond(trial) == 1 || stairCond(trial) == 2 % random
        MeanSize = stair(stairCond(trial)).xCurrent;
        params.randSym(trial) = 0;
    else %symmetry
        MeanSize = stair(stairCond(trial)).xCurrent;
        params.randSym(trial) = 1;
    end
    [params, MeanSize] = makeStimuli(params, trial, MeanSize, dp, center);
    
    isResponse = 0;
    while ~isResponse
        [keyIsDown, secs, keyCode] = KbCheck(-1);
        if keyCode(kb.leftKey) && ~params.randLoc(trial) || keyCode(kb.rightKey) && params.randLoc(trial)
            %result.response(trial) = 1;
            %result.intensity(trial) = MeanSize;
            if stair(stairCond(trial)).xCurrent > 1
                response = 1;
                if params.randSym(trial)
                    stair(stairCond(trial)) = PAL_AMPM_updatePM(stair(stairCond(trial)), response);
                    result(stairCond(trial)).trialIndex(SymIndex) = trial;
                    result(stairCond(trial)).Location(SymIndex) = params.randLoc(trial); % true: test right visual field
                    SymIndex = SymIndex + 1;
                else
                    stair(stairCond(trial)) = PAL_AMPM_updatePM(stair(stairCond(trial)), response);
                    result(stairCond(trial)).trialIndex(RanIndex) = trial;
                    result(stairCond(trial)).Location(RanIndex) = params.randLoc(trial); % true: test right visual field
                    RanIndex = RanIndex + 1;
                end
            else
                response = 0;
                if params.randSym(trial)
                    stair(stairCond(trial)) = PAL_AMPM_updatePM(stair(stairCond(trial)), response);
                    result(stairCond(trial)).trialIndex(SymIndex) = trial;
                    result(stairCond(trial)).Location(SymIndex) = params.randLoc(trial); % true: test right visual field
                    SymIndex = SymIndex + 1;
                else
                    stair(stairCond(trial)) = PAL_AMPM_updatePM(stair(stairCond(trial)), response);
                    result(stairCond(trial)).trialIndex(RanIndex) = trial;
                    result(stairCond(trial)).Location(RanIndex) = params.randLoc(trial); % true: test right visual field
                    RanIndex = RanIndex + 1;
                end
            end
            isResponse = 1;
        elseif keyCode(kb.leftKey) && params.randLoc(trial) || keyCode(kb.rightKey) && ~params.randLoc(trial)
            %result.response(trial) = 0;
            %result.intensity(trial) = MeanSize;
            if stair(stairCond(trial)).xCurrent >= 1
                response = 0;
                if params.randSym(trial)
                    stair(stairCond(trial)) = PAL_AMPM_updatePM(stair(stairCond(trial)), response);
                    result(stairCond(trial)).trialIndex(SymIndex) = trial;
                    result(stairCond(trial)).Location(SymIndex) = params.randLoc(trial); % true: test right visual field
                    SymIndex = SymIndex + 1;
                else
                    stair(stairCond(trial)) = PAL_AMPM_updatePM(stair(stairCond(trial)), response);
                    result(stairCond(trial)).trialIndex(RanIndex) = trial;
                    result(stairCond(trial)).Location(RanIndex) = params.randLoc(trial); % true: test right visual field
                    RanIndex = RanIndex + 1;
                end
            else
                response = 1;
                if params.randSym(trial)
                    stair(stairCond(trial)) = PAL_AMPM_updatePM(stair(stairCond(trial)), response);
                    result(stairCond(trial)).trialIndex(SymIndex) = trial;
                    result(stairCond(trial)).Location(SymIndex) = params.randLoc(trial); % true: test right visual field
                    SymIndex = SymIndex + 1;
                else
                    stair(stairCond(trial)) = PAL_AMPM_updatePM(stair(stairCond(trial)), response);
                    result(stairCond(trial)).trialIndex(RanIndex) = trial;
                    result(stairCond(trial)).Location(RanIndex) = params.randLoc(trial); % true: test right visual field
                    RanIndex = RanIndex + 1;
                end
            end
            isResponse = 1;
        elseif keyCode(kb.escKey)
            ListenChar(1);
            ShowCursor;
            sca;
            break
        end
    end
    
    if response == 1
        make_fixation(dp,center,params.fixSize,[0 .8 0]);
    elseif response == 0
        make_fixation(dp,center,params.fixSize,[.8 0 0]);
    end
    
    Screen('Flip', dp.wPtr);
    WaitSecs(.8);
end

for i = 1 : params.nCon
    result(i).threshold = stair(i).threshold;
    result(i).intensity = stair(i).x;
    result(i).response = stair(i).response;
    result(i).slope = stair(i).slope;
    result(i).stimRange = stair(i).stimRange;
end

if ~exist('results')
    mkdir('results')
end

cd('results')
expVerN = num2str(expVer);
fn = sprintf('%s_%s_%s.mat',subID, expVerN, datestr(now,'yyyymmdd_HHMMSS'));
save(fn,'result','params');

str = double('잠시만 기다려주세요.');
strBounds = Screen('TextBounds', dp.wPtr, str);
Screen(dp.wPtr,'TextSize', 30);
Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-120);
make_fixation(dp,center,params.fixSize,params.fixColor)
Screen('Flip', dp.wPtr);
WaitSecs(2);

cd('..')
