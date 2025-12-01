AssertOpenGL;
KbName('UnifyKeyNames');

subID = strtrim(input('Enter subject ID (e.g., sub01): ', 's'));
if isempty(subID)
    error('Subject ID is required.');
end

gender = upper(strtrim(input('성별을 입력하세요 (M/F): ', 's')));
ageStr = strtrim(input('만 나이를 숫자로 입력하세요: ', 's'));
age    = str2double(ageStr);
hand   = upper(strtrim(input('주로 사용하는 손을 입력하세요 (R/L): ', 's')));
colorVision = upper(strtrim(input('시력 이상 또는 색약/색맹 여부 Y / N: ', 's')));

planFile = fullfile('sampling', sprintf('%s_sampling.mat', subID));
if ~exist(planFile, 'file')
    error('Plan file %s not found. Run sampling2.m first.', planFile);
end

planData = load(planFile, 'plan');
if ~isfield(planData, 'plan')
    error('Plan file %s does not contain a plan struct.', planFile);
end
plan = planData.plan;

if ~isfield(plan, 'display') || ~isstruct(plan.display)
    error('Plan file %s is missing display settings. Regenerate using sampling2.m.', planFile);
end
dp = plan.display;
if ~isfield(dp, 'screenNum') || isempty(dp.screenNum)
    dp.screenNum = max(Screen('Screens'));
end

if ~isfield(plan, 'dotParams') || ~isstruct(plan.dotParams)
    error('Plan file %s is missing dot configuration. Regenerate using sampling2.m.', planFile);
end
params = plan.dotParams;
if ~isfield(params, 'perceptualExponent') || isempty(params.perceptualExponent)
    params.perceptualExponent = plan.meta.expo;
end
if ~isfield(params, 'meanDiffLevels') || isempty(params.meanDiffLevels)
    params.meanDiffLevels = plan.meta.diffLevels;
end

if ~isfield(plan, 'motionParams') || ~isstruct(plan.motionParams)
    error('Plan file %s is missing motion configuration. Regenerate using sampling2.m.', planFile);
end
motionParams = plan.motionParams;

if ~isfield(plan, 'timingParams') || ~isstruct(plan.timingParams)
    error('Plan file %s is missing timing configuration. Regenerate using sampling2.m.', planFile);
end
timingParams = plan.timingParams;

if ~isfield(plan, 'gridConfig') || ~isstruct(plan.gridConfig)
    error('Plan file %s is missing grid configuration. Regenerate using sampling2.m.', planFile);
end
gridConfig = plan.gridConfig;

% gridLayout: use from plan if present, else rebuild
if isfield(plan, 'gridLayout') && ~isempty(plan.gridLayout)
    gridLayout = plan.gridLayout;
else
    gridLayout = buildCentralGridLayout(gridConfig);
end

if isfield(timingParams, 'breakDurationSec') && ~isempty(timingParams.breakDurationSec)
    breakDurationSec = timingParams.breakDurationSec;
else
    breakDurationSec = 20;
end

%% ----------------- Keyboard & results struct -----------------
kbConfig = struct();
kbConfig.useKbQueueCheck = 0;
kb = init_keyboard(kbConfig);

results = struct();
results.meta = plan.meta;
results.meta.subID       = subID;
results.meta.generatedAt = datestr(now, 'yyyymmdd_HHMMSS');
results.meta.gender      = gender;
results.meta.age         = age;
results.meta.hand        = hand;
results.meta.colorVision = colorVision;
results.meta.display      = plan.display;
results.meta.dotParams    = plan.dotParams;
results.meta.motionParams = plan.motionParams;
results.meta.timingParams = plan.timingParams;
results.meta.gridConfig   = plan.gridConfig;

saveFileName = sprintf('results_%s_%s.mat', subID);

%% ----------------- Main try/catch -----------------
try
    dp = OpenWindow(dp);
    HideCursor;
    ListenChar(2);

    epsDeg = 1 / dp.ppd;
    fprintf('Estimated pixel pitch: %.4f°/pixel.\n', epsDeg);
    if params.gToleranceDeg < 0.5 * epsDeg
        fprintf(['[경고] gTolerance가 디스플레이 상도(%.4f°)보다 매우 작습니다. ' ...
            '계산상 문제는 없으나 물리적으로 구분이 어려울 수 있습니다.\n'], epsDeg);
    end

    Screen('FillRect', dp.wPtr, dp.bkColor);
    Screen('Flip', dp.wPtr);

    % --------- Start screen: choose practice vs main ---------
    mode = showStartScreen(dp, kb); % 'main', 'practice', 'abort'
    if strcmp(mode, 'abort')
        error('Experiment aborted at start screen.');
    end

    % --------- Practice blocks (optional) ---------
    if strcmp(mode, 'practice')
        if isfield(plan, 'practiceBlocks') && ~isempty(plan.practiceBlocks)
            fprintf('Running practice blocks...\n');
            results.practiceTrials = runPracticeBlocks(dp, kb, plan, params, gridLayout, gridConfig, motionParams, timingParams);
            showPracticeEndScreen(dp, kb);
        else
            fprintf('No practiceBlocks found in plan. Skipping practice.\n');
        end
    end

    % --------- Main experiment ---------
    if isfield(plan.meta, 'NMain') && ~isempty(plan.meta.NMain)
        totalTrials = plan.meta.NMain;
    else
        totalTrials = plan.meta.N;
    end
    results.trials(totalTrials, 1) = initTrialResultStruct();

    trialIdx        = 1;
    abortExperiment = false;

    for b = 1:numel(plan.blocks)
        block = plan.blocks(b);
        showBlockIntro(dp, sprintf('블록 %d / %d (%s)', b, numel(plan.blocks), block.label));

        for t = 1:numel(block.trials)
            trialSpec = block.trials(t);
            if abortExperiment
                break;
            end

            [stimulusPair, stimInfo] = buildStimulusFromPlan(trialSpec, params, gridLayout, gridConfig);

            if timingParams.fixationMs > 0
                abortExperiment = presentFixation(dp, timingParams.fixationMs, kb); % default white
                if abortExperiment
                    break;
                end
            end

            abortExperiment = presentStimulus(dp, stimulusPair.t1, trialSpec.isMoving(1), motionParams, timingParams.stimDurationMs, kb);
            if abortExperiment, break; end

            abortExperiment = presentBlank(dp, timingParams.isiDurationMs, kb);
            if abortExperiment, break; end

            abortExperiment = presentStimulus(dp, stimulusPair.t2, trialSpec.isMoving(2), motionParams, timingParams.stimDurationMs, kb);
            if abortExperiment, break; end

            response = collectResponse(dp, kb, timingParams.postTrialMs);
            if response.wasAborted
                abortExperiment = true;
                break;
            end

            trialResult = populateTrialResult(trialIdx, trialSpec, stimInfo, response, params);
            trialResult.isPractice = false;
            results.trials(trialIdx) = trialResult;
            trialIdx = trialIdx + 1;
        end

        if abortExperiment
            break;
        end

        if b < numel(plan.blocks)
            abortExperiment = presentBreakScreen(dp, kb, breakDurationSec);
            if abortExperiment
                break;
            end
        end
    end

    completed = trialIdx - 1;
    results.trials = results.trials(1:completed);

    presentBlank(dp, 0, kb);

    Screen('CloseAll');
    ListenChar(0);
    ShowCursor;
    if kb.useKbQueueCheck
        KbQueueRelease;
    end
    RestrictKeysForKbCheck([]);

    save(saveFileName, 'results');
    fprintf('Results saved to %s\n', saveFileName);

catch ME
    Screen('CloseAll');
    ListenChar(0);
    ShowCursor;
    if exist('kb','var') && isfield(kb,'useKbQueueCheck') && kb.useKbQueueCheck
        KbQueueRelease;
    end
    RestrictKeysForKbCheck([]);
    fprintf('Error occurred: %s\n', ME.message);
    crashFile = sprintf('results_%s_%s_crash.mat', subID, datestr(now, 'yyyymmdd_HHMMSS'));
    if exist('results','var')
        save(crashFile, 'results');
        fprintf('Partial results saved to %s\n', crashFile);
    end
    rethrow(ME);
end


%% =================== Helper functions ===================================


function resultStruct = initTrialResultStruct()
resultStruct = struct( ...
    'trialIndex', [], ...
    'comboLabel', '', ...
    'methodPair', '', ...
    'freqRatio', '', ...
    'equalFreq', true, ...
    'diffLevelTarget', [], ...
    'signS', [], ...
    'largerSide', '', ...
    'isMoving', [], ...
    't1MeanPS', [], ...
    't2MeanPS', [], ...
    'psDiffObserved', [], ...
    'diffLevelObservedPS', [], ...
    't1MeanDeg', [], ...
    't2MeanDeg', [], ...
    'expoUsed', [], ...
    'jitterPS_applied', [], ...
    'responseKey', '', ...
    'responseChoice', '', ...
    'responseRtMs', NaN, ...
    'didRespond', false, ...
    'correct', NaN, ...
    'targetPSMeans', [], ...
    'methodPairLabel', '', ...
    'comboMotion', '', ...
    'psTargetsAfterJitter', [], ...
    'isPractice', false); 
end

function showBlockIntro(dp, labelText)
if nargin < 2 || isempty(labelText)
    labelText = '';
end
Screen('FillRect', dp.wPtr, dp.bkColor);
if isfield(dp, 'textFont') && ~isempty(dp.textFont)
    Screen('TextFont', dp.wPtr, dp.textFont);
end
if isfield(dp, 'textSize') && ~isempty(dp.textSize)
    Screen('TextSize', dp.wPtr, dp.textSize);
end
if ~isempty(labelText)
    Screen('DrawText', dp.wPtr, double(labelText), dp.cx - 200, dp.cy, dp.textColor);
end
Screen('Flip', dp.wPtr);
WaitSecs(1.0);
end

function mode = showStartScreen(dp, kb)
% mode: 'main', 'practice', 'abort'

if kb.useKbQueueCheck
    KbQueueFlush;
else
    KbReleaseWait;
end

lines = {
    '이 실험은 총 3 block으로 진행됩니다.';
    ' ';
    '각 block이 끝날 때마다 실험자를 불러';
    ' ';
    '다음 block 설명을 듣고, 실험을 진행해주세요.';
    ' ';
    '스페이스바를 누르면 본 실험을 시작합니다.';
    ' ';
    'P 키를 누르면 연습 시행을 시작합니다.';
    ' ';
};

Screen('FillRect', dp.wPtr, dp.bkColor);
if isfield(dp, 'textFont') && ~isempty(dp.textFont)
    Screen('TextFont', dp.wPtr, dp.textFont);
end
if isfield(dp, 'textSize') && ~isempty(dp.textSize)
    Screen('TextSize', dp.wPtr, dp.textSize);
end

firstFrame = true;
vbl = 0;

while true
    Screen('FillRect', dp.wPtr, dp.bkColor);
    y = dp.cy - numel(lines)*15;
    for i = 1:numel(lines)
        bounds = Screen('TextBounds', dp.wPtr, double(lines{i}));
        x = dp.cx - RectWidth(bounds)/2;
        Screen('DrawText', dp.wPtr, double(lines{i}), x, y, dp.textColor);
        y = y + 30;
    end

    if firstFrame
        vbl = Screen('Flip', dp.wPtr);
        firstFrame = false;
    else
        vbl = Screen('Flip', dp.wPtr, vbl + 0.5*dp.ifi);
    end

    if kb.useKbQueueCheck
        [pressed, firstPress] = KbQueueCheck;
        if pressed
            if any(firstPress(kb.escKey)),   mode = 'abort';    return; end
            if any(firstPress(kb.spaceKey)), mode = 'main';     return; end
            if any(firstPress(kb.pKey)),     mode = 'practice'; return; end
        end
    else
        [keyIsDown, ~, keyCode] = KbCheck(-1);
        if keyIsDown
            if keyCode(kb.escKey),   mode = 'abort';    return; end
            if keyCode(kb.spaceKey), mode = 'main';     return; end
            if keyCode(kb.pKey),     mode = 'practice'; return; end
        end
    end
end
end

function abort = presentFixation(dp, durationMs, kb, fixColor)
% Colored fixation cross; default white if fixColor omitted.
if nargin < 4 || isempty(fixColor)
    fixColor = [1 1 1];
end

abort = false;
numFrames = max(1, round((durationMs/1000) / dp.ifi));
firstFrame = true;
vbl = 0;

for frameIdx = 1:numFrames %#ok<NASGU>
    Screen('FillRect', dp.wPtr, dp.bkColor);
    % draw fixation here (color version)
    fixSize   = 20;
    lineWidth = 2;
    Screen('DrawLine', dp.wPtr, fixColor, dp.cx-fixSize/2, dp.cy, dp.cx+fixSize/2, dp.cy, lineWidth);
    Screen('DrawLine', dp.wPtr, fixColor, dp.cx, dp.cy-fixSize/2, dp.cx, dp.cy+fixSize/2, lineWidth);

    if firstFrame
        vbl = Screen('Flip', dp.wPtr);
        firstFrame = false;
    else
        vbl = Screen('Flip', dp.wPtr, vbl + 0.5 * dp.ifi);
    end
    if shouldAbort(kb)
        abort = true;
        break;
    end
end
end

function practiceResults = runPracticeBlocks(dp, kb, plan, params, gridLayout, gridConfig, motionParams, timingParams)

practiceBlocks = plan.practiceBlocks;
nPracticeTotal = sum(arrayfun(@(b) numel(b.trials), practiceBlocks));

practiceResults(nPracticeTotal,1) = initTrialResultStruct();

trialIdx        = 1;
abortExperiment = false;

fixColorNext = [1 1 1]; % first fixation: white

for b = 1:numel(practiceBlocks)
    block = practiceBlocks(b);
    showBlockIntro(dp, sprintf('연습 블록 %d / %d (%s)', b, numel(practiceBlocks), block.label));

    for t = 1:numel(block.trials)
        if abortExperiment, break; end

        trialSpec = block.trials(t);
        [stimulusPair, stimInfo] = buildStimulusFromPlan(trialSpec, params, gridLayout, gridConfig);

        % feedback fixation (color based on previous trial correctness)
        if timingParams.fixationMs > 0
            abortExperiment = presentFixation(dp, timingParams.fixationMs, kb, fixColorNext);
            if abortExperiment, break; end
        end

        abortExperiment = presentStimulus(dp, stimulusPair.t1, trialSpec.isMoving(1), motionParams, timingParams.stimDurationMs, kb);
        if abortExperiment, break; end

        abortExperiment = presentBlank(dp, timingParams.isiDurationMs, kb);
        if abortExperiment, break; end

        abortExperiment = presentStimulus(dp, stimulusPair.t2, trialSpec.isMoving(2), motionParams, timingParams.stimDurationMs, kb);
        if abortExperiment, break; end

        response = collectResponse(dp, kb, timingParams.postTrialMs);
        if response.wasAborted
            abortExperiment = true;
            break;
        end

        tr = populateTrialResult(trialIdx, trialSpec, stimInfo, response, params);
        tr.isPractice = true;
        practiceResults(trialIdx) = tr;

        % Decide next fixation color by correctness
        if tr.didRespond && ~isnan(tr.correct)
            if tr.correct == 1
                fixColorNext = [0 1 0]; % green
            else
                fixColorNext = [1 0 0]; % red
            end
        else
            fixColorNext = [1 1 1];     % white for miss/NaN
        end

        trialIdx = trialIdx + 1;
    end

    if abortExperiment, break; end
end

practiceResults = practiceResults(1:trialIdx-1);
end

function showPracticeEndScreen(dp, kb)
msg = '연습 시행이 끝났습니다. 스페이스바를 누르면 본 실험을 시작합니다.';

Screen('FillRect', dp.wPtr, dp.bkColor);
if isfield(dp, 'textFont') && ~isempty(dp.textFont)
    Screen('TextFont', dp.wPtr, dp.textFont);
end
if isfield(dp, 'textSize') && ~isempty(dp.textSize)
    Screen('TextSize', dp.wPtr, dp.textSize);
end

bounds = Screen('TextBounds', dp.wPtr, double(msg));
x = dp.cx - RectWidth(bounds)/2;
y = dp.cy;

Screen('DrawText', dp.wPtr, double(msg), x, y, dp.textColor);
Screen('Flip', dp.wPtr);

while true
    if kb.useKbQueueCheck
        [pressed, firstPress] = KbQueueCheck;
        if pressed && any(firstPress(kb.spaceKey)), return; end
    else
        [keyIsDown, ~, keyCode] = KbCheck(-1);
        if keyIsDown && keyCode(kb.spaceKey), return; end
    end
end
end

function [stimPair, info] = buildStimulusFromPlan(trialSpec, params, layout, gridConfig)
expo = params.perceptualExponent;
freqCounts = getFrequencyRatioCounts(trialSpec.freqRatio);
eqCounts   = params.equalFreqCounts;

freqBase = make_frequency_set(freqCounts, params);
eqBase   = make_equalFreq_set(eqCounts, params);

if strcmp(trialSpec.methodPair, 'FREQ_vs_EQ')
    baseT1 = freqBase;
    baseT2 = eqBase;
else
    baseT1 = eqBase;
    baseT2 = freqBase;
end

ratioTarget = 1 + trialSpec.signS * trialSpec.diffLevel;
mu1Target   = meanPS(baseT1, expo);
mu2Target   = mu1Target * ratioTarget;

jitterRangePS = computeJitterRangePS(mu1Target, mu2Target, params);
[mu1Target, mu2Target] = jitterMeanPairPS(mu1Target, mu2Target, jitterRangePS, trialSpec.diffLevel, trialSpec.signS);
mu1TargetJittered = mu1Target;
mu2TargetJittered = mu2Target;

[t1SizesDeg, t1MeanPS] = synthesizeSetFromBase(baseT1, mu1Target, params);
ratioTarget = 1 + trialSpec.signS * trialSpec.diffLevel;
mu2Target   = t1MeanPS * ratioTarget;
[t2SizesDeg, t2MeanPS] = synthesizeSetFromBase(baseT2, mu2Target, params);

stimPair = struct();
stimPair.t1 = createStimulusFromSizes(t1SizesDeg, layout, gridConfig);
stimPair.t2 = createStimulusFromSizes(t2SizesDeg, layout, gridConfig);

info = struct();
info.t1MeanPS = t1MeanPS;
info.t2MeanPS = t2MeanPS;
info.t1MeanDeg = mean(t1SizesDeg);
info.t2MeanDeg = mean(t2SizesDeg);
info.psDiffObserved = (t2MeanPS / t1MeanPS) - 1;
info.diffLevelObservedPS = info.psDiffObserved * 100;
info.expoUsed  = expo;
info.jitterPS  = jitterRangePS;
info.targetPSMeans        = [mu1TargetJittered, mu2TargetJittered];
info.psTargetsAfterJitter = [mu1TargetJittered, mu2Target];
info.largerSide           = ternary(trialSpec.signS >= 0, 'T2', 'T1');
end

function [sizesDeg, finalMeanPS] = synthesizeSetFromBase(baseSizesDeg, targetMeanPS, params)
% test13.m의 synthesizeSetFromBase 복사
expo   = params.perceptualExponent;
basePS = toPS(baseSizesDeg, expo);
scale  = targetMeanPS / mean(basePS);
basePS = basePS * scale;

minPS = toPS(params.minSizeDeg, expo);
maxPS = toPS(params.maxSizeDeg, expo);
tolerancePS = computePSTolerance(params.gToleranceDeg, targetMeanPS, expo);
if tolerancePS <= 0
    tolerancePS = 1e-5;
end

maxAttempts = 1000;
noiseStdPS  = params.jitterStdRatio * targetMeanPS;
bestPS      = basePS;
bestDiff    = abs(mean(basePS) - targetMeanPS);

for attempt = 1:maxAttempts
    if attempt == round(maxAttempts * 0.5)
        noiseStdPS = noiseStdPS * 0.5;
    elseif attempt == round(maxAttempts * 0.75)
        noiseStdPS = noiseStdPS * 0.25;
    end

    noise = noiseStdPS .* randn(size(basePS));
    noise = noise - mean(noise);

    candidatePS = basePS + noise;
    candidatePS = min(max(candidatePS, minPS), maxPS);
    candidatePS = adjustMeanPSWithinBounds(candidatePS, targetMeanPS, minPS, maxPS, tolerancePS);

    candidateDeg = fromPS(candidatePS, expo);
    candidateDeg = min(max(candidateDeg, params.minSizeDeg), params.maxSizeDeg);
    candidatePS  = toPS(candidateDeg, expo);
    candidatePS  = adjustMeanPSWithinBounds(candidatePS, targetMeanPS, minPS, maxPS, tolerancePS);

    currentDiff = abs(mean(candidatePS) - targetMeanPS);
    if currentDiff < bestDiff
        bestDiff = currentDiff;
        bestPS   = candidatePS;
    end

    if currentDiff <= tolerancePS
        bestPS = candidatePS;
        break;
    end
end

finalPS   = bestPS;
sizesDeg  = fromPS(finalPS, expo);
sizesDeg  = min(max(sizesDeg, params.minSizeDeg), params.maxSizeDeg);
finalPS   = toPS(sizesDeg, expo);
finalMeanPS = mean(finalPS);
end

function tolerancePS = computePSTolerance(tolDeg, targetMeanPS, expo)
targetDeg = fromPS(targetMeanPS, expo);
tolerancePS = abs(toPS(targetDeg + tolDeg, expo) - targetMeanPS);
end

function projected = adjustMeanPSWithinBounds(psValues, targetMeanPS, minPS, maxPS, tolerance)
numVals   = numel(psValues);
maxIter   = 50;
projected = psValues;
for iter = 1:maxIter
    meanDiff = targetMeanPS - mean(projected);
    if abs(meanDiff) <= tolerance
        break;
    end

    totalDiff = meanDiff * numVals;
    if totalDiff > 0
        adjustable = find(projected < maxPS - eps);
        if isempty(adjustable), break; end
        slack      = maxPS - projected(adjustable);
        totalSlack = sum(slack);
        if totalSlack <= 0, break; end
        factor     = min(1, totalDiff / totalSlack);
        increments = factor .* slack;
        projected(adjustable) = projected(adjustable) + increments;
    else
        adjustable = find(projected > minPS + eps);
        if isempty(adjustable), break; end
        slack      = projected(adjustable) - minPS;
        totalSlack = sum(slack);
        if totalSlack <= 0, break; end
        factor     = min(1, -totalDiff / totalSlack);
        decrements = factor .* slack;
        projected(adjustable) = projected(adjustable) - decrements;
    end
    projected = min(max(projected, minPS), maxPS);
end

if abs(mean(projected) - targetMeanPS) > tolerance
    offset   = targetMeanPS - mean(projected);
    projected = projected + offset;
    projected = min(max(projected, minPS), maxPS);
end
end

function jitterPS = computeJitterRangePS(mu1PS, mu2PS, params)
expo    = params.perceptualExponent;
meanDeg = mean([fromPS(mu1PS, expo), fromPS(mu2PS, expo)]);
meanDeg = max(meanDeg, params.minSizeDeg);
jitterPS = abs(toPS(meanDeg + params.meanJitterDeg, expo) - toPS(meanDeg, expo));
end

function [mu1Adj, mu2Adj] = jitterMeanPairPS(mu1PS, mu2PS, jitterRangePS, diffLevel, signS)
if jitterRangePS <= 0
    mu1Adj = mu1PS;
    mu2Adj = mu2PS;
    return;
end

jitter1 = (rand * 2 - 1) * jitterRangePS;
jitter2 = (rand * 2 - 1) * jitterRangePS;
mu1Adj  = max(eps, mu1PS + jitter1);
mu2Adj  = max(eps, mu2PS + jitter2);

ratioTarget  = 1 + signS * diffLevel;
projectedMu2 = mu1Adj * ratioTarget;
projectedMu1 = mu2Adj / ratioTarget;
if abs(projectedMu2 - mu2Adj) <= abs(projectedMu1 - mu1Adj)
    mu2Adj = max(eps, projectedMu2);
else
    mu1Adj = max(eps, projectedMu1);
end
end

function sizes = make_frequency_set(counts, params)
sizes = [repmat(params.smallSizeDeg, 1, counts(1)), ...
         repmat(params.largeSizeDeg, 1, counts(2))];
sizes = sizes(randperm(numel(sizes)));
end

function sizes = make_equalFreq_set(counts, params)
base  = [repmat(params.smallSizeDeg, 1, counts(1)), ...
         repmat(params.largeSizeDeg, 1, counts(2))];
sizes = base(randperm(numel(base)));
end

function counts = getFrequencyRatioCounts(label)
switch label
    case '6:2'
        counts = [6 2];
    case '5:3'
        counts = [5 3];
    case '3:5'
        counts = [3 5];
    case '2:6'
        counts = [2 6];
    otherwise
        error('Unknown frequency ratio label: %s', label);
end
end

function stimSet = createStimulusFromSizes(dotSizesDeg, layout, gridConfig)
stimSet.dotSizeDeg = dotSizesDeg;
[stimSet.xPosDeg, stimSet.yPosDeg] = placeDotsOnGrid(dotSizesDeg, layout, gridConfig);

stimSet.outerAperture = layout.outerAperture;
stimSet.innerAperture = layout.innerAperture;
stimSet.halfSizeDeg   = stimSet.dotSizeDeg / 2;

stimSet.topEdgesDeg    = stimSet.outerAperture.topDeg    + stimSet.halfSizeDeg;
stimSet.bottomEdgesDeg = stimSet.outerAperture.bottomDeg - stimSet.halfSizeDeg;
stimSet.leftEdgesDeg   = stimSet.outerAperture.leftDeg   + stimSet.halfSizeDeg;
stimSet.rightEdgesDeg  = stimSet.outerAperture.rightDeg  - stimSet.halfSizeDeg;

stimSet.verticalSpanDeg   = stimSet.bottomEdgesDeg   - stimSet.topEdgesDeg;
stimSet.horizontalSpanDeg = stimSet.rightEdgesDeg    - stimSet.leftEdgesDeg;

stimSet.innerTopEdgesDeg    = stimSet.innerAperture.topDeg    + stimSet.halfSizeDeg;
stimSet.innerBottomEdgesDeg = stimSet.innerAperture.bottomDeg - stimSet.halfSizeDeg;
stimSet.innerLeftEdgesDeg   = stimSet.innerAperture.leftDeg   + stimSet.halfSizeDeg;
stimSet.innerRightEdgesDeg  = stimSet.innerAperture.rightDeg  - stimSet.halfSizeDeg;

if any(stimSet.innerBottomEdgesDeg <= stimSet.innerTopEdgesDeg) || ...
   any(stimSet.innerRightEdgesDeg  <= stimSet.innerLeftEdgesDeg)
    error('Inner aperture must be larger than dot diameters.');
end

if any(stimSet.verticalSpanDeg <= 0) || any(stimSet.horizontalSpanDeg <= 0)
    error('Outer aperture size must exceed the diameter of every dot.');
end
end

function [xPosDeg, yPosDeg] = placeDotsOnGrid(dotSizesDeg, layout, gridConfig)
numDots       = numel(dotSizesDeg);
occupiedCells = false(numel(layout.gridX), 1);

xPosDeg = zeros(1, numDots);
yPosDeg = zeros(1, numDots);

for dotIdx = 1:numDots
    placed    = false;
    cellOrder = randperm(numel(occupiedCells));
    currentSize = dotSizesDeg(dotIdx);
    halfSize    = currentSize / 2;

    cellJitterX = layout.cellWidthDeg  / 2 - (halfSize + gridConfig.safetyMarginDeg);
    cellJitterY = layout.cellHeightDeg / 2 - (halfSize + gridConfig.safetyMarginDeg);

    if cellJitterX <= 0 || cellJitterY <= 0
        error('Cell dimensions are too small for dot %d.', dotIdx);
    end

    for candidate = cellOrder
        if occupiedCells(candidate), continue; end

        baseX = layout.gridX(candidate);
        baseY = layout.gridY(candidate);

        leftLimit   = baseX - (layout.innerAperture.leftDeg  + halfSize + gridConfig.safetyMarginDeg);
        rightLimit  = (layout.innerAperture.rightDeg - halfSize - gridConfig.safetyMarginDeg) - baseX;
        topLimit    = baseY - (layout.innerAperture.topDeg   + halfSize + gridConfig.safetyMarginDeg);
        bottomLimit = (layout.innerAperture.bottomDeg - halfSize - gridConfig.safetyMarginDeg) - baseY;

        maxJitterX = min([gridConfig.maxJitterDeg, cellJitterX, leftLimit, rightLimit]);
        maxJitterY = min([gridConfig.maxJitterDeg, cellJitterY, topLimit, bottomLimit]);

        if maxJitterX <= 0 || maxJitterY <= 0
            continue;
        end

        for attempt = 1:gridConfig.maxAttemptsPerCell
            jitterX = (rand * 2 - 1) * maxJitterX;
            jitterY = (rand * 2 - 1) * maxJitterY;

            candidateX = baseX + jitterX;
            candidateY = baseY + jitterY;

            if candidateX - halfSize < layout.innerAperture.leftDeg  + gridConfig.safetyMarginDeg, continue; end
            if candidateX + halfSize > layout.innerAperture.rightDeg - gridConfig.safetyMarginDeg, continue; end
            if candidateY - halfSize < layout.innerAperture.topDeg   + gridConfig.safetyMarginDeg, continue; end
            if candidateY + halfSize > layout.innerAperture.bottomDeg- gridConfig.safetyMarginDeg, continue; end

            xPosDeg(dotIdx) = candidateX;
            yPosDeg(dotIdx) = candidateY;
            occupiedCells(candidate) = true;
            placed = true;
            break;
        end

        if placed, break; end
    end

    if ~placed
        error('Could not place dot %d within constraints.', dotIdx);
    end
end
end

function abort = presentStimulus(dp, stim, isMoving, motionParams, stimDurationMs, kb)
abort = false;
numFrames = max(1, round((stimDurationMs/1000) / dp.ifi));

sizePix = stim.dotSizeDeg * dp.ppd;

if isMoving
    currentXY = [stim.xPosDeg; stim.yPosDeg];
else
    staticXY = [stim.xPosDeg; stim.yPosDeg] * dp.ppd;
    staticXY(1, :) = staticXY(1, :) + dp.cx;
    staticXY(2, :) = staticXY(2, :) + dp.cy;
end

Screen('FillRect', dp.wPtr, dp.bkColor);
vbl = Screen('Flip', dp.wPtr);

for frameIdx = 1:numFrames
    if isMoving
        if frameIdx > 1
            [stim.xPosDeg, stim.yPosDeg] = updatePositions(stim.xPosDeg, stim.yPosDeg, motionParams, stim, dp);
            currentXY = [stim.xPosDeg; stim.yPosDeg];
        end
        xyPix = currentXY * dp.ppd;
        xyPix(1, :) = xyPix(1, :) + dp.cx;
        xyPix(2, :) = xyPix(2, :) + dp.cy;
    else
        xyPix = staticXY;
    end

    Screen('FillRect', dp.wPtr, dp.bkColor);
    Screen('DrawDots', dp.wPtr, xyPix, sizePix, [1 1 1], [0 0], 2);

    vbl = Screen('Flip', dp.wPtr, vbl + 0.5 * dp.ifi);

    if shouldAbort(kb)
        abort = true;
        break;
    end
end
end

function abort = presentBlank(dp, durationMs, kb)
abort = false;
Screen('FillRect', dp.wPtr, dp.bkColor);

if durationMs <= 0
    Screen('Flip', dp.wPtr);
    abort = shouldAbort(kb);
    return;
end

numFrames = max(1, round((durationMs/1000) / dp.ifi));
vbl = Screen('Flip', dp.wPtr);

for frameIdx = 1:numFrames %#ok<NASGU>
    Screen('FillRect', dp.wPtr, dp.bkColor);
    vbl = Screen('Flip', dp.wPtr, vbl + 0.5 * dp.ifi);
    if shouldAbort(kb)
        abort = true;
        break;
    end
end
end

function abort = presentBreakScreen(dp, kb, waitSeconds)
if nargin < 3 || isempty(waitSeconds)
    waitSeconds = 20;
end

abort = false;

if kb.useKbQueueCheck
    KbQueueFlush;
else
    KbReleaseWait;
end

if isfield(dp, 'textFont') && ~isempty(dp.textFont)
    Screen('TextFont', dp.wPtr, dp.textFont);
end
if isfield(dp, 'textSize') && ~isempty(dp.textSize)
    Screen('TextSize', dp.wPtr, dp.textSize);
end

messageLine = '쉬는시간 입니다. 스페이스바를 누르면 이어서 시작합니다.';
messageLine = double(messageLine);

allowResponseTime = GetSecs + waitSeconds;
firstFrame = true;
vbl = 0;

while true
    Screen('FillRect', dp.wPtr, dp.bkColor);
    bounds = Screen('TextBounds', dp.wPtr, messageLine);
    textX = dp.cx - RectWidth(bounds) / 2;
    textY = dp.cy - RectHeight(bounds) / 2;
    Screen('DrawText', dp.wPtr, messageLine, textX, textY, dp.textColor);

    if firstFrame
        vbl = Screen('Flip', dp.wPtr);
        firstFrame = false;
    else
        vbl = Screen('Flip', dp.wPtr, vbl + 0.5 * dp.ifi);
    end

    if kb.useKbQueueCheck
        [pressed, firstPress] = KbQueueCheck;
        if pressed
            if any(firstPress(kb.escKey))
                abort = true;
                break;
            elseif GetSecs >= allowResponseTime && any(firstPress(kb.spaceKey))
                break;
            end
        end
    else
        [keyIsDown, ~, keyCode] = KbCheck(-1);
        if keyIsDown
            if keyCode(kb.escKey)
                abort = true;
                break;
            elseif GetSecs >= allowResponseTime && keyCode(kb.spaceKey)
                break;
            end
        end
    end
end

if ~abort
    if kb.useKbQueueCheck
        KbQueueFlush;
    else
        KbReleaseWait;
    end
end

if ~firstFrame
    Screen('FillRect', dp.wPtr, dp.bkColor);
    Screen('Flip', dp.wPtr, vbl + 0.5 * dp.ifi);
end
end

function response = collectResponse(dp, kb, durationMs)
response = struct('wasAborted', false, 'didRespond', false, 'keyName', '', 'keyCode', [], 'rt', NaN);

instructionLines = dp.responseInstructions;
if isempty(instructionLines)
    instructionLines = {
        '1번 키: T1 평균이 더 큽니다';
        '2번 키: T2 평균이 더 큽니다'
    };
end

for lineIdx = 1:numel(instructionLines)
    currentLine = instructionLines{lineIdx};
    if isnumeric(currentLine)
        instructionLines{lineIdx} = currentLine;
    else
        if isstring(currentLine)
            currentLine = char(currentLine);
        elseif iscell(currentLine)
            currentLine = char(currentLine{:});
        end
        if ischar(currentLine)
            instructionLines{lineIdx} = double(currentLine);
        else
            instructionLines{lineIdx} = double(string(currentLine));
        end
    end
end

startTime = GetSecs;
deadline  = startTime + durationMs / 1000;

if kb.useKbQueueCheck
    KbQueueFlush;
else
    KbReleaseWait;
end

if isfield(dp, 'textFont') && ~isempty(dp.textFont)
    Screen('TextFont', dp.wPtr, dp.textFont);
end
if isfield(dp, 'textSize') && ~isempty(dp.textSize)
    Screen('TextSize', dp.wPtr, dp.textSize);
end

textSize = Screen('TextSize', dp.wPtr);
if textSize <= 0
    textSize = 24;
    Screen('TextSize', dp.wPtr, textSize);
end

if isfield(dp, 'textLineSpacingMultiplier') && ~isempty(dp.textLineSpacingMultiplier)
    lineSpacingMultiplier = dp.textLineSpacingMultiplier;
else
    lineSpacingMultiplier = 1.2;
end

lineSpacingPx = max(1, round(lineSpacingMultiplier * textSize));
baseY = dp.cy - 1.5 * dp.ppd;

firstFrame = true;
vbl = 0;

while GetSecs < deadline && ~response.wasAborted && ~response.didRespond
    Screen('FillRect', dp.wPtr, dp.bkColor);

    for lineIdx = 1:numel(instructionLines)
        lineY = baseY + (lineIdx - 1) * lineSpacingPx;
        bounds = Screen('TextBounds', dp.wPtr, instructionLines{lineIdx});
        textX = dp.cx - RectWidth(bounds) / 2;
        Screen('DrawText', dp.wPtr, instructionLines{lineIdx}, textX, lineY, dp.textColor);
    end

    if firstFrame
        vbl = Screen('Flip', dp.wPtr);
        firstFrame = false;
    else
        vbl = Screen('Flip', dp.wPtr, vbl + 0.5 * dp.ifi);
    end

    if kb.useKbQueueCheck
        [pressed, firstPress] = KbQueueCheck;
        if pressed
            if any(firstPress(kb.escKey))
                response.wasAborted = true;
                break;
            else
                % 1-key
                if any(firstPress(kb.oneKey) > 0)
                    response.didRespond = true;
                    response.keyCode    = kb.oneKey;
                    response.keyName    = '1';
                    rtTimes = firstPress(kb.oneKey(firstPress(kb.oneKey) > 0));
                    response.rt = min(rtTimes) - startTime;
                % 2-key
                elseif any(firstPress(kb.twoKey) > 0)
                    response.didRespond = true;
                    response.keyCode    = kb.twoKey;
                    response.keyName    = '2';
                    rtTimes = firstPress(kb.twoKey(firstPress(kb.twoKey) > 0));
                    response.rt = min(rtTimes) - startTime;
                end
            end
        end
    else
        [keyIsDown, keyTime, keyCode] = KbCheck(-1);
        if keyIsDown
            if keyCode(kb.escKey)
                response.wasAborted = true;
                break;
            elseif any(keyCode(kb.oneKey))
                response.didRespond = true;
                response.keyCode    = kb.oneKey;
                response.keyName    = '1';
                response.rt         = keyTime - startTime;
            elseif any(keyCode(kb.twoKey))
                response.didRespond = true;
                response.keyCode    = kb.twoKey;
                response.keyName    = '2';
                response.rt         = keyTime - startTime;
            end
        end
    end
end

if ~firstFrame
    Screen('FillRect', dp.wPtr, dp.bkColor);
    Screen('Flip', dp.wPtr, vbl + 0.5 * dp.ifi);
end
end

function [xPosDeg, yPosDeg] = updatePositions(xPosDeg, yPosDeg, motionParams, stim, dp)
stepDeg = motionParams.speedDegPerSec / dp.frameRate;

direction = lower(motionParams.direction);
switch direction
    case 'up'
        yPosDeg = yPosDeg - stepDeg;
        yPosDeg = stim.topEdgesDeg + mod(yPosDeg - stim.topEdgesDeg, stim.verticalSpanDeg);
    case 'down'
        yPosDeg = yPosDeg + stepDeg;
        yPosDeg = stim.topEdgesDeg + mod(yPosDeg - stim.topEdgesDeg, stim.verticalSpanDeg);
    case 'left'
        xPosDeg = xPosDeg - stepDeg;
        xPosDeg = stim.leftEdgesDeg + mod(xPosDeg - stim.leftEdgesDeg, stim.horizontalSpanDeg);
    case 'right'
        xPosDeg = xPosDeg + stepDeg;
        xPosDeg = stim.leftEdgesDeg + mod(xPosDeg - stim.leftEdgesDeg, stim.horizontalSpanDeg);
    otherwise
        error('Unknown motion direction: %s', motionParams.direction);
end
end

function trialResult = populateTrialResult(trialIndex, trialSpec, stimInfo, response, params)
trialResult = initTrialResultStruct();
trialResult.trialIndex = trialIndex;
trialResult.comboLabel = trialSpec.comboLabel;
trialResult.comboMotion = trialSpec.motionCombo;
trialResult.methodPair = trialSpec.methodPair;
trialResult.methodPairLabel = trialSpec.methodPair;
trialResult.freqRatio = trialSpec.freqRatio;
trialResult.equalFreq = trialSpec.equalFreq;
trialResult.diffLevelTarget = trialSpec.diffLevel;
trialResult.signS = trialSpec.signS;
trialResult.largerSide = stimInfo.largerSide;
trialResult.isMoving = trialSpec.isMoving;
trialResult.t1MeanPS = stimInfo.t1MeanPS;
trialResult.t2MeanPS = stimInfo.t2MeanPS;
trialResult.psDiffObserved = stimInfo.psDiffObserved;
trialResult.diffLevelObservedPS = stimInfo.diffLevelObservedPS;
trialResult.t1MeanDeg = stimInfo.t1MeanDeg;
trialResult.t2MeanDeg = stimInfo.t2MeanDeg;
trialResult.expoUsed = params.perceptualExponent;
trialResult.jitterPS_applied = stimInfo.jitterPS;
trialResult.targetPSMeans = stimInfo.targetPSMeans;
trialResult.psTargetsAfterJitter = stimInfo.psTargetsAfterJitter;

if response.didRespond
    trialResult.didRespond    = true;
    trialResult.responseKey   = response.keyName;  % '1' or '2'
    trialResult.responseRtMs  = response.rt * 1000;
    if strcmpi(response.keyName, '1')
        trialResult.responseChoice = 'T1';
    elseif strcmpi(response.keyName, '2')
        trialResult.responseChoice = 'T2';
    else
        trialResult.responseChoice = '';
    end
else
    trialResult.didRespond    = false;
    trialResult.responseKey   = '';
    trialResult.responseChoice = '';
    trialResult.responseRtMs  = NaN;
end

if trialResult.didRespond
    if strcmp(trialResult.responseChoice, 'T1')
        trialResult.correct = stimInfo.t1MeanPS > stimInfo.t2MeanPS;
    elseif strcmp(trialResult.responseChoice, 'T2')
        trialResult.correct = stimInfo.t2MeanPS > stimInfo.t1MeanPS;
    else
        trialResult.correct = NaN;
    end
else
    trialResult.correct = NaN;
end
end

function layout = buildCentralGridLayout(gridConfig)
inner.centerDeg = [0, 0];
inner.widthDeg  = gridConfig.windowWidthDeg;
inner.heightDeg = gridConfig.windowHeightDeg;
inner           = updateApertureEdges(inner);

outer = inner;
outer.widthDeg  = inner.widthDeg  + 2 * gridConfig.outerPaddingDeg;
outer.heightDeg = inner.heightDeg + 2 * gridConfig.outerPaddingDeg;
outer           = updateApertureEdges(outer);

cellWidthDeg  = inner.widthDeg  / gridConfig.cols;
cellHeightDeg = inner.heightDeg / gridConfig.rows;

colCenters = inner.leftDeg + (0.5:1:gridConfig.cols-0.5) * cellWidthDeg;
rowCenters = inner.topDeg  + (0.5:1:gridConfig.rows-0.5) * cellHeightDeg;
[gridX, gridY] = meshgrid(colCenters, rowCenters);

layout.innerAperture = inner;
layout.outerAperture = outer;
layout.gridX         = gridX(:);
layout.gridY         = gridY(:);
layout.cellWidthDeg  = cellWidthDeg;
layout.cellHeightDeg = cellHeightDeg;
end

function aperture = updateApertureEdges(aperture)
halfWidth  = aperture.widthDeg / 2;
halfHeight = aperture.heightDeg / 2;

aperture.leftDeg   = aperture.centerDeg(1) - halfWidth;
aperture.rightDeg  = aperture.centerDeg(1) + halfWidth;
aperture.topDeg    = aperture.centerDeg(2) - halfHeight;
aperture.bottomDeg = aperture.centerDeg(2) + halfHeight;
end

function val = ternary(cond, trueVal, falseVal)
if cond
    val = trueVal;
else
    val = falseVal;
end
end

function psVals = toPS(xDeg, expo)
psVals = xDeg .^ expo;
end

function xDeg = fromPS(psVals, expo)
xDeg = psVals .^ (1/expo);
end

function m = meanPS(xDeg, expo)
m = mean(xDeg .^ expo);
end

function abort = shouldAbort(kb)
[keyIsDown, ~, keyCode] = KbCheck(-1);
abort = keyIsDown && any(keyCode(kb.escKey));
end
