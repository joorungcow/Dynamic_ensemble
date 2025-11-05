AssertOpenGL;
KbName('UnifyKeyNames');

kb = struct();
kb.useKbQueueCheck = 1;
kb = init_keyboard(kb);

%% Display settings
dp.screenNum = max(Screen('Screens'));

dp.dist   = 55;                 % 관찰자와 화면 사이 거리(cm)
dp.width  = 60;                 % 사용 중인 디스플레이의 가로 폭(cm)
dp.bkColor   = [0.5 0.5 0.5];
dp.textColor = [1 1 1];
dp.textFont  = 'Malgun Gothic';
dp.textSize  = 20;
dp.textLineSpacingMultiplier = 3;
dp.responseInstructions = {
    double('왼쪽 방향키: T1 평균이 더 큽니다');
    double('오른쪽 방향키: T2 평균이 더 큽니다')
};
dp.skipChecks = 0;

saveFileName = sprintf('results_test9_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));

try
    dp = OpenWindow(dp);
    HideCursor;
    ListenChar(2);

    %% Dot configuration
    dotParams.smallSizeDeg = 0.7;              % 작은 점 집단의 평균 지름(시야각)
    dotParams.largeSizeDeg = 1.3;              % 큰 점 집단의 평균 지름(시야각)
    dotParams.targetMeanDeg     = 1.0;         % 목표 평균 지름(시야각)
    dotParams.meanJitterDeg     = 0.05;        % 목표 평균에 추가되는 랜덤 지터 범위(시야각)
    dotParams.minSizeDeg        = 0.4;         % 생성 가능한 점의 최소 지름(시야각)
    dotParams.maxSizeDeg        = 1.8;         % 생성 가능한 점의 최대 지름(시야각)
    dotParams.gToleranceDeg     = 0.001;       % 점 생성 시 허용 오차(시야각)
    dotParams.jitterStdRatio    = 0.15;        % 지터 표준편차에 대한 비율
    dotParams.perceptualExponent = 0.76;       % 지각적 크기 변환에 사용하는 지수 값
    dotParams.meanDiffLevels    = [0.06 0.12 0.18 0.24 0.30 0.36]; % 두 자극 간 평균 차이 수준(시야각) 0.06 0.12 0.18 0.24 0.30 0.36
    dotParams.safetyMarginDeg   = 0.05;        % 점이 경계에 겹치지 않도록 확보하는 안전 여유(시야각)

    ratioAssignments = {
        struct('label','S6L2_vs_S2L6','t1Counts',[6 2],'t2Counts',[2 6]); % 작은/큰 점 개수 조합(6:2 vs 2:6)
        struct('label','S5L3_vs_S3L5','t1Counts',[5 3],'t2Counts',[3 5]); % 작은/큰 점 개수 조합(5:3 vs 3:5)
        struct('label','S3L5_vs_S5L3','t1Counts',[3 5],'t2Counts',[5 3]); % 작은/큰 점 개수 조합(3:5 vs 5:3)
        struct('label','S2L6_vs_S6L2','t1Counts',[2 6],'t2Counts',[6 2])  % 작은/큰 점 개수 조합(2:6 vs 6:2)
    };

    numDots = sum(ratioAssignments{1}.t1Counts); % 한 자극에 표시되는 전체 점 개수

    %% Motion configuration
    motionParams.direction        = 'up';    % 움직이는 자극의 이동 방향
    motionParams.speedDegPerSec   = 1.5;     % 움직이는 자극의 속도(시야각/초)

    %% Timing configuration (ms)
    timingParams.fixationMs     = 300;  % 응시점 제시 시간(ms)
    timingParams.stimDurationMs = 500;  % 각 자극 제시 시간(ms)
    timingParams.isiDurationMs  = 200; % 자극 사이 공백 지속 시간(ms)
    timingParams.postTrialMs    = 1000; % 반응 대기 및 안내 문구 표시 시간(ms)

    comboRepeats = 1;                   % 자극 조합 반복 횟수
    stimCombos = {'MM','SM','MS','SS'}; % 움직임/정지 조합(M: 움직임, S: 정지)

    gridConfig.rows = 6;
    gridConfig.cols = 6;
    gridConfig.windowWidthDeg = 6 * 2;       % 자극 제시 창의 가로 크기(시야각)
    gridConfig.windowHeightDeg = 6 * 2;      % 자극 제시 창의 세로 크기(시야각)
    gridConfig.maxJitterDeg = 0.12;             % 셀 중심 기준 최대 지터(시야각)
    gridConfig.maxAttemptsPerCell = 50;
    gridConfig.safetyMarginDeg = dotParams.safetyMarginDeg;
    gridConfig.outerPaddingDeg = dotParams.maxSizeDeg / 2 + dotParams.safetyMarginDeg;

    gridLayout = buildCentralGridLayout(gridConfig);

    epsDeg = 1 / dp.ppd; % 픽셀당 시야각 환산값
    fprintf('Estimated pixel pitch: %.4f°/pixel.\n', epsDeg);
    if dotParams.gToleranceDeg < 0.5 * epsDeg
        fprintf(['[경고] gTolerance가 디스플레이 상도(%.4f°)보다 매우 작습니다. ' ...
                 '계산상 문제는 없으나 물리적으로 구분이 어려울 수 있습니다.\n'], epsDeg);
    end

    Screen('FillRect', dp.wPtr, dp.bkColor);
    Screen('Flip', dp.wPtr);

    condList = buildTrialConditions(stimCombos, ratioAssignments, dotParams.meanDiffLevels, comboRepeats);
    order = randperm(numel(condList));

    result.trials(numel(order), 1) = struct( ...
        'condIdx', [], ...
        'comboLabel', '', ...
        'ratioIdx', [], ...
        'diffLevel', [], ...
        't1TargetMeanDeg', [], ...
        't2TargetMeanDeg', [], ...
        'responseKey', '', ...
        'responseChoice', '', ...
        'responseRtMs', NaN, ...
        'didRespond', false, ...
        'correct', NaN);

    abortExperiment = false;
    trialCounter = 1;

    for idx = order
        cond = condList(idx);

        if abortExperiment
            break;
        end

        comboLabel = cond.comboLabel;
        t1IsMoving = comboLabel(1) == 'M';
        t2IsMoving = comboLabel(2) == 'M';

        [t1TargetMeanDeg, t2TargetMeanDeg, t1Counts, t2Counts] = ...
            computeTargetMeans(dotParams, ratioAssignments{cond.ratioIdx}, cond.diffLevel);

        t1TargetMeanDeg = jitterMean(t1TargetMeanDeg, dotParams);
        t2TargetMeanDeg = jitterMean(t2TargetMeanDeg, dotParams);

        % Fixation screen before T1 onset
        if timingParams.fixationMs > 0
            numFixFrames = max(1, round((timingParams.fixationMs/1000) / dp.ifi));
            firstFixationFrame = true;
            for frameIdx = 1:numFixFrames
                Screen('FillRect', dp.wPtr, dp.bkColor);
                make_fixation(dp);
                if firstFixationFrame
                    vbl = Screen('Flip', dp.wPtr);
                    firstFixationFrame = false;
                else
                    vbl = Screen('Flip', dp.wPtr, vbl + 0.5 * dp.ifi);
                end
                if shouldAbort(kb)
                    abortExperiment = true;
                    break;
                end
            end
            if abortExperiment
                break;
            end
        end

        abortExperiment = showStimulus(numDots, t1TargetMeanDeg, dotParams, gridLayout, gridConfig, t1Counts, t1IsMoving, motionParams, timingParams.stimDurationMs, dp, kb);
        if abortExperiment
            break;
        end

        abortExperiment = presentBlank(dp, timingParams.isiDurationMs, kb);
        if abortExperiment
            break;
        end

        abortExperiment = showStimulus(numDots, t2TargetMeanDeg, dotParams, gridLayout, gridConfig, t2Counts, t2IsMoving, motionParams, timingParams.stimDurationMs, dp, kb);
        if abortExperiment
            break;
        end

        response = collectResponse(dp, kb, timingParams.postTrialMs);
        if response.wasAborted
            abortExperiment = true;
            break;
        end

        trialResult = struct();
        trialResult.condIdx = idx;
        trialResult.comboLabel = cond.comboLabel;
        trialResult.ratioIdx = cond.ratioIdx;
        trialResult.diffLevel = cond.diffLevel;
        trialResult.t1TargetMeanDeg = t1TargetMeanDeg;
        trialResult.t2TargetMeanDeg = t2TargetMeanDeg;
        trialResult.responseKey = response.keyName;
        trialResult.didRespond = response.didRespond;
        trialResult.responseRtMs = response.rt * 1000;

        if response.didRespond
            if strcmpi(response.keyName, 'LeftArrow')
                trialResult.responseChoice = 'T1';
            elseif strcmpi(response.keyName, 'RightArrow')
                trialResult.responseChoice = 'T2';
            else
                trialResult.responseChoice = '';
            end
        else
            trialResult.responseChoice = '';
        end

        if response.didRespond
            if strcmp(trialResult.responseChoice, 'T1')
                trialResult.correct = t1TargetMeanDeg > t2TargetMeanDeg;
            elseif strcmp(trialResult.responseChoice, 'T2')
                trialResult.correct = t2TargetMeanDeg > t1TargetMeanDeg;
            else
                trialResult.correct = NaN;
            end
        else
            trialResult.correct = NaN;
        end

        result.trials(trialCounter) = trialResult;
        trialCounter = trialCounter + 1;
    end

    result.trials = result.trials(1:trialCounter-1);

    presentBlank(dp, 0, kb);

    Screen('CloseAll');
    ListenChar(0);
    ShowCursor;
    if kb.useKbQueueCheck
        KbQueueRelease;
    end
    RestrictKeysForKbCheck([]);
    save(saveFileName, 'result');
catch ME
    Screen('CloseAll');
    ListenChar(0);
    ShowCursor;
    if exist('kb','var') && isfield(kb,'useKbQueueCheck') && kb.useKbQueueCheck
        KbQueueRelease;
    end
    RestrictKeysForKbCheck([]);
    rethrow(ME);
end

%% Local functions -------------------------------------------------------
function layout = buildCentralGridLayout(gridConfig)
inner.centerDeg = [0, 0];
inner.widthDeg  = gridConfig.windowWidthDeg;
inner.heightDeg = gridConfig.windowHeightDeg;
inner = updateApertureEdges(inner);

outer = inner;
outer.widthDeg  = inner.widthDeg  + 2 * gridConfig.outerPaddingDeg;
outer.heightDeg = inner.heightDeg + 2 * gridConfig.outerPaddingDeg;
outer = updateApertureEdges(outer);

cellWidthDeg = inner.widthDeg / gridConfig.cols;
cellHeightDeg = inner.heightDeg / gridConfig.rows;

colCenters = inner.leftDeg + (0.5:1:gridConfig.cols-0.5) * cellWidthDeg;
rowCenters = inner.topDeg  + (0.5:1:gridConfig.rows-0.5) * cellHeightDeg;
[gridX, gridY] = meshgrid(colCenters, rowCenters);

layout.innerAperture = inner;
layout.outerAperture = outer;
layout.gridX = gridX(:);
layout.gridY = gridY(:);
layout.cellWidthDeg = cellWidthDeg;
layout.cellHeightDeg = cellHeightDeg;
end

function abort = showStimulus(numDots, targetMeanDeg, params, layout, gridConfig, ratioCounts, isMoving, motionParams, stimDurationMs, dp, kb)
stim = createStimulus(numDots, targetMeanDeg, params, layout, gridConfig, ratioCounts);
abort = presentStimulus(dp, stim, isMoving, motionParams, stimDurationMs, kb);
end

function stim = createStimulus(numDots, targetMeanDeg, params, layout, gridConfig, ratioCounts)
stim.dotSizeDeg = generateDotSizes(numDots, targetMeanDeg, params, ratioCounts);
[stim.xPosDeg, stim.yPosDeg] = placeDotsOnGrid(stim.dotSizeDeg, layout, gridConfig);

stim.outerAperture = layout.outerAperture;
stim.innerAperture = layout.innerAperture;
stim.halfSizeDeg = stim.dotSizeDeg / 2;

stim.topEdgesDeg    = stim.outerAperture.topDeg    + stim.halfSizeDeg;
stim.bottomEdgesDeg = stim.outerAperture.bottomDeg - stim.halfSizeDeg;
stim.leftEdgesDeg   = stim.outerAperture.leftDeg   + stim.halfSizeDeg;
stim.rightEdgesDeg  = stim.outerAperture.rightDeg  - stim.halfSizeDeg;

stim.verticalSpanDeg   = stim.bottomEdgesDeg - stim.topEdgesDeg;
stim.horizontalSpanDeg = stim.rightEdgesDeg - stim.leftEdgesDeg;

stim.innerTopEdgesDeg    = stim.innerAperture.topDeg    + stim.halfSizeDeg;
stim.innerBottomEdgesDeg = stim.innerAperture.bottomDeg - stim.halfSizeDeg;
stim.innerLeftEdgesDeg   = stim.innerAperture.leftDeg   + stim.halfSizeDeg;
stim.innerRightEdgesDeg  = stim.innerAperture.rightDeg  - stim.halfSizeDeg;

if any(stim.innerBottomEdgesDeg <= stim.innerTopEdgesDeg) || any(stim.innerRightEdgesDeg <= stim.innerLeftEdgesDeg)
    error('Inner aperture must be larger than dot diameters.');
end

if any(stim.verticalSpanDeg <= 0) || any(stim.horizontalSpanDeg <= 0)
    error('Outer aperture size must exceed the diameter of every dot.');
end
end

function dotSizesDeg = generateDotSizes(numDots, targetMeanDeg, params, ratioCounts)
maxAttempts = 1000;
baseSizes = [repmat(params.smallSizeDeg,1,ratioCounts(1)), repmat(params.largeSizeDeg,1,ratioCounts(2))];

if numel(baseSizes) ~= numDots
    error('Ratio counts must sum to %d dots.', numDots);
end

scaleToTarget = targetMeanDeg / mean(baseSizes);
baseSizes = baseSizes * scaleToTarget;

jitterStdDeg = params.jitterStdRatio * targetMeanDeg;

for attempt = 1:maxAttempts %#ok<NASGU>
    noise = jitterStdDeg .* randn(1, numDots);
    noise = noise - mean(noise);

    dotSizesDeg = baseSizes + noise;
    dotSizesDeg = min(max(dotSizesDeg, params.minSizeDeg), params.maxSizeDeg);

    currentMean = mean(dotSizesDeg);
    delta = targetMeanDeg - currentMean;
    dotSizesDeg = dotSizesDeg + delta/numDots;

    if any(dotSizesDeg < params.minSizeDeg) || any(dotSizesDeg > params.maxSizeDeg)
        continue;
    end

    if abs(mean(dotSizesDeg) - targetMeanDeg) <= params.gToleranceDeg
        return;
    end
end

error('Could not generate dot sizes within tolerance. Adjust parameters.');
end

function [xPosDeg, yPosDeg] = placeDotsOnGrid(dotSizesDeg, layout, gridConfig)
numDots = numel(dotSizesDeg);
occupiedCells = false(numel(layout.gridX), 1);

xPosDeg = zeros(1, numDots);
yPosDeg = zeros(1, numDots);

for dotIdx = 1:numDots
    placed = false;
    cellOrder = randperm(numel(occupiedCells));
    currentSize = dotSizesDeg(dotIdx);
    halfSize = currentSize / 2;

    cellJitterX = layout.cellWidthDeg / 2 - (halfSize + gridConfig.safetyMarginDeg);
    cellJitterY = layout.cellHeightDeg / 2 - (halfSize + gridConfig.safetyMarginDeg);

    if cellJitterX <= 0 || cellJitterY <= 0
        error('Cell dimensions are too small for dot %d.', dotIdx);
    end

    for candidate = cellOrder
        if occupiedCells(candidate)
            continue;
        end

        baseX = layout.gridX(candidate);
        baseY = layout.gridY(candidate);

        leftLimit   = baseX - (layout.innerAperture.leftDeg  + halfSize + gridConfig.safetyMarginDeg);
        rightLimit  = (layout.innerAperture.rightDeg  - halfSize - gridConfig.safetyMarginDeg) - baseX;
        topLimit    = baseY - (layout.innerAperture.topDeg    + halfSize + gridConfig.safetyMarginDeg);
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

            if dotIdx > 1
                dx = xPosDeg(1:dotIdx-1) - candidateX;
                dy = yPosDeg(1:dotIdx-1) - candidateY;
                minDist = (dotSizesDeg(1:dotIdx-1) + currentSize)/2 + gridConfig.safetyMarginDeg;

                if any((dx.^2 + dy.^2) < (minDist.^2))
                    continue;
                end
            end

            xPosDeg(dotIdx) = candidateX;
            yPosDeg(dotIdx) = candidateY;
            occupiedCells(candidate) = true;
            placed = true;
            break;
        end

        if placed
            break;
        end
    end

    if ~placed
        error('Failed to place dot %d within 6x6 grid without overlap.', dotIdx);
    end
end
end

function condList = buildTrialConditions(stimCombos, ratioAssignments, diffLevels, repeats)
condIdx = 1;
for rep = 1:repeats
    for ratioIdx = 1:numel(ratioAssignments)
        for diffIdx = 1:numel(diffLevels)
            for comboIdx = 1:numel(stimCombos)
                condList(condIdx).comboLabel = stimCombos{comboIdx}; %#ok<AGROW>
                condList(condIdx).ratioIdx   = ratioIdx;
                condList(condIdx).diffLevel  = diffLevels(diffIdx);
                condIdx = condIdx + 1;
            end
        end
    end
end
end

function [t1TargetMeanDeg, t2TargetMeanDeg, t1Counts, t2Counts] = computeTargetMeans(params, ratioStruct, diffLevel)
smallDeg = params.smallSizeDeg;
largeDeg = params.largeSizeDeg;
expo     = params.perceptualExponent;

t1Counts = ratioStruct.t1Counts;
t2Counts = ratioStruct.t2Counts;

baseT1 = [repmat(smallDeg,1,t1Counts(1)), repmat(largeDeg,1,t1Counts(2))];
baseT2 = [repmat(smallDeg,1,t2Counts(1)), repmat(largeDeg,1,t2Counts(2))];

meanT1deg = mean(baseT1);
meanT2deg = mean(baseT2);

meanT1ps = mean(baseT1.^expo);
meanT2ps = mean(baseT2.^expo);

diffSign = sign(meanT2ps - meanT1ps);
if diffSign == 0
    diffSign = 1;
end

t1TargetPS = meanT1ps;
t2TargetPS = t1TargetPS * (1 + diffSign * diffLevel);

scaleT1 = (t1TargetPS / meanT1ps)^(1/expo);
scaleT2 = (t2TargetPS / meanT2ps)^(1/expo);

t1TargetMeanDeg = meanT1deg * scaleT1;
t2TargetMeanDeg = meanT2deg * scaleT2;

minDeg = params.minSizeDeg;
maxDeg = params.maxSizeDeg;

t1TargetMeanDeg = min(max(t1TargetMeanDeg, minDeg), maxDeg);
t2TargetMeanDeg = min(max(t2TargetMeanDeg, minDeg), maxDeg);
end

function meanDeg = jitterMean(baseMeanDeg, params)
meanDeg = baseMeanDeg + params.meanJitterDeg * (2 * rand - 1);
meanDeg = min(max(meanDeg, params.minSizeDeg), params.maxSizeDeg);
end

function abort = presentStimulus(dp, stim, isMoving, motionParams, stimDurationMs, kb)
abort = false;
numFrames = max(1, round((stimDurationMs/1000) / dp.ifi));

Screen('FillRect', dp.wPtr, dp.bkColor);
vbl = Screen('Flip', dp.wPtr);

for frameIdx = 1:numFrames
    if isMoving && frameIdx > 1
        [stim.xPosDeg, stim.yPosDeg] = updatePositions(stim.xPosDeg, stim.yPosDeg, motionParams, stim, dp);
    end

    Screen('FillRect', dp.wPtr, dp.bkColor);

    xyPix = convertCentersToPixels(stim.xPosDeg, stim.yPosDeg, dp);
    sizePix = stim.dotSizeDeg * dp.ppd;

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

function xyPix = convertCentersToPixels(xDeg, yDeg, dp)
xyPix = [xDeg; yDeg] * dp.ppd;
xyPix(1, :) = xyPix(1, :) + dp.cx;
xyPix(2, :) = xyPix(2, :) + dp.cy;
end

function aperture = updateApertureEdges(aperture)
halfWidth  = aperture.widthDeg / 2;
halfHeight = aperture.heightDeg / 2;

aperture.leftDeg   = aperture.centerDeg(1) - halfWidth;
aperture.rightDeg  = aperture.centerDeg(1) + halfWidth;
aperture.topDeg    = aperture.centerDeg(2) - halfHeight;
aperture.bottomDeg = aperture.centerDeg(2) + halfHeight;
end

function response = collectResponse(dp, kb, durationMs)
response = struct('didRespond', false, ...
                  'keyCode', NaN, ...
                  'keyName', '', ...
                  'rt', NaN, ...
                  'wasAborted', false);

if durationMs <= 0
    return;
end

if isfield(dp, 'responseInstructions') && ~isempty(dp.responseInstructions)
    if iscell(dp.responseInstructions)
        instructionLines = dp.responseInstructions(:);
    else
        instructionLines = cellstr(dp.responseInstructions);
    end
else
    instructionLines = {
        '왼쪽 방향키: T1 평균이 더 큽니다';
        '오른쪽 방향키: T2 평균이 더 큽니다'
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
deadline = startTime + durationMs / 1000;

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
            elseif firstPress(kb.leftKey) > 0
                response.didRespond = true;
                response.keyCode = kb.leftKey;
                response.keyName = 'LeftArrow';
                response.rt = firstPress(kb.leftKey) - startTime;
            elseif firstPress(kb.rightKey) > 0
                response.didRespond = true;
                response.keyCode = kb.rightKey;
                response.keyName = 'RightArrow';
                response.rt = firstPress(kb.rightKey) - startTime;
            end
        end
    else
        [keyIsDown, keyTime, keyCode] = KbCheck(-1);
        if keyIsDown
            if keyCode(kb.escKey)
                response.wasAborted = true;
                break;
            elseif keyCode(kb.leftKey)
                response.didRespond = true;
                response.keyCode = kb.leftKey;
                response.keyName = 'LeftArrow';
                response.rt = keyTime - startTime;
            elseif keyCode(kb.rightKey)
                response.didRespond = true;
                response.keyCode = kb.rightKey;
                response.keyName = 'RightArrow';
                response.rt = keyTime - startTime;
            end
        end
    end
end

if ~firstFrame
    Screen('FillRect', dp.wPtr, dp.bkColor);
    Screen('Flip', dp.wPtr, vbl + 0.5 * dp.ifi);
end

end

function abort = shouldAbort(kb)
[keyIsDown, ~, keyCode] = KbCheck(-1);
abort = keyIsDown && any(keyCode(kb.escKey));
end