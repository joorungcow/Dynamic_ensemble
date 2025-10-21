AssertOpenGL;
KbName('UnifyKeyNames');

kb = struct();
kb.useKbQueueCheck = 0;
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
dp.skipChecks = 1;

saveFileName = sprintf('results_test8_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));

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
    dotParams.meanDiffLevels    = [0.06 0.12 0.18 0.24 0.30 0.36]; % 두 자극 간 평균 차이 수준(시야각)
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
    timingParams.isiDurationMs  = 1000; % 자극 사이 공백 지속 시간(ms)
    timingParams.postTrialMs    = 1000; % 반응 대기 및 안내 문구 표시 시간(ms)

    comboRepeats = 1;                   % 자극 조합 반복 횟수
    stimCombos = {'MM','SM','MS','SS'}; % 움직임/정지 조합(M: 움직임, S: 정지)

    gridConfig.rows = 6;
    gridConfig.cols = 6;
    gridConfig.windowWidthDeg = 6 * 2.65;       % 자극 제시 창의 가로 크기(시야각)
    gridConfig.windowHeightDeg = 6 * 2.65;      % 자극 제시 창의 세로 크기(시야각)
    gridConfig.maxJitterDeg = 0.12;             % 셀 중심 기준 최대 지터(시야각)
    gridConfig.maxAttemptsPerCell = 50;
    gridConfig.safetyMarginDeg = dotParams.safetyMarginDeg;
    gridConfig.outerPaddingDeg = dotParams.maxSizeDeg / 2 + dotParams.safetyMarginDeg;

    gridLayout.central = buildCentralGridLayout(gridConfig);

    epsDeg = 1 / dp.ppd; % 픽셀당 시야각 환산값
    fprintf('Estimated pixel pitch: %.4f°/pixel.\n', epsDeg);
    if dotParams.gToleranceDeg < 0.5 * epsDeg
        fprintf(['[경고] gTolerance가 디스플레이 해상도(%.4f°)보다 매우 작습니다. ' ...
                 '계산상 문제는 없으나 물리적으로 구분이 어려울 수 있습니다.\n'], epsDeg);
    end

    Screen('FillRect', dp.wPtr, dp.bkColor);
    Screen('Flip', dp.wPtr);

    condList = design.buildTrialConditions(stimCombos, ratioAssignments, dotParams.meanDiffLevels, comboRepeats);
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
            design.computeTargetMeans(dotParams, ratioAssignments{cond.ratioIdx}, cond.diffLevel);

        t1TargetMeanDeg = design.jitterMean(t1TargetMeanDeg, dotParams);
        t2TargetMeanDeg = design.jitterMean(t2TargetMeanDeg, dotParams);

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
                if input.shouldAbort(kb)
                    abortExperiment = true;
                    break;
                end
            end
            if abortExperiment
                break;
            end
        end

        stim1 = stim.createStimulusStruct(numDots, t1TargetMeanDeg, dotParams, gridLayout.central.outerAperture, gridLayout.central.innerAperture, t1Counts);
        [stim1.xPosDeg, stim1.yPosDeg] = placeDotsOnGrid(stim1.dotSizeDeg, gridLayout.central, gridConfig);
        abortExperiment = render.presentStimulus(dp, stim1, t1IsMoving, motionParams, timingParams.stimDurationMs, kb);
        if abortExperiment
            break;
        end

        abortExperiment = render.presentBlank(dp, timingParams.isiDurationMs, kb);
        if abortExperiment
            break;
        end

        stim2 = stim.createStimulusStruct(numDots, t2TargetMeanDeg, dotParams, gridLayout.central.outerAperture, gridLayout.central.innerAperture, t2Counts);
        [stim2.xPosDeg, stim2.yPosDeg] = placeDotsOnGrid(stim2.dotSizeDeg, gridLayout.central, gridConfig);
        abortExperiment = render.presentStimulus(dp, stim2, t2IsMoving, motionParams, timingParams.stimDurationMs, kb);
        if abortExperiment
            break;
        end

        response = input.collectResponse(dp, kb, timingParams.postTrialMs);
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

    render.presentBlank(dp, 0, kb);

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
inner = geom.updateApertureEdges(inner);

outer = inner;
outer.widthDeg  = inner.widthDeg  + 2 * gridConfig.outerPaddingDeg;
outer.heightDeg = inner.heightDeg + 2 * gridConfig.outerPaddingDeg;
outer = geom.updateApertureEdges(outer);

cellWidthDeg = inner.widthDeg / gridConfig.cols;
cellHeightDeg = inner.heightDeg / gridConfig.rows;

colCenters = inner.leftDeg + (0.5:1:gridConfig.cols-0.5) * cellWidthDeg;
rowCenters = inner.topDeg  + (0.5:1:gridConfig.rows-0.5) * cellHeightDeg;
[gridX, gridY] = meshgrid(colCenters, rowCenters);

layout.innerAperture = inner;
layout.outerAperture = outer;
layout.gridX = gridX(:);
layout.gridY = gridY(:);
end

function [xPosDeg, yPosDeg] = placeDotsOnGrid(dotSizesDeg, layout, gridConfig)
numDots = numel(dotSizesDeg);
occupiedCells = false(numel(layout.gridX), 1);

cellWidthDeg = layout.innerAperture.widthDeg / gridConfig.cols;
cellHeightDeg = layout.innerAperture.heightDeg / gridConfig.rows;

xPosDeg = zeros(1, numDots);
yPosDeg = zeros(1, numDots);

for dotIdx = 1:numDots
    placed = false;
    cellOrder = randperm(numel(occupiedCells));
    currentSize = dotSizesDeg(dotIdx);
    halfSize = currentSize / 2;

    cellJitterX = cellWidthDeg / 2 - (halfSize + gridConfig.safetyMarginDeg);
    cellJitterY = cellHeightDeg / 2 - (halfSize + gridConfig.safetyMarginDeg);

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