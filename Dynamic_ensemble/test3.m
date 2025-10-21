Screen('Preference','SkipSyncTests',1);

PsychDefaultSetup(2);         % 0~1 색 범위, 키맵 통일
AssertOpenGL;                 % PTB OpenGL 경로/드라이버 체크
KbName('UnifyKeyNames');

%% Display settings (params)
dp.screenNum = max(Screen('Screens'));

dp.dist   = 55;   % 관찰 거리 (cm)
dp.width  = 60;   % 화면 가로 물리폭 (cm)
dp.bkColor   = 0.5;            % 배경색 (회색)
dp.textColor = [0 0 0];        % 텍스트 색 (검정)

% 해상도/주사율/비율/ppd 미리 계산(창 열기 없이 가능)
d = Screen('Resolution', dp.screenNum);
dp.resolution  = [d.width, d.height];        % [W H] px
dp.frameRate   = d.hz;                        % Hz

% 픽셀-각도 변환(ppd, pixels per degree)
dp.ppd = dp.resolution(1) / ((2*atan(dp.width/(2*dp.dist)))*180/pi);

rect = [];
%rect = [0 0 800 600];


try
    [dp.wPtr, dp.wRect] = PsychImaging('OpenWindow', dp.screenNum, dp.bkColor, rect, [], [], 0);

    % --- after openwindow : real params update ---
    Screen('BlendFunction', dp.wPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    dp.ifi        = Screen('GetFlipInterval', dp.wPtr);
    dp.frameRate  = round(1/dp.ifi);        % 실제 주사율
    dp.resolution = dp.wRect([3 4]);        % 실제 창 크기
    [dp.cx, dp.cy]= RectCenter(dp.wRect);

    %% dot configuration (degrees)
    dotParams.smallSizeDeg = 0.7;   % 작은 원 지름 (°)
    dotParams.largeSizeDeg = 1.3;   % 큰 원 지름 (°)
    dotParams.targetMeanDeg     = 1.0;   % 기본 목표 평균 (°)
    dotParams.meanJitterDeg     = 0.05;  % trial 간 평균 난수 범위 (°)
    dotParams.minSizeDeg        = 0.4;   % 최소 크기 제한 (°)
    dotParams.maxSizeDeg        = 1.8;   % 최대 크기 제한 (°)
    dotParams.gToleranceDeg     = 0.001;  % 허용 오차 (°)
    dotParams.jitterStdRatio    = 0.15;  % 개별 dot size 산출 시 표준편차 비율
    dotParams.perceptualExponent = 0.76; % 지각 척도 지수
    dotParams.meanDiffLevels    = [0.06 0.12 0.18 0.24 0.30 0.36]; % 지각 척도 상 평균 차이 단계
    dotParams.safetyMarginDeg   = 0.05;  % 비겹침 검사 시 여유 (°)

    % small:large 비율 조합 (합계 8개 dot)
    ratioAssignments = {
        struct('label','S6L2_vs_S2L6','t1Counts',[6 2],'t2Counts',[2 6]);
        struct('label','S5L3_vs_S3L5','t1Counts',[5 3],'t2Counts',[3 5]);
        struct('label','S3L5_vs_S5L3','t1Counts',[3 5],'t2Counts',[5 3]);
        struct('label','S2L6_vs_S6L2','t1Counts',[2 6],'t2Counts',[6 2])
    };

    numDots = sum(ratioAssignments{1}.t1Counts);

    % Motion configuration (deg/sec)
    motionParams.direction        = 'up';    % 'up', 'down', 'left', 'right'
    motionParams.speedDegPerSec   = 3 ;     % 이동 속도 (°/sec)

    % Timing configuration (ms)
    timingParams.stimDurationMs    = 520;   % T1/T2 자극 제시 시간
    timingParams.isiDurationMs     = 1000;   % T1-T2 사이 공백
    timingParams.postTrialMs       = 1000;   % T2 이후 공백

    % 자극 조합 반복 횟수 (모든 조합 × n)
    comboRepeats = 1;

    % 자극 이동 여부 조합 (T1/T2)
    stimCombos = {'MM','SM','MS','SS'};

    % aperture(사각형) 파라미터 (°)
    aperture.centerDeg = [0 0];          % 화면 중앙 기준 (°)
    aperture.widthDeg  = 12;             % 가로 폭 (°)
    aperture.heightDeg = 9;              % 세로 높이 (°)
    aperture = updateApertureEdges(aperture);

    % 물리적 디스플레이 해상도 경고
    epsDeg = 1 / dp.ppd;
    fprintf('Estimated pixel pitch: %.4f°/pixel.\n', epsDeg);
    if dotParams.gToleranceDeg < 0.5 * epsDeg
        fprintf(['[경고] gTolerance가 디스플레이 해상도(%.4f°)보다 매우 작습니다. ' ...
                 '계산상 문제는 없으나 물리적으로 구분이 어려울 수 있습니다.\n'], epsDeg);
    end

    % 초기 화면을 배경으로 클리어
    Screen('FillRect', dp.wPtr, dp.bkColor);
    Screen('Flip', dp.wPtr);

    % trial 조합 생성 및 섞음
    condList = buildTrialConditions(stimCombos, ratioAssignments, dotParams.meanDiffLevels, comboRepeats);
    order = randperm(numel(condList));

    abortExperiment = false;

    for idx = order
        cond = condList(idx);

        if abortExperiment
            break;
        end

        comboLabel = cond.comboLabel;
        t1IsMoving = comboLabel(1) == 'M';
        t2IsMoving = comboLabel(2) == 'M';

        % 목표 평균 계산 (°)
        [t1TargetMeanDeg, t2TargetMeanDeg, t1Counts, t2Counts] = ...
            computeTargetMeans(dotParams, ratioAssignments{cond.ratioIdx}, cond.diffLevel);

        % trial별 평균에 작은 난수 가중치 추가 (°)
        t1TargetMeanDeg = jitterMean(t1TargetMeanDeg, dotParams);
        t2TargetMeanDeg = jitterMean(t2TargetMeanDeg, dotParams);

        stim1 = createStimulusStruct(numDots, t1TargetMeanDeg, dotParams, aperture, t1Counts);
        abortExperiment = presentStimulus(dp, stim1, t1IsMoving, motionParams, timingParams.stimDurationMs);
        if abortExperiment
            break;
        end

        abortExperiment = presentBlank(dp, timingParams.isiDurationMs);
        if abortExperiment
            break;
        end

        stim2 = createStimulusStruct(numDots, t2TargetMeanDeg, dotParams, aperture, t2Counts);
        abortExperiment = presentStimulus(dp, stim2, t2IsMoving, motionParams, timingParams.stimDurationMs);
        if abortExperiment
            break;
        end

        abortExperiment = presentBlank(dp, timingParams.postTrialMs);
    end

    presentBlank(dp, timingParams.postTrialMs);

    Screen('CloseAll');
catch ME
    Screen('CloseAll');
    rethrow(ME);
end

%% --- Local help functions ---
% ApertureEdges 안에서만 stimuli 제시
function aperture = updateApertureEdges(aperture)
halfWidth  = aperture.widthDeg / 2;
halfHeight = aperture.heightDeg / 2;
aperture.leftDeg   = aperture.centerDeg(1) - halfWidth;
aperture.rightDeg  = aperture.centerDeg(1) + halfWidth;
aperture.topDeg    = aperture.centerDeg(2) - halfHeight;
aperture.bottomDeg = aperture.centerDeg(2) + halfHeight;
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

function meanDeg = jitterMean(baseMeanDeg, params)
meanDeg = baseMeanDeg + params.meanJitterDeg * (2*rand - 1);
meanDeg = min(max(meanDeg, params.minSizeDeg), params.maxSizeDeg);
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

% Target 평균이 최소/최대 범위를 벗어나지 않도록 scale 계산
scaleT1 = (t1TargetPS / meanT1ps)^(1/expo);
scaleT2 = (t2TargetPS / meanT2ps)^(1/expo);

t1TargetMeanDeg = meanT1deg * scaleT1;
t2TargetMeanDeg = meanT2deg * scaleT2;

% 범위 체크 (°)
minDeg = params.minSizeDeg;
maxDeg = params.maxSizeDeg;

if t1TargetMeanDeg < minDeg
    t1TargetMeanDeg = minDeg;
elseif t1TargetMeanDeg > maxDeg
    t1TargetMeanDeg = maxDeg;
end

if t2TargetMeanDeg < minDeg
    t2TargetMeanDeg = minDeg;
elseif t2TargetMeanDeg > maxDeg
    t2TargetMeanDeg = maxDeg;
end
end

function stim = createStimulusStruct(numDots, targetMeanDeg, params, aperture, ratioCounts)
dotSizeDeg = generateDotSizes(numDots, targetMeanDeg, params, ratioCounts);
[xPosDeg, yPosDeg] = initializeNonOverlappingPositions(numDots, dotSizeDeg, aperture, params.safetyMarginDeg);

stim.dotSizeDeg = dotSizeDeg;
stim.xPosDeg = xPosDeg;
stim.yPosDeg = yPosDeg;
stim.aperture = aperture;
stim.halfSizeDeg = dotSizeDeg / 2;
stim.topEdgesDeg    = aperture.topDeg    + stim.halfSizeDeg;
stim.bottomEdgesDeg = aperture.bottomDeg - stim.halfSizeDeg;
stim.leftEdgesDeg   = aperture.leftDeg   + stim.halfSizeDeg;
stim.rightEdgesDeg  = aperture.rightDeg  - stim.halfSizeDeg;
stim.verticalSpanDeg   = stim.bottomEdgesDeg - stim.topEdgesDeg;
stim.horizontalSpanDeg = stim.rightEdgesDeg - stim.leftEdgesDeg;

if any(stim.verticalSpanDeg <= 0) || any(stim.horizontalSpanDeg <= 0)
    error('Aperture size must exceed the diameter of every dot.');
end
end

function abort = presentStimulus(dp, stim, isMoving, motionParams, stimDurationMs)
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

    if checkForEscape()
        abort = true;
        break;
    end
end
end

function abort = presentBlank(dp, durationMs)
abort = false;

Screen('FillRect', dp.wPtr, dp.bkColor);

if durationMs <= 0
    Screen('Flip', dp.wPtr);
    abort = checkForEscape();
    return;
end

numFrames = max(1, round((durationMs/1000) / dp.ifi));
vbl = Screen('Flip', dp.wPtr);

for frameIdx = 1:numFrames
    Screen('FillRect', dp.wPtr, dp.bkColor);
    vbl = Screen('Flip', dp.wPtr, vbl + 0.5 * dp.ifi);

    if checkForEscape()
        abort = true;
        break;
    end
end
end

function [xPosDeg, yPosDeg] = updatePositions(xPosDeg, yPosDeg, motionParams, stim, dp)
stepDeg = motionParams.speedDegPerSec / dp.frameRate;

switch lower(motionParams.direction)
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

function pressed = checkForEscape()
[keyIsDown, ~, keyCode] = KbCheck;
pressed = keyIsDown && keyCode(KbName('ESCAPE'));
end

function sizesDeg = generateDotSizes(numDots, targetMeanDeg, params, ratioCounts)
maxAttempts = 1000;
baseSizes = [repmat(params.smallSizeDeg,1,ratioCounts(1)), repmat(params.largeSizeDeg,1,ratioCounts(2))];

if numel(baseSizes) ~= numDots
    error('Ratio counts must sum to %d dots.', numDots);
end

scaleToTarget = targetMeanDeg / mean(baseSizes);
baseSizes = baseSizes * scaleToTarget;

jitterStdDeg = params.jitterStdRatio * targetMeanDeg;

for attempt = 1:maxAttempts
    noise = jitterStdDeg .* randn(1, numDots);
    noise = noise - mean(noise);  % 평균 보존

    sizesDeg = baseSizes + noise;
    sizesDeg = min(max(sizesDeg, params.minSizeDeg), params.maxSizeDeg);

    currentMean = mean(sizesDeg);
    delta = targetMeanDeg - currentMean;
    sizesDeg = sizesDeg + delta/numDots;

    if any(sizesDeg < params.minSizeDeg) || any(sizesDeg > params.maxSizeDeg)
        continue;
    end

    if abs(mean(sizesDeg) - targetMeanDeg) <= params.gToleranceDeg
        return;
    end
end

error('Could not generate dot sizes within tolerance. Adjust parameters.');
end

function [xPosDeg, yPosDeg] = initializeNonOverlappingPositions(numDots, dotSizeDeg, aperture, safetyMarginDeg)
xPosDeg = zeros(1, numDots);
yPosDeg = zeros(1, numDots);

for ii = 1:numDots
    currentSize = dotSizeDeg(ii);
    halfSize = currentSize / 2;

    xBounds = [aperture.leftDeg + halfSize, aperture.rightDeg - halfSize];
    yBounds = [aperture.topDeg + halfSize,  aperture.bottomDeg - halfSize];

    [xPosDeg(ii), yPosDeg(ii)] = sampleNonOverlappingPosition( ...
        xPosDeg(1:ii-1), yPosDeg(1:ii-1), dotSizeDeg(1:ii-1), currentSize, xBounds, yBounds, safetyMarginDeg);
end
end

function [xCandidate, yCandidate] = sampleNonOverlappingPosition(existingX, existingY, existingSizes, currentSize, xBounds, yBounds, safetyMarginDeg)
maxAttempts = 1000;
widthDeg  = xBounds(2) - xBounds(1);
heightDeg = yBounds(2) - yBounds(1);

if widthDeg <= 0 || heightDeg <= 0
    error('Aperture is too small for the requested dot size.');
end

for attempt = 1:maxAttempts
    xCandidate = xBounds(1) + rand * widthDeg;
    yCandidate = yBounds(1) + rand * heightDeg;

    if isempty(existingX)
        return;
    end

    dx = existingX - xCandidate;
    dy = existingY - yCandidate;
    minDist = (existingSizes + currentSize)/2 + safetyMarginDeg;

    if all((dx.^2 + dy.^2) >= (minDist.^2))
        return;
    end
end

error('Failed to place dot without overlap after %d attempts.', maxAttempts);
end

function xyPix = convertCentersToPixels(xDeg, yDeg, dp)
xyPix = [xDeg; yDeg] * dp.ppd;
xyPix(1,:) = xyPix(1,:) + dp.cx;
xyPix(2,:) = xyPix(2,:) + dp.cy;
end