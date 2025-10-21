Screen('Preference','SkipSyncTests',1);

PsychDefaultSetup(2);         % 0~1 색 범위, 키맵 통일
AssertOpenGL;                 % PTB OpenGL 경로/드라이버 체크
KbName('UnifyKeyNames');

%% Display settings (params)
dp.screenNum = max(Screen('Screens'));

dp.dist   = 55;   
dp.width  = 60;   
dp.bkColor   = 0.5;            
dp.textColor = [0 0 0];        

% 해상도/주사율/비율/ppd
d = Screen('Resolution', dp.screenNum);
dp.resolution  = [d.width, d.height];        % [W H] px
dp.frameRate   = d.hz;                        % Hz
%dp.aspectRatio = dp.resolution(2)/dp.resolution(1);

% 픽셀-각도 변환(ppd, pixels per degree)
dp.ppd = dp.resolution(1) / ((2*atan(dp.width/(2*dp.dist)))*180/pi);

rect = []; %if want FullScreen
%rect = [0 0 800 600];

try
    [dp.wPtr, dp.wRect] = PsychImaging('OpenWindow', dp.screenNum, dp.bkColor, rect, [], [], 0);

    % --- after openwindow : real params 값 update ---
    Screen('BlendFunction', dp.wPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    dp.ifi        = Screen('GetFlipInterval', dp.wPtr);
    dp.frameRate  = round(1/dp.ifi);        % 실제 주사율
    dp.resolution = dp.wRect([3 4]);        % 실제 창 크기
    [dp.cx, dp.cy]= RectCenter(dp.wRect);

    %% dot configuration
    numDots = 16;

    % Motion configuration
    motionParams.direction        = 'up';   % 'up', 'down', 'left', 'right'
    motionParams.speedPixPerFrame = 0.5;    % 이동 속도 (픽셀/프레임)

    % Timing configuration (ms 단위)
    timingParams.stimDurationMs    = 500;   % T1/T2 자극 제시 시간
    timingParams.isiDurationMs     = 200;   % T1-T2 사이 공백
    timingParams.postTrialMs       = 500;   % T2 이후 공백

    % 자극 조합 반복 횟수 (4 * n 트라이얼 수행)
    comboRepeats = 1;  % n 값을 원하는 만큼 변경

    % Dot size control parameter
    dotSizeParams.targetMean     = 12;   % 기본 평균 크기 (px)
    dotSizeParams.meanJitter     = 3;    % trial 시작 시 평균을 얼마나 랜덤 가중할지 (±값)
    dotSizeParams.minSize        = 5;    % 최소 크기 제한 (px)
    dotSizeParams.maxSize        = 10;   % 최대 크기 제한 (px)
    dotSizeParams.gTolerance     = 0.01; % 실제 평균이 목표 평균에서 벗어날 수 있는 허용 오차(g 값)
    dotSizeParams.jitterStdRatio = 0.35; % 개별 도트 크기 산출 시 사용할 표준편차 비율

    % Define the rectangular area centered on the screen
    areaWidth  = 600;
    areaHeight = 500;
    areaRect   = CenterRectOnPoint([0 0 areaWidth areaHeight], dp.cx, dp.cy);

    % 색: 전체 동일 (0~1 범위, PsychDefaultSetup(2) 기준)
    dotColor = [1 1 1];

    % 초기 화면을 배경으로 클리어
    Screen('FillRect', dp.wPtr, dp.bkColor);
    Screen('Flip', dp.wPtr);

    stimCombos = {'MM','SM','MS','SS'};
    abortExperiment = false;

    for repIdx = 1:comboRepeats
        for comboIdx = 1:numel(stimCombos)
            if abortExperiment
                break;
            end

            comboLabel = stimCombos{comboIdx};
            t1IsMoving = comboLabel(1) == 'M';
            t2IsMoving = comboLabel(2) == 'M';

            targetMeanDotSize = dotSizeParams.targetMean + dotSizeParams.meanJitter * (2*rand - 1);
            targetMeanDotSize = min(max(targetMeanDotSize, dotSizeParams.minSize), dotSizeParams.maxSize);

            stim1 = createStimulusStruct(numDots, targetMeanDotSize, dotSizeParams, areaRect);
            abortExperiment = presentStimulus(dp, stim1, t1IsMoving, motionParams, dotColor, timingParams.stimDurationMs);
            if abortExperiment
                break;
            end

            abortExperiment = presentBlank(dp, timingParams.isiDurationMs);
            if abortExperiment
                break;
            end

            targetMeanDotSize = dotSizeParams.targetMean + dotSizeParams.meanJitter * (2*rand - 1);
            targetMeanDotSize = min(max(targetMeanDotSize, dotSizeParams.minSize), dotSizeParams.maxSize);

            stim2 = createStimulusStruct(numDots, targetMeanDotSize, dotSizeParams, areaRect);
            abortExperiment = presentStimulus(dp, stim2, t2IsMoving, motionParams, dotColor, timingParams.stimDurationMs);
            if abortExperiment
                break;
            end

            abortExperiment = presentBlank(dp, timingParams.postTrialMs);
        end

        if abortExperiment
            break;
        end
    end

    presentBlank(dp, timingParams.postTrialMs);

    Screen('CloseAll');
catch ME
    Screen('CloseAll');
    rethrow(ME);
end

%% --- Local functions ---
function stim = createStimulusStruct(numDots, targetMean, params, areaRect)
dotSize = generateDotSizes(numDots, targetMean, params);
[xPos, yPos] = initializeNonOverlappingPositions(numDots, dotSize, areaRect);

stim.dotSize = dotSize;
stim.xPos = xPos;
stim.yPos = yPos;
stim.areaRect = areaRect;
stim.topEdges    = areaRect(2) + dotSize/2;
stim.bottomEdges = areaRect(4) - dotSize/2;
stim.leftEdges   = areaRect(1) + dotSize/2;
stim.rightEdges  = areaRect(3) - dotSize/2;
stim.verticalSpan   = stim.bottomEdges - stim.topEdges;
stim.horizontalSpan = stim.rightEdges - stim.leftEdges;

if any(stim.verticalSpan <= 0) || any(stim.horizontalSpan <= 0)
    error('Aperture size must exceed the diameter of every dot.');
end
end

function abort = presentStimulus(dp, stim, isMoving, motionParams, dotColor, stimDurationMs)
abort = false;
numFrames = max(1, round((stimDurationMs/1000) / dp.ifi));

Screen('FillRect', dp.wPtr, dp.bkColor);
vbl = Screen('Flip', dp.wPtr);

for frameIdx = 1:numFrames
    if isMoving && frameIdx > 1
        [stim.xPos, stim.yPos] = updatePositions(stim.xPos, stim.yPos, motionParams, stim);
    end

    Screen('FillRect', dp.wPtr, dp.bkColor);
    xy = [stim.xPos; stim.yPos];
    Screen('DrawDots', dp.wPtr, xy, stim.dotSize, dotColor, [0 0], 2);
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

function [xPos, yPos] = updatePositions(xPos, yPos, motionParams, stim)
speed = motionParams.speedPixPerFrame;
switch lower(motionParams.direction)
    case 'up'
        yPos = yPos - speed;
        yPos = stim.topEdges + mod(yPos - stim.topEdges, stim.verticalSpan);
    case 'down'
        yPos = yPos + speed;
        yPos = stim.topEdges + mod(yPos - stim.topEdges, stim.verticalSpan);
    case 'left'
        xPos = xPos - speed;
        xPos = stim.leftEdges + mod(xPos - stim.leftEdges, stim.horizontalSpan);
    case 'right'
        xPos = xPos + speed;
        xPos = stim.leftEdges + mod(xPos - stim.leftEdges, stim.horizontalSpan);
    otherwise
        error('Unknown motion direction: %s', motionParams.direction);
end
end

function pressed = checkForEscape()
[keyIsDown, ~, keyCode] = KbCheck;
pressed = keyIsDown && keyCode(KbName('ESCAPE'));
end

function sizes = generateDotSizes(numDots, targetMean, params)
% 개별 도트 크기를 생성하되 평균이 gTolerance 이내로 유지되도록 반복
maxAttempts = 1000;
jitterStd = params.jitterStdRatio * targetMean;

for attempt = 1:maxAttempts
    sizes = targetMean + jitterStd .* randn(1, numDots);
    sizes = min(max(sizes, params.minSize), params.maxSize);
    sizes = round(sizes);

    if abs(mean(sizes) - targetMean) <= params.gTolerance
        return;
    end
end

error('Could not generate dot sizes within tolerance. Adjust parameters.');
end

function [xPos, yPos] = initializeNonOverlappingPositions(numDots, dotSize, areaRect)
% 영역 내에서 도트끼리 겹치지 않도록 초기 위치를 배정
xPos = zeros(1, numDots);
yPos = zeros(1, numDots);

for ii = 1:numDots
    currentSize = dotSize(ii);
     halfSize = currentSize / 2;
    xBounds = [areaRect(1) + halfSize, areaRect(3) - halfSize];
    yBounds = [areaRect(2) + halfSize, areaRect(4) - halfSize];

    [xPos(ii), yPos(ii)] = sampleNonOverlappingPosition( ...
        xPos(1:ii-1), yPos(1:ii-1), dotSize(1:ii-1), currentSize, xBounds, yBounds);
end
end

function [xCandidate, yCandidate] = sampleNonOverlappingPosition(existingX, existingY, existingSizes, currentSize, xBounds, yBounds)
% 기존 도트들과 겹치지 않는 위치를 무작위
maxAttempts = 1000;
width = xBounds(2) - xBounds(1);
height = yBounds(2) - yBounds(1);

if width <= 0 || height <= 0
    error('Aperture is too small for the requested dot size.');
end

for attempt = 1:maxAttempts
    xCandidate = xBounds(1) + rand * width;
    yCandidate = yBounds(1) + rand * height;

    if isempty(existingX)
        return;
    end

    dx = existingX - xCandidate;
    dy = existingY - yCandidate;
    minDist = (existingSizes + currentSize) / 2;

    if all((dx.^2 + dy.^2) >= (minDist.^2))
        return;
    end
end
end