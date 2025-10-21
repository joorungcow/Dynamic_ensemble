Screen('Preference','SkipSyncTests',1);

PsychDefaultSetup(2);         % 0~1 색 범위, 키맵 통일
AssertOpenGL;                 % PTB OpenGL 경로/드라이버 체크

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
%dp.aspectRatio = dp.resolution(2)/dp.resolution(1);

% 픽셀-각도 변환(ppd, pixels per degree)
dp.ppd = dp.resolution(1) / ((2*atan(dp.width/(2*dp.dist)))*180/pi);

rect = []; %if want FullScreen
%rect = [0 0 800 600];
[dp.wPtr, dp.wRect] = PsychImaging('OpenWindow', dp.screenNum, dp.bkColor, rect, [], [], 0);

% --- after openwindow : real params 값 update ---
Screen('BlendFunction', dp.wPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
dp.ifi        = Screen('GetFlipInterval', dp.wPtr);
dp.frameRate  = round(1/dp.ifi);        % 실제 주사율
dp.resolution = dp.wRect([3 4]);        % 실제 창 크기
[dp.cx, dp.cy]= RectCenter(dp.wRect);

%% dot configuration
numDots = 16;

% coherence 1: 모든 도트가 동일한 위 방향 속도를 가짐
speedPixPerFrame = 0.5;  % 프레임당 위쪽 이동 픽셀 수 (양수)

% Dot size control parameter
dotSizeParams.targetMean     = 12;   % 기본 평균 크기 (px)
dotSizeParams.meanJitter     = 3;    % trial 시작 시 평균을 얼마나 랜덤 가중할지 (±값)
dotSizeParams.minSize        = 4;    % 최소 크기 제한 (px)
dotSizeParams.maxSize        = 25;   % 최대 크기 제한 (px)
dotSizeParams.gTolerance     = 0.05;  % 실제 평균이 목표 평균에서 벗어날 수 있는 허용 오차(g 값)
dotSizeParams.jitterStdRatio = 0.35; % 개별 도트 크기 산출 시 사용할 표준편차 비율

% Randomize the target mean slightly and clip it within the specified bounds
targetMeanDotSize = dotSizeParams.targetMean + dotSizeParams.meanJitter * (2*rand - 1);
targetMeanDotSize = min(max(targetMeanDotSize, dotSizeParams.minSize), dotSizeParams.maxSize);

% Generate a set of dot sizes that meet the specified conditions
dotSize = generateDotSizes(numDots, targetMeanDotSize, dotSizeParams);

% Define the rectangular area centered on the screen
areaWidth  = 400;  % 
areaHeight = 300; % 
areaRect   = CenterRectOnPoint([0 0 areaWidth areaHeight], dp.cx, dp.cy);

% Calculate the height of the aperture based on the defined rectangle
apertureHeight = areaRect(4) - areaRect(2);

% Initialize non-overlapping positions for the dots within the defined area
[xPos, yPos] = initializeNonOverlappingPositions(numDots, dotSize, areaRect);

% Calculate boundary values for wrapping within the area
topEdges     = areaRect(2) + dotSize/2;
bottomEdges  = areaRect(4) - dotSize/2;
verticalSpan = bottomEdges - topEdges;

% Check if the vertical span is valid; if not, throw an error
if any(verticalSpan <= 0)
    error('Aperture height (%.1f px) must exceed every dot diameter.', apertureHeight);
end

% 색: 전체 동일 (0~1 범위, PsychDefaultSetup(2) 기준)
dotColor = [1 1 1];

% 애니메이션 루프 설정
animationDuration = 5; % 초
numFrames = round(animationDuration / dp.ifi);

vbl = Screen('Flip', dp.wPtr);
for frameIdx = 1:numFrames
    %#ok<NASGU> % frameIdx는 디버그 시 참조 가능
    % 위치 업데이트 (위쪽으로 이동)
    yPos = yPos - speedPixPerFrame;

    % 위쪽 경계에 닿은 도트는 초기 샘플링 위치 관계를 유지한 채 아래로 wrap
    yPos = topEdges + mod(yPos - topEdges, verticalSpan);

    % 좌표는 ScreenDrawDots에 2xN 행렬로 전달
    xy = [xPos; yPos];

    % 도트 그리기 (절대좌표 사용을 위해 center 인수로 [0 0] 전달)
    Screen('DrawDots', dp.wPtr, xy, dotSize, dotColor, [0 0], 2);

    vbl = Screen('Flip', dp.wPtr, vbl + 0.5 * dp.ifi);

    % Esc 키 입력 시 종료
    [keyIsDown, ~, keyCode] = KbCheck;
    if keyIsDown && keyCode(KbName('ESCAPE'))
        break;
    end
end

WaitSecs(0.5);

Screen('CloseAll');

%% --- Local functions ---
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
% 기존 도트들과 겹치지 않는 위치를 무작위로 플링
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

error('Failed to place dot without overlap after %d attempts.', maxAttempts);
end