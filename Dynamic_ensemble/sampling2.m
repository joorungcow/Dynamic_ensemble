function sampling2()
%SAMPLING2 Generate counterbalanced trial plan (main + practice) and save to file.
%
%   sampling2() prompts for a subject identifier and block order string
%   (permutation of the letters M, S, R) and creates a plan struct whose
%   total number of trials is determined by the stimulus combination
%   factors (combo × frequency ratio × diff level × repeat count). The
%   plan is saved to ./sampling/<subID>_sampling2.mat.
%
%   This version also prepares practice blocks using a separate repeat
%   count (repeatCountP) so that the practice and main sessions share the
%   same sampling logic.

repeatCount  = 30; % 본 실험에서 각 조합 반복 횟수
repeatCountP = 1;  % 연습 시행에서 각 조합 반복 횟수 (기본 1 → 72 조건 1회씩)

subID = strtrim(input('Enter subject ID (e.g., sub01): ', 's'));
if isempty(subID)
    error('Subject ID must be provided.');
end

blockOrderStr = upper(strtrim(input('Enter block order (SMR, SRM, MSR, MRS, RSM, or RMS): ', 's')));
blockLetters = parseBlockOrder(blockOrderStr);

expo = 0.76;
diffLevels = [0.06 0.12 0.18 0.24 0.30 0.36];

freqRatios = struct(...
    'label', {'6:2','5:3','3:5','2:6'}, ...
    'counts', { [6 2], [5 3], [3 5], [2 6] });

plan = struct();
plan.meta = struct();
plan.meta.subID = subID;
plan.meta.blockOrderStr = blockOrderStr;
plan.meta.repeatCount = repeatCount;
plan.meta.expo = expo;
plan.meta.diffLevels = diffLevels;
plan.meta.freqRatioLabels = {freqRatios.label};
plan.meta.equalFreqLabel = '4:4';

%% Display settings saved with the plan
displayParams = struct();
displayParams.screenNum = [];
displayParams.dist = 55;
displayParams.width = 60;
displayParams.bkColor = [0.5 0.5 0.5];
displayParams.textColor = [1 1 1];
displayParams.textFont = 'Malgun Gothic';
displayParams.textSize = 24;
displayParams.textLineSpacingMultiplier = 3;
displayParams.responseInstructions = {
    '1번 키: T1 평균이 더 큽니다';
    '2번 키: T2 평균이 더 큽니다'
};
displayParams.skipChecks = 1;

%% Dot configuration saved with the plan
dotParams = struct();
dotParams.smallSizeDeg = 0.7;
dotParams.largeSizeDeg = 1.3;
dotParams.minSizeDeg = 0.4;
dotParams.maxSizeDeg = 1.8;
dotParams.meanJitterDeg = 0.01;
dotParams.meanDiffTolerance = 0.001;
dotParams.gToleranceDeg = 0.001;
dotParams.jitterStdRatio = 0.05;
dotParams.perceptualExponent = expo;
dotParams.meanDiffLevels = diffLevels;
dotParams.safetyMarginDeg = 0.05;
dotParams.equalFreqCounts = [4 4];
dotParams.useHalfRange = false;

%% Motion configuration saved with the plan
motionParams = struct();
motionParams.direction = 'up';
motionParams.speedDegPerSec = 1.5;

%% Timing configuration saved with the plan (milliseconds unless noted)
timingParams.postTrialMs = 1000;
timingParams.breakDurationSec = 20;

%% Grid layout configuration saved with the plan
gridConfig = struct();
gridConfig.rows = 6;
gridConfig.cols = 6;
gridConfig.windowWidthDeg = 12;
gridConfig.windowHeightDeg = 12;
gridConfig.maxJitterDeg = 0.12;
gridConfig.maxAttemptsPerCell = 50;
gridConfig.safetyMarginDeg = dotParams.safetyMarginDeg;
gridConfig.outerPaddingDeg = dotParams.maxSizeDeg/2 + dotParams.safetyMarginDeg;

plan.display = displayParams;
plan.dotParams = dotParams;
plan.motionParams = motionParams;
plan.timingParams = timingParams;
plan.gridConfig = gridConfig;

gridLayout = buildCentralGridLayout(gridConfig);
plan.gridLayout = gridLayout;

signVals = [+1, -1];
methodPairs = {'FREQ_vs_EQ','EQ_vs_FREQ'};

plan.blocks         = repmat(struct('label', '', 'trials', []), 1, numel(blockLetters));
plan.practiceBlocks = repmat(struct('label', '', 'trials', []), 1, numel(blockLetters));

totalTrialsMain     = 0;
totalTrialsPractice = 0;

for b = 1:numel(blockLetters)
    letter = blockLetters{b};
    switch letter
        case 'M'
            blockLabel = 'MM';
            comboLabels = {'MM'};
        case 'S'
            blockLabel = 'SS';
            comboLabels = {'SS'};
        case 'R'
            blockLabel = 'R';
            comboLabels = {'SM','MS'};
        otherwise
            error('Unknown block letter %s.', letter);
    end

    trialsPractice = buildBlockTrials(blockLabel, comboLabels, diffLevels, freqRatios, ...
        signVals, methodPairs, repeatCountP, dotParams, gridLayout, gridConfig);
    plan.practiceBlocks(b).label  = blockLabel;
    plan.practiceBlocks(b).trials = trialsPractice;
    totalTrialsPractice = totalTrialsPractice + numel(trialsPractice);

    trialsMain = buildBlockTrials(blockLabel, comboLabels, diffLevels, freqRatios, ...
        signVals, methodPairs, repeatCount, dotParams, gridLayout, gridConfig);
    plan.blocks(b).label  = blockLabel;
    plan.blocks(b).trials = trialsMain;
    totalTrialsMain = totalTrialsMain + numel(trialsMain);
end

plan.meta.NMain      = totalTrialsMain;
plan.meta.NPractice  = totalTrialsPractice;
plan.meta.N          = totalTrialsMain;        % test13m 기존 코드 호환용
plan.meta.repeatMain = repeatCount;
plan.meta.repeatPrac = repeatCountP;

if ~exist('sampling', 'dir')
    mkdir('sampling');
end

saveFile = fullfile('sampling', sprintf('%s_sampling2.mat', subID));
save(saveFile, 'plan');
fprintf('Plan saved to %s\n', saveFile);
end

function trials = buildBlockTrials(blockLabel, comboLabels, diffLevels, freqRatios, signVals, methodPairs, repeatCount, dotParams, gridLayout, gridConfig)
% Build the trial list for a single block according to balance rules and repeat count.

freqRatioLabels = {freqRatios.label};

baseCombos = struct('freqRatio', {}, 'diffLevel', {});
idx = 1;
for r = 1:numel(freqRatios)
    for d = 1:numel(diffLevels)
        for rep = 1:repeatCount
            baseCombos(idx).freqRatio = freqRatioLabels{r};
            baseCombos(idx).diffLevel = diffLevels(d);
            idx = idx + 1;
        end
    end
end

numTrials = numel(baseCombos);
if strcmp(blockLabel, 'R') && mod(numTrials, 2) ~= 0
    error('Repeat count must yield an even number of trials for the mixed block.');
end

signSeq = repmat(signVals, 1, ceil(numTrials / numel(signVals)));
signSeq = signSeq(1:numTrials);
signSeq = signSeq(randperm(numTrials));

methodSeq = repmat(methodPairs, 1, ceil(numTrials / numel(methodPairs)));
methodSeq = methodSeq(1:numTrials);
methodSeq = methodSeq(randperm(numTrials));

comboPool = baseCombos(randperm(numTrials));

switch blockLabel
    case {'MM','SS'}
        comboSeq = repmat(comboLabels, 1, ceil(numTrials / numel(comboLabels)));
        comboSeq = comboSeq(1:numTrials);
for attempt = 1:maxAttempts
    noise = noiseStdPS .* randn(size(basePS));
    noise = noise - mean(noise);

    candidatePS = basePS + noise;
    candidatePS = min(max(candidatePS, minPS), maxPS);
    candidatePS = adjustMeanPSWithinBounds(candidatePS, targetMeanPS, minPS, maxPS, tolerancePS);

    candidateDeg = fromPS(candidatePS, expo);
    candidateDeg = min(max(candidateDeg, params.minSizeDeg), params.maxSizeDeg);
    candidatePS = toPS(candidateDeg, expo);
    candidatePS = adjustMeanPSWithinBounds(candidatePS, targetMeanPS, minPS, maxPS, tolerancePS);

    currentDiff = abs(mean(candidatePS) - targetMeanPS);
    if currentDiff < bestDiff
        bestDiff = currentDiff;
        bestPS = candidatePS;
    end

    if currentDiff <= tolerancePS
        bestPS = candidatePS;
        break;
    end
end

finalPS = bestPS;
sizesDeg = fromPS(finalPS, expo);
sizesDeg = min(max(sizesDeg, params.minSizeDeg), params.maxSizeDeg);
finalPS = toPS(sizesDeg, expo);
finalMeanPS = mean(finalPS);
end

function tolerancePS = computePSTolerance(tolDeg, targetMeanPS, expo)
targetDeg = fromPS(targetMeanPS, expo);
tolerancePS = abs(toPS(targetDeg + tolDeg, expo) - targetMeanPS);
end

function projected = adjustMeanPSWithinBounds(psValues, targetMeanPS, minPS, maxPS, tolerance)
numVals = numel(psValues);
maxIter = 50;
projected = psValues;
for iter = 1:maxIter
    meanDiff = targetMeanPS - mean(projected);
    if abs(meanDiff) <= tolerance
        break;
    end

    totalDiff = meanDiff * numVals;
    if totalDiff > 0
        adjustable = find(projected < maxPS - eps);
        if isempty(adjustable)
            break;
        end
        slack = maxPS - projected(adjustable);
        totalSlack = sum(slack);
        if totalSlack <= 0
            break;
        end
        factor = min(1, totalDiff / totalSlack);
        increments = factor .* slack;
        projected(adjustable) = projected(adjustable) + increments;
    else
        adjustable = find(projected > minPS + eps);
        if isempty(adjustable)
            break;
        end
        slack = projected(adjustable) - minPS;
        totalSlack = sum(slack);
        if totalSlack <= 0
            break;
        end
        factor = min(1, -totalDiff / totalSlack);
        decrements = factor .* slack;
        projected(adjustable) = projected(adjustable) - decrements;
    end
    projected = min(max(projected, minPS), maxPS);
end
end