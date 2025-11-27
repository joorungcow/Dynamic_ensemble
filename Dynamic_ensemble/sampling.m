function sampling()
%SAMPLING Generate counterbalanced trial plan and save to file.
%
%   sampling() prompts for a subject identifier and block order string
%   (permutation of the letters M, S, R) and creates a plan struct whose
%   total number of trials is determined by the stimulus combination
%   factors (combo × frequency ratio × diff level × repeat count). The
%   plan is saved to ./sampling/<subID>_sampling.mat.
%
%   The script can be converted to a function handle for batch
%   generation, but the default usage is to run it interactively.

repeatCount = 1; % Set the repeat count per combination here (edit as needed).

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
    '왼쪽 방향키: T1 평균이 더 큽니다';
    '오른쪽 방향키: T2 평균이 더 큽니다'
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
timingParams = struct();
timingParams.fixationMs = 300;
timingParams.stimDurationMs = 500;
timingParams.isiDurationMs = 1000;
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

plan.blocks = repmat(struct('label', '', 'trials', []), 1, numel(blockLetters));

totalTrials = 0;

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

    trials = buildBlockTrials(blockLabel, comboLabels, diffLevels, freqRatios, signVals, methodPairs, repeatCount, dotParams, gridLayout, gridConfig);
    plan.blocks(b).label = blockLabel;
    plan.blocks(b).trials = trials;
    totalTrials = totalTrials + numel(trials);
end

plan.meta.N = totalTrials;

if ~exist('sampling', 'dir')
    mkdir('sampling');
end

saveFile = fullfile('sampling', sprintf('%s_sampling.mat', subID));
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
        comboSeq = comboSeq(randperm(numTrials));

        firstInfo = augmentComboInfo(comboPool(1), signSeq(1), methodSeq{1});
        firstTrial = finalizeTrialStimulus(composeTrialStruct(comboSeq{1}, firstInfo), dotParams, gridLayout, gridConfig);
        trials = repmat(firstTrial, 1, numTrials);

        for t = 1:numTrials
            if t == 1
                continue;
            end
            comboInfo = augmentComboInfo(comboPool(t), signSeq(t), methodSeq{t});
            baseTrial = composeTrialStruct(comboSeq{t}, comboInfo);
            trials(t) = finalizeTrialStimulus(baseTrial, dotParams, gridLayout, gridConfig);
        end
    case 'R'
        labels = [repmat({'SM'}, 1, numTrials/2), repmat({'MS'}, 1, numTrials/2)];
        labels = labels(randperm(numTrials));

        firstInfo = augmentComboInfo(comboPool(1), signSeq(1), methodSeq{1});
        firstTrial = finalizeTrialStimulus(composeTrialStruct(labels{1}, firstInfo), dotParams, gridLayout, gridConfig);
        trials = repmat(firstTrial, 1, numTrials);

        for t = 1:numTrials
            if t == 1
                continue;
            end
            comboInfo = augmentComboInfo(comboPool(t), signSeq(t), methodSeq{t});
            baseTrial = composeTrialStruct(labels{t}, comboInfo);
            trials(t) = finalizeTrialStimulus(baseTrial, dotParams, gridLayout, gridConfig);
        end
    otherwise
        error('Unsupported block label %s', blockLabel);
end
end

function trial = composeTrialStruct(comboLabel, comboInfo)
trial = struct();
trial.comboLabel = comboLabel;
trial.motionCombo = comboLabel;
trial.methodPair = comboInfo.methodPair;
trial.freqRatio = comboInfo.freqRatio;
trial.equalFreq = true;
trial.diffLevel = comboInfo.diffLevel;
trial.signS = comboInfo.signS;
trial.isMoving = [comboLabel(1) == 'M', comboLabel(2) == 'M'];
trial.numDots = 8; % fixed (6:2 vs 4:4 etc.)
end

function comboInfo = augmentComboInfo(comboInfo, signVal, methodPair)
comboInfo.signS = signVal;
comboInfo.methodPair = methodPair;
end

function letters = parseBlockOrder(orderStr)
validOrders = {'SMR','SRM','MSR','MRS','RSM','RMS'};

if numel(orderStr) ~= 3
    error('Block order must contain exactly three letters.');
end

if ~ismember(orderStr, validOrders)
    error('Block order must be one of: %s.', char(strjoin(validOrders, ', ')));
end

letters = cellfun(@(c) string(c), cellstr(orderStr.'));
letters = cellfun(@char, letters, 'UniformOutput', false);
end

function trial = finalizeTrialStimulus(trial, params, layout, gridConfig)
[stimPair, stimInfo] = buildStimulusForTrial(trial, params, layout, gridConfig);
trial.t1Stimulus = stimPair.t1;
trial.t2Stimulus = stimPair.t2;
trial.stimInfo = stimInfo;
trial.t1MeanPS = stimInfo.t1MeanPS;
trial.t2MeanPS = stimInfo.t2MeanPS;
trial.t1MeanDeg = stimInfo.t1MeanDeg;
trial.t2MeanDeg = stimInfo.t2MeanDeg;
trial.jitterPS_applied = stimInfo.jitterPS;
trial.targetPSMeans = stimInfo.targetPSMeans;
trial.psTargetsAfterJitter = stimInfo.psTargetsAfterJitter;
trial.largerSide = stimInfo.largerSide;
end

function [stimPair, info] = buildStimulusForTrial(trialSpec, params, layout, gridConfig)
expo = params.perceptualExponent;
freqCounts = getFrequencyRatioCounts(trialSpec.freqRatio);
eqCounts = params.equalFreqCounts;

freqBase = make_frequency_set(freqCounts, params);
eqBase = make_equalFreq_set(eqCounts, params);

if strcmp(trialSpec.methodPair, 'FREQ_vs_EQ')
    baseT1 = freqBase;
    baseT2 = eqBase;
else
    baseT1 = eqBase;
    baseT2 = freqBase;
end

ratioTarget = 1 + trialSpec.signS * trialSpec.diffLevel;
mu1Target = meanPS(baseT1, expo);
mu2Target = mu1Target * ratioTarget;

jitterRangePS = computeJitterRangePS(mu1Target, mu2Target, params);
[mu1Target, mu2Target] = jitterMeanPairPS(mu1Target, mu2Target, jitterRangePS, trialSpec.diffLevel, trialSpec.signS);
mu1TargetJittered = mu1Target;
mu2TargetJittered = mu2Target;

[t1SizesDeg, t1MeanPS] = synthesizeSetFromBase(baseT1, mu1Target, params);
ratioTarget = 1 + trialSpec.signS * trialSpec.diffLevel;
mu2Target = t1MeanPS * ratioTarget;
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
info.expoUsed = expo;
info.jitterPS = jitterRangePS;
info.targetPSMeans = [mu1TargetJittered, mu2TargetJittered];
info.psTargetsAfterJitter = [mu1TargetJittered, mu2Target];
info.largerSide = ternary(trialSpec.signS >= 0, 'T2', 'T1');
end

function [sizesDeg, finalMeanPS] = synthesizeSetFromBase(baseSizesDeg, targetMeanPS, params)
expo = params.perceptualExponent;
basePS = toPS(baseSizesDeg, expo);
scale = targetMeanPS / mean(basePS);
basePS = basePS * scale;

minPS = toPS(params.minSizeDeg, expo);
maxPS = toPS(params.maxSizeDeg, expo);
tolerancePS = computePSTolerance(params.gToleranceDeg, targetMeanPS, expo);
if tolerancePS <= 0
    tolerancePS = 1e-5;
end

maxAttempts = 1000;
noiseStdPS = params.jitterStdRatio * targetMeanPS;
bestPS = basePS;
bestDiff = abs(mean(basePS) - targetMeanPS);

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

if abs(mean(projected) - targetMeanPS) > tolerance
    offset = targetMeanPS - mean(projected);
    projected = projected + offset;
    projected = min(max(projected, minPS), maxPS);
end
end

function jitterPS = computeJitterRangePS(mu1PS, mu2PS, params)
expo = params.perceptualExponent;
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
mu1Adj = max(eps, mu1PS + jitter1);
mu2Adj = max(eps, mu2PS + jitter2);

ratioTarget = 1 + signS * diffLevel;
projectedMu2 = mu1Adj * ratioTarget;
projectedMu1 = mu2Adj / ratioTarget;
if abs(projectedMu2 - mu2Adj) <= abs(projectedMu1 - mu1Adj)
    mu2Adj = max(eps, projectedMu2);
else
    mu1Adj = max(eps, projectedMu1);
end
end

function sizes = make_frequency_set(counts, params)
sizes = [repmat(params.smallSizeDeg, 1, counts(1)), repmat(params.largeSizeDeg, 1, counts(2))];
sizes = sizes(randperm(numel(sizes)));
end

function sizes = make_equalFreq_set(counts, params)
base = [repmat(params.smallSizeDeg, 1, counts(1)), repmat(params.largeSizeDeg, 1, counts(2))];
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
stimSet.halfSizeDeg = stimSet.dotSizeDeg / 2;

stimSet.topEdgesDeg = stimSet.outerAperture.topDeg + stimSet.halfSizeDeg;
stimSet.bottomEdgesDeg = stimSet.outerAperture.bottomDeg - stimSet.halfSizeDeg;
stimSet.leftEdgesDeg = stimSet.outerAperture.leftDeg + stimSet.halfSizeDeg;
stimSet.rightEdgesDeg = stimSet.outerAperture.rightDeg - stimSet.halfSizeDeg;

stimSet.verticalSpanDeg = stimSet.bottomEdgesDeg - stimSet.topEdgesDeg;
stimSet.horizontalSpanDeg = stimSet.rightEdgesDeg - stimSet.leftEdgesDeg;

stimSet.innerTopEdgesDeg = stimSet.innerAperture.topDeg + stimSet.halfSizeDeg;
stimSet.innerBottomEdgesDeg = stimSet.innerAperture.bottomDeg - stimSet.halfSizeDeg;
stimSet.innerLeftEdgesDeg = stimSet.innerAperture.leftDeg + stimSet.halfSizeDeg;
stimSet.innerRightEdgesDeg = stimSet.innerAperture.rightDeg - stimSet.halfSizeDeg;

if any(stimSet.innerBottomEdgesDeg <= stimSet.innerTopEdgesDeg) || any(stimSet.innerRightEdgesDeg <= stimSet.innerLeftEdgesDeg)
    error('Inner aperture must be larger than dot diameters.');
end

if any(stimSet.verticalSpanDeg <= 0) || any(stimSet.horizontalSpanDeg <= 0)
    error('Outer aperture size must exceed the diameter of every dot.');
end
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

        leftLimit = baseX - (layout.innerAperture.leftDeg + halfSize + gridConfig.safetyMarginDeg);
        rightLimit = (layout.innerAperture.rightDeg - halfSize - gridConfig.safetyMarginDeg) - baseX;
        topLimit = baseY - (layout.innerAperture.topDeg + halfSize + gridConfig.safetyMarginDeg);
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

            if candidateX - halfSize < layout.innerAperture.leftDeg + gridConfig.safetyMarginDeg
                continue;
            end
            if candidateX + halfSize > layout.innerAperture.rightDeg - gridConfig.safetyMarginDeg
                continue;
            end
            if candidateY - halfSize < layout.innerAperture.topDeg + gridConfig.safetyMarginDeg
                continue;
            end
            if candidateY + halfSize > layout.innerAperture.bottomDeg - gridConfig.safetyMarginDeg
                continue;
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
        error('Could not place dot %d within constraints.', dotIdx);
    end
end
end

function layout = buildCentralGridLayout(gridConfig)
inner.centerDeg = [0, 0];
inner.widthDeg = gridConfig.windowWidthDeg;
inner.heightDeg = gridConfig.windowHeightDeg;
inner = updateApertureEdges(inner);

outer = inner;
outer.widthDeg = inner.widthDeg + 2 * gridConfig.outerPaddingDeg;
outer.heightDeg = inner.heightDeg + 2 * gridConfig.outerPaddingDeg;
outer = updateApertureEdges(outer);

cellWidthDeg = inner.widthDeg / gridConfig.cols;
cellHeightDeg = inner.heightDeg / gridConfig.rows;

colCenters = inner.leftDeg + (0.5:1:gridConfig.cols-0.5) * cellWidthDeg;
rowCenters = inner.topDeg + (0.5:1:gridConfig.rows-0.5) * cellHeightDeg;
[gridX, gridY] = meshgrid(colCenters, rowCenters);

layout.innerAperture = inner;
layout.outerAperture = outer;
layout.gridX = gridX(:);
layout.gridY = gridY(:);
layout.cellWidthDeg = cellWidthDeg;
layout.cellHeightDeg = cellHeightDeg;
end

function aperture = updateApertureEdges(aperture)
halfWidth = aperture.widthDeg / 2;
halfHeight = aperture.heightDeg / 2;

aperture.leftDeg = aperture.centerDeg(1) - halfWidth;
aperture.rightDeg = aperture.centerDeg(1) + halfWidth;
aperture.topDeg = aperture.centerDeg(2) - halfHeight;
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