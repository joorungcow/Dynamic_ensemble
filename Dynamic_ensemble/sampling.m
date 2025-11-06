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
dotParams.meanJitterDeg = 0.05;
dotParams.meanDiffTolerance = 0.005;
dotParams.gToleranceDeg = 0.001;
dotParams.jitterStdRatio = 0.15;
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

    trials = buildBlockTrials(blockLabel, comboLabels, diffLevels, freqRatios, signVals, methodPairs, repeatCount);
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

function trials = buildBlockTrials(blockLabel, comboLabels, diffLevels, freqRatios, signVals, methodPairs, repeatCount)
% Build the trial list for a single block according to balance rules and repeat count.

freqRatioLabels = {freqRatios.label};

baseCombos = struct('freqRatio', {}, 'diffLevel', {});
idx = 1;
for r = 1:numel(freqRatios)
    for d = 1:numel(diffLevels)
        for rep = 1:repeatCount
            baseCombos(idx).freqRatio = freqRatioLabels{r}; %#ok<AGROW>
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

        trialTemplate = composeTrialStruct(comboSeq{1}, augmentComboInfo(comboPool(1), signSeq(1), methodSeq{1}));
        trials = repmat(trialTemplate, 1, numTrials);

        for t = 1:numTrials
            comboInfo = augmentComboInfo(comboPool(t), signSeq(t), methodSeq{t});
            trials(t) = composeTrialStruct(comboSeq{t}, comboInfo);
        end
    case 'R'
        labels = [repmat({'SM'}, 1, numTrials/2), repmat({'MS'}, 1, numTrials/2)];
        labels = labels(randperm(numTrials));

        trialTemplate = composeTrialStruct(labels{1}, augmentComboInfo(comboPool(1), signSeq(1), methodSeq{1}));
        trials = repmat(trialTemplate, 1, numTrials);

        for t = 1:numTrials
            comboInfo = augmentComboInfo(comboPool(t), signSeq(t), methodSeq{t});
            trials(t) = composeTrialStruct(labels{t}, comboInfo);
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