clear variables
lab = 1;
params.debugFullDots = 0;
params.forcedBreak = 1;

if ~lab
    Screen('Preference', 'SkipSyncTests', 1);
end

% create random sequence
rng('shuffle'); % Generate a random seed based on the current time
params.randSeed = rng; % Save it

screens = Screen('Screens');
dp.screenNum = max(screens);
dp.bkColor = [.5 .5 .5];
dp.width = 60;
dp.dist = 55;

subID = input('Subject ID: ', 's');
params.subNum = input('Subject Number: ');
params.gender = input('Gender(male 1 female 0 others 9): ');
params.age = input('Age: ');
params.hand = input('Hand(right 1 left 0): ');

dp = OpenWindow(dp);
ListenChar(0);
HideCursor;

% to linearize the monitor
if lab
    load('BenQ_20220808.mat'); 
    gammaTable = repmat(luminanceGamma',1,3);
    oldGamma = Screen('LoadNormalizedGammaTable', dp.wPtr, gammaTable);
end

kb.useKbQueueCheck = 1;
kb.keyList = {'space','LeftArrow','rightArrow','ESCAPE'};
kb = init_keyboard(kb);

params.sSize = 1.4; % in visual angle

%to use psychometric function Parademes
params.nCon = 4; %random, symmetry each *2
params.NumTrials = 25 * params.nCon; 

% params.dist2Fix = 6;
% params.duration = .5; % .25

if rem(params.subNum,2) == 0
    params.duration = [repmat([.5;.25],12/2,1)];
else
    params.duration = [repmat([.25;.5],12/2,1)];
end

params.fixColor = [0 0 0];
params.fixSize = .7*dp.ppd; %origin: .5

center.cx = dp.cx;
center.cy = dp.cy;

%start exp
expVer = 1;
try
    iRun = 0;
    practice_symmetry(params,center,dp,kb,subID,expVer,iRun);
    for iRun = 1:10
        run_symmetry(params,center,dp,kb,subID,expVer,iRun);
    end
catch ME
    %Screen('LoadNormalizedGammaTable', dp.wPtr, repmat(linspace(0,1, 256)',1,3));
    sca;
    rethrow(ME);
    ListenChar(1);
    ShowCursor;
end

if lab
    Screen('LoadNormalizedGammaTable', dp.wPtr, oldGamma);
end

sca;
ListenChar(1);
ShowCursor;

