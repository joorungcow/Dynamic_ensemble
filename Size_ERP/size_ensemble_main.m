function size_ensemble_main
try
clear;
s = RandStream('mt19937ar','Seed','Shuffle');
RandStream.setGlobalStream(s);
KbName('UnifyKeyNames');

%% initial settings
dirname = 'EnsembleResults';
if ~isfolder(dirname)
    mkdir(dirname);
end

% get participants' info
p_info = struct();
p_info.randstream = s;
p_info.id = input('participant id: ', 's');
p_info.age = input('age: ');
p_info.biol_sex = input('biological sex: ', 's');
p_info.vision = input('normal (color) vision: ', 's');
p_info.hand = input('left-or a right-hander: ', 's');
p_info.expStart = datetime("now");

% file name
fn = fullfile(dirname, strcat(p_info.id, '_size_ensemble.mat'));

%% set parameters
setsize = [16; 24]; % edited on Nov.01.2023
weberf = [-0.244; -0.068; 0; 0.073; 0.323]; % Weber Fraction from Baek & Chong (2020) AP&P
nreps = 40;
tempmat = [sort(repmat(setsize, [length(weberf), 1])), repmat(weberf, [length(setsize), 1])];
expmat = [];
for iReps = 1: nreps
    expmat = [expmat; Shuffle(tempmat, 2)];
end

expseq = round(1:length(setsize)*length(weberf)*nreps);

expmat = [expseq', expmat];
stimmat = struct(); % save info of individual circles

% standard circles
nstdstim = 20; % number of circles
stdsize = 40; % initial value for the mean size of the items in the standard set (in radius, pixel) 

% timing
Fixdur = 0.5;
ensembledur = 0.5;
ITIDur = 1.5;

% key settings
onekey = [KbName('1'), KbName('1!')];
twokey = [KbName('2'), KbName('2@')];

%% EEG settings
daq_id = DaqDeviceIndex;
if ~isempty(daq_id), DaqDConfigPort(daq_id, 0, 0); end

% trigers
expstart = 101;
trialend = 99;
expend = 109;
taskstart = 77;
corrresp = 191;
incorrresp = 199;
neutresp = 195;
press1 = 1;
press2 = 2;

%% Open Screen
whichscreen = max(Screen('Screens'));
% Screen('Preference', 'Verbosity', 0);
Screen('Preference', 'SkipSyncTests', 0);
Screen('Preference','TextEncodingLocale','UTF-8');
Screen('Preference', 'TextRenderer', 0);
% Screen('Preference', 'VisualDebugLevel', 0);
%Screen('Preference', 'DefaultFontName', 'Arial');
Screen('Preference', 'DefaultFontSize', 28);
if strcmp(getenv('computername'), 'DESKTOP-OKJ1LHN') % 317-1
    [w, rect] = Screen('OpenWindow', whichscreen, [128 128 128]);
else
    [w, rect] = Screen('OpenWindow', whichscreen, [128 128 128], [0 0 1920 1080]);
end

Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
HideCursor;
[cx, cy] = RectCenter(rect);

%%% display information (change the values depending on the experimental settings)
dp.width = 60;
dp.dist = 98;
dp.ifi = Screen('GetFlipInterval', w);
dp.frameRate = round(1/dp.ifi); %Hz
dp.resolution = [rect(3), rect(4)];
dp.cx = cx;
dp.cy = cy;
dp.ppd = dp.resolution(1)/((2*atan(dp.width/2/dp.dist))*180/pi);
dp.wPtr = w;

bk = [0, 0, 0]; % black
wt = [255, 255, 255]; % white

[gridx, gridy] = make_grid2(dp.resolution);

% 1st save point
save(fn, 'p_info', 'dp', 'expmat', 'expseq', 'stimmat');

%% instruction
if ~isempty(daq_id), err = DaqDOut(daq_id, 0, expstart); end %%% send trigger
str = double('실험에 참여해주셔서 감사합니다.');
strBounds = Screen('TextBounds', w, str);
Screen('DrawText', w, str, cx-strBounds(3)/2, (cy-strBounds(4)/2)-400, bk);
 
str = double('본 실험에서는 원 집합 두 개가 순차적으로 제시됩니다.');
strBounds = Screen('TextBounds', w, str);
Screen('DrawText', w, str, cx-strBounds(3)/2, (cy-strBounds(4)/2)-300, bk);

str = double('제시되는 자극을 보고 평균 크기가 더 커보이는 집합이 어떤 것인지 반응해주십시오.');
strBounds = Screen('TextBounds', w, str);
Screen('DrawText', w, str, cx-strBounds(3)/2, (cy-strBounds(4)/2)-200, bk);

str = double('첫 번째로 제시된 원 집합의 평균 크기가 더 커보이면 1번,');
strBounds = Screen('TextBounds', w, str);
Screen('DrawText', w, str, dp.cx-strBounds(3)/2, cy-strBounds(4)/2-100, bk);

str = double('두 번째로 제시된 원 집합의 평균 크기가 더 커보이면 2번을 눌러주십시오.');
strBounds = Screen('TextBounds', w, str);
Screen('DrawText', w, str, cx-strBounds(3)/2, cy-strBounds(4)/2, bk);

str = double('스페이스 바를 누르면 실험이 시작됩니다.');
strBounds = Screen('TextBounds', w, str);
Screen('DrawText', w, str, cx-strBounds(3)/2, cy-strBounds(4)/2+200, bk);

Screen('Flip', w);

KbWait;
while(true)
    [~, ~,keyCode] = KbCheck;

    if any(keyCode(KbName('space')))
        do_practice = 0;
        break;
    elseif keyCode(KbName('p'))
        do_practice = 1;
        break;
    elseif any(keyCode(KbName('ESCAPE')))
        ShowCursor;
        sca;
        break;
    end
end

% %% practice routine (not ready yet)
% if do_practice
%     DrawFormattedText(w, 'Practice session will be added soon...', 'center', 'center', bk);
%     DrawFormattedText(w, 'Press the space bar.', 'center', cy+100, bk);
%     Screen('Flip', w);
% 
%     KbWait;
%     while(true)
%         [~, ~,keyCode] = KbCheck;
% 
%         if any(keyCode(KbName('space')))
%             break;
%         elseif any(keyCode(KbName('ESCAPE')))
%             ShowCursor;
%             sca;
%             break;
%         end
%     end
% end

%% draw
for i = expseq
   
    %%% label exp conds
    thisss = expmat(i, 2);
    thiswf = expmat(i, 3);
    
    %%% get stimulus locations
    % for test
    % grididx = randperm(length(grid), thisss);
    % tlocs = grid(:, grididx);
    grididx = Shuffle(find(~isnan(gridx)));
    grididx = grididx(1:thisss);
    tlocs = [gridx(grididx)'; gridy(grididx)']; 

    % for standard
    % grididx = randperm(length(grid), nstdstim);
    % slocs = grid(:, grididx);
    grididx = Shuffle(find(~isnan(gridx)));
    grididx = grididx(1:nstdstim);
    slocs = [gridx(grididx)'; gridy(grididx)']; 
   
    %%% ensemble stimuli (consider screen resolution and grid slot size)  
    [SindivR, TindivR, SmeanA_subj, TmeanA_subj] = set_stim_size(thisss, thiswf, nstdstim, stdsize);
    expmat(i, 4) = SmeanA_subj; % perceived mean area of the standard set
    expmat(i, 5) = TmeanA_subj; % perceived mean area of the test set

    stimmat(i).order = expmat(i, 1);
    stimmat(i).tlocs = tlocs;
    stimmat(i).slocs = slocs;
    stimmat(i).SindivR = SindivR;
    stimmat(i).TindivR = TindivR;

    %%% Make checkerboard textures
    [Simg, Timg] = check_texture(SindivR, TindivR, dp.ppd);

    if mod(i, 2) == 1
        for st = 1:nstdstim
            Stex(1, st) = Screen('MakeTexture', w, Simg(st).img1);
        end

        for tt = 1:thisss
            Ttex(1, tt) = Screen('MakeTexture', w, Timg(tt).img1);
        end        
    else
        for st = 1:nstdstim
            Stex(1, st) = Screen('MakeTexture', w, Simg(st).img2);
        end

        for tt = 1:thisss
            Ttex(1, tt) = Screen('MakeTexture', w, Timg(tt).img2);
        end 
    end

    %%% TRIGGER
    switch thiswf
        case -0.244
            TRIG = 1;
        case -0.068
            TRIG = 2;
        case 0
            TRIG = 3;
        case 0.073
            TRIG = 4;
        case 0.323
            TRIG = 5;
    end

    TRIG = (thisss/8 - 1) * 10 + TRIG; % 10^1 = set size, 10^0 = weber fraction
    stdset = 100 + TRIG; 

    %%% fixation
    draw_fixation(w, cx, cy, bk);
    t0 = Screen('Flip', w);
    while (GetSecs - t0) < Fixdur
        ;
    end
    % t1 = GetSecs;
    % t1d = t1 - t0

    %%% draw circles (test)
    for j = 1:thisss
        jitter = -10 + (10 + 10) .* rand();
        draw_fixation(w, cx, cy, bk);
        Screen('DrawTexture', w, Ttex(j), [], [tlocs(1, j) - TindivR(j) + jitter, tlocs(2, j) - TindivR(j) + jitter, tlocs(1, j) + TindivR(j) + jitter, tlocs(2, j) + TindivR(j) + jitter]);
        % Screen('FillRect', w, wt, [rect(3) - 100, rect(4) - 100, rect(3), rect(4)]);
        %Screen('FrameOval', w, bk, [tlocs(1, j) - TindivR(j) + jitter, tlocs(2, j) - TindivR(j) + jitter, tlocs(1, j) + TindivR(j) + jitter, tlocs(2, j) + TindivR(j) + jitter]); 
    end
    Screen('Flip', w);
    if ~isempty(daq_id), err = DaqDOut(daq_id, 0, TRIG); end %%% send trigger
    while (GetSecs - t0) < (Fixdur + ensembledur)
        ;
    end
    % t2 = GetSecs;
    % t2d = t2 - t1

    %%% fixation
    while GetSecs - t0 < (2 * Fixdur + ensembledur) - 0.5 * dp.ifi
        draw_fixation(w, cx, cy, bk);
        Screen('Flip', w);
    end
    while (GetSecs - t0) < (2 * Fixdur + ensembledur)
        ;
    end
    % t3 = GetSecs;
    % t3d = t3 - t2

    %%% draw circles (standard)
    for n = 1:nstdstim
        jitter = -10 + (10 + 10) .* rand();
        draw_fixation(w, cx, cy, bk);
        Screen('DrawTexture', w, Stex(n), [], [slocs(1, n) - SindivR(n) + jitter, slocs(2, n) - SindivR(n) + jitter, slocs(1, n) + SindivR(n) + jitter, slocs(2, n) + SindivR(n) + jitter]);
        % Screen('FillRect', w, wt, [rect(3) - 100, rect(4) - 100, rect(3), rect(4)]);
        % Screen('FrameOval', w, bk, [slocs(1, n) - SindivR(n) + jitter, slocs(2, n) - SindivR(n) + jitter, slocs(1, n) + SindivR(n) + jitter, slocs(2, n) + SindivR(n) + jitter]); 
    end
    Screen('Flip', w);
    if ~isempty(daq_id), err = DaqDOut(daq_id, 0, stdset); end %%% send trigger
    while (GetSecs - t0) < (2 * Fixdur +  2 * ensembledur)
        ;
    end
    % t4 = GetSecs;
    % t4d = t4 - t3
    
%     while GetSecs - t0 < (3 * Fixdur + 2 * ensembledur) - 0.5 * dp.ifi
%         draw_fixation(w, cx, cy, bk);
%         Screen('Flip', w);
%     end
%     while (GetSecs - t0) < (3 * Fixdur + 2 * ensembledur)
%         ;
%     end
    
    %%% task
    str = double('첫 번째로 제시된 자극의 평균크기가 더 커보이면 1번,');
    strBounds = Screen('TextBounds', w, str);
    Screen('DrawText', w, str, cx-strBounds(3)/2, (cy-strBounds(4)/2)-50, bk);

    str = double('두 번째로 제시된 자극의 평균크기가 더 커보이면 2번을 눌러주십시오.');
    strBounds = Screen('TextBounds', w, str);
    Screen('DrawText', w, str, cx-strBounds(3)/2, (cy-strBounds(4)/2)+50, bk);
    vbl = Screen('Flip', w);
    if ~isempty(daq_id), err = DaqDOut(daq_id, 0, taskstart); end %%% send trigger

    %%% get key response
    KbWait;
    while(true)
        [~, ~, keyCode] = KbCheck();
        RT = GetSecs - vbl;

        if any(keyCode(onekey))
            response = 1; % interval 1 (test) is bigger
            if ~isempty(daq_id), err = DaqDOut(daq_id, 0, press1); end %%% send trigger
            break;
        elseif any(keyCode(twokey))
            response = 2; % interval 2 (standard) is bigger
            if ~isempty(daq_id), err = DaqDOut(daq_id, 0, press2); end %%% send trigger
            break;
        end
    
        if any(keyCode(KbName('ESCAPE')))
            ShowCursor;
            sca;
            break;
        end
    end
    WaitSecs(.05); % Just wait a little for triggers
    %%% get accuracy
    if thiswf < 0 % test < standard
        if response == 1
            accuracy = 999;
            if ~isempty(daq_id), err = DaqDOut(daq_id, 0, incorrresp); end 
        else
            accuracy = 111;
            if ~isempty(daq_id), err = DaqDOut(daq_id, 0, corrresp); end 
        end
    elseif thiswf > 0 % standard < test
        if response == 1
            accuracy = 111;
            if ~isempty(daq_id), err = DaqDOut(daq_id, 0, corrresp); end
        else
            accuracy = 999;
            if ~isempty(daq_id), err = DaqDOut(daq_id, 0, incorrresp); end 
        end
    elseif thiswf == 0 % test and standard are same.
        accuracy = 555;
        if ~isempty(daq_id), err = DaqDOut(daq_id, 0, neutresp); end 
    end
    WaitSecs(.05); % Just wait a little for triggers
    if ~isempty(daq_id), err = DaqDOut(daq_id, 0, trialend); end 

    %%% update emat
    expmat(i, 6) = response;
    expmat(i, 7) = RT;
    expmat(i, 8) = accuracy;

    %%% save the result files (append)
    save(fn, 'p_info', 'dp', 'expmat', 'expseq', 'stimmat', '-append');

    %%% ITI
    draw_fixation(w, cx, cy, bk);
    vbl = Screen('Flip', w);
    while (GetSecs - vbl) < ITIDur
        ;
    end


    %%% break
    bt = length(weberf) * length(setsize) * 5;
     if mod(expmat(i, 1), bt) == 0 && expmat(i, 1) < length(expmat(:, 1))
        str = double('쉬는 시간입니다.');
        strBounds = Screen('TextBounds', w, str);
        Screen('DrawText', w, str, cx-strBounds(3)/2, cy-strBounds(4)/2, bk);
        Screen('Flip', w);
        WaitSecs(10);

        str = double('준비가 되면 스페이스 바를 눌러 다음 세션을 진행해주십시오.');
        strBounds = Screen('TextBounds', w, str);
        Screen('DrawText', w, str, cx-strBounds(3)/2, cy-strBounds(4)/2, bk);
        Screen('Flip', w);

        KbWait;
        while (true)
            [~, ~,keyCode] = KbCheck();

            if any(keyCode(KbName('space')))
                break;
            elseif any(keyCode(KbName('ESCAPE')))
                ShowCursor;
                sca;
                break;
            end
        end         
     end
end


%% byebye
str = double('실험에 참여해주셔서 감사합니다.');
strBounds = Screen('TextBounds', w, str);
Screen('DrawText', w, str, cx-strBounds(3)/2, cy-strBounds(4)/2, bk);
Screen('Flip', w);
if ~isempty(daq_id), err = DaqDOut(daq_id, 0, expend); end %%% send trigger
p_info.expEnd = datetime("now");
save(fn, 'p_info', 'dp', 'expmat', 'expseq', 'stimmat', '-append'); % for the last time

WaitSecs(3);
ShowCursor;
sca;

catch exception
    Screen('CloseAll');
    rethrow(exception);
end

end

function draw_fixation(w, cx, cy, col)
    Screen('DrawLine', w, col, cx-20, cy, cx+20, cy, 5); 
    Screen('DrawLine', w, col, cx, cy-20, cx, cy+20, 5);
end
