clear variables
Screen('Preference', 'SkipSyncTests', 1);

screens = Screen('Screens');
dp.screenNum = max(screens);
dp.bkColor = [.5 .5 .5];
dp.width = 60;
dp.dist = 55;

dp = OpenWindow(dp);
ListenChar(0);
HideCursor; 
kb = init_keyboard;

try
    sMatrix = 5;
    
    sRect = [0 0 sMatrix*dp.ppd sMatrix*dp.ppd];
    byun = sMatrix*dp.ppd;
    [x, y] = meshgrid(-byun /2:byun/2, -byun/2:byun/2);
    imStandard = .5*ones(size(x,1),size(y,1));
    id = x.^2 + y.^2 < (-byun/2)^2;
    imStandard(id) = 1; 
    sTexture = Screen('MakeTexture', dp.wPtr, imStandard);

    Screen('DrawTextures', dp.wPtr, sTexture, sRect, CenterRectOnPoint(sRect, dp.cx, dp.cy), 0);
    Screen('Flip', dp.wPtr);

    isResponse = 0;
    while ~isResponse
        [keyIsDown, secs, keyCode] = KbCheck(-1);
        if keyCode(kb.spaceKey)
            isResponse = 1;
        end

    end

catch ME
    sca;
    rethrow(ME);
    ListenChar(1);
    ShowCursor;
end

sca;
ListenChar(1);
ShowCursor;
