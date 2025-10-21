function [params, MeanSize] = makeStimuli(params, trial, MeanSize, dp, center, duraIndx)  
%grid setting (자극 제시 위치 계산)
    %params.virtualRectsize = (MeanSize + 2*0.15*MeanSize) * 1.4 * dp.ppd; % distance between circles, 2.4, 2.7
    params.virtualRectsize = 2.65 * dp.ppd;
    params.distance = 2 * dp.ppd; % origin: 3 degree, from fixation to circles
    
    for row = 1: 6 % origin: 1:4
        for column = 1: 6 % origin: 1:3 
            params.virtualCentersL{row,column} = [dp.cx - params.distance - (6-column)*params.virtualRectsize - params.virtualRectsize/2, ...
                dp.cy - 2*params.virtualRectsize - params.virtualRectsize/2 + (row-1)*params.virtualRectsize];
            params.virtualCentersR{row,column} = [dp.cx + params.distance + (column-1)*params.virtualRectsize + params.virtualRectsize/2, ...
                dp.cy - 2*params.virtualRectsize - params.virtualRectsize/2 + (row-1)*params.virtualRectsize];
        end
    end
    params.lCenters = params.virtualCentersL(:);
    params.rCenters = params.virtualCentersR(:);
    
%% each size
    dotNum = 8; % input num/2
    
    if params.sSizeSdBig(trial)
        mu1 = params.sSize;
        mu2 = MeanSize;
    else
        mu1 = MeanSize;
        mu2 = params.sSize;
    end

    g = 0.0001; %(허용된 오차범위)
    sd1 = .4;
    sd2 = .2;
    minNum = 0.6; %(크기의 최대, 최소 제한)
    maxNum = 2;
    
    isFound = 0;
    while ~isFound 
        x1 = normrnd(mu1,sd1,dotNum,1);
        if mean(x1) < mu1+g && mean(x1) > mu1-g && min(x1) > minNum && max(x1) < maxNum
            isFound = 1;
        end
    end
    isFound = 0;
    while ~isFound 
        x2 = normrnd(mu2,sd2,dotNum,1);
        if mean(x2) < mu2+g && mean(x2) > mu2-g && min(x2) > minNum && max(x2) < maxNum
            isFound = 1;
        end
    end
    
    if params.sSizeSdBig(trial)
        sSize = x1';
        tSize = x2';
    else
        sSize = x2';
        tSize = x1';
    end

    if params.debugFullDots
        fullSize = maxNum; 
        sSize = ones(1,18);
        tSize = ones(1,18);
        for i = 1:18
            sSize(i) = fullSize;
            tSize(i) = fullSize;
        end 
    end
    
    %dotSizes(trial,data,standard(1) or test(2))
    params.dotSizes(trial,1:8,1:2) = [sSize;tSize]'; %Col1 = sSizes, Col2 = tSizes


%% previous
%     k=[.1, .125, .15];
%     knum=Sample(1:length(k));
%     params.sSizeKnum(trial) = knum;
%     
% % disk number = 8*2 = 16
%     sSize = [params.sSize-2*k(knum)*params.sSize, params.sSize-1.7*k(knum)*params.sSize, params.sSize-1.4*k(knum)*params.sSize, params.sSize-1.1*k(knum)*params.sSize,...
%         params.sSize+1.1*k(knum)*params.sSize, params.sSize+1.4*k(knum)*params.sSize, params.sSize+1.7*k(knum)*params.sSize, params.sSize+2*k(knum)*params.sSize];
%     
%     knum=Sample(1:length(k));
%     params.tSizeKnum(trial) = knum;
% 
%     tSize = [MeanSize-2*k(knum)*MeanSize MeanSize-1.7*k(knum)*MeanSize MeanSize-1.4*k(knum)*MeanSize MeanSize-1.1*k(knum)*MeanSize ...
%         MeanSize+1.1*k(knum)*MeanSize MeanSize+1.4*k(knum)*MeanSize MeanSize+1.7*k(knum)*MeanSize MeanSize+2*k(knum)*MeanSize];

% % % full dots
% %     fullSize = 1.4+2*.15*1.4; %MeanSize Maximum+2*k(knum maximum)*MeanSize Maximum
% %     sSize = ones(1,18);
% %     tSize = ones(1,18);
% %     for i = 1:18
% %         sSize(i) = fullSize;
% %         tSize(i) = fullSize;
% %     end 
    
%%    
    % make matrix
    sMatrix = zeros(6,3);
    finishNum = length(sSize);
    startNum = 1;
    while startNum <= finishNum
        fRand = floor(rand*6)+1;
        sRand = floor(rand*3)+1;
        if sMatrix(fRand, sRand) == 0
            sMatrix(fRand, sRand) = sSize(startNum);
            startNum = startNum + 1;
        end    
    end
    sMatrix = [sMatrix, fliplr(sMatrix)];
 
    tMatrix = zeros(6,3);
    finishNum = length(tSize);
    startNum = 1;
    while startNum <= finishNum
        fRand = floor(rand*6)+1;
        sRand = floor(rand*3)+1;
        if tMatrix(fRand, sRand) == 0
            tMatrix(fRand, sRand) = tSize(startNum);
            startNum = startNum + 1;
        end    
    end
    tMatrix = [tMatrix, fliplr(tMatrix)];

    
%%
% Make Standard Array
    for ii = 1: length(sMatrix)*length(sMatrix)
        sRect{ii} = [0 0 sMatrix(ii)*dp.ppd sMatrix(ii)*dp.ppd];
    end
    
    % Make Test Array
    for ii = 1: length(tMatrix)*length(tMatrix)
        tRect{ii} = round([0 0 tMatrix(ii)*dp.ppd tMatrix(ii)*dp.ppd]);
    end
%%
    % Positioning
    tIndex1 = randperm(length(tMatrix)*length(tMatrix));
    sIndex1 = randperm(length(sMatrix)*length(sMatrix));
    tIndex2 = 1:length(tMatrix)*length(tMatrix);

    temp.jitter = 10;
    tempJR = zeros(6,6);
    tempJL = zeros(6,6);
    %random jitter
    for i = 1: 36
        rDir = ceil(rand(1)+.5);
        jitter = Sample(5:temp.jitter);
        if params.debugFullDots
            jitter = 10;
        end
        temp.rCenters{i} = params.rCenters{i} + (-1)^rDir*jitter;
        tempJR(i) = (-1)^rDir*jitter;
        rDir = ceil(rand(1)+.5);
        jitter = Sample(5:temp.jitter);
        temp.lCenters{i} = params.lCenters{i} + (-1)^rDir*jitter;
        tempJL(i) = (-1)^rDir*jitter;
    end
    %symmetry jitter, cf)tMatrix(7)==tMatrix(1,2)
    tempJR_S = reshape(tempJR(1:18),[6,3]);
    tempJR_S = fliplr(tempJR_S);
    tempJL_S = reshape(tempJL(1:18),[6,3]);
    tempJL_S = fliplr(tempJL_S);
    temp.sym_rCenters = temp.rCenters(1,1:18);
    temp.sym_lCenters = temp.lCenters(1,1:18);
    for i = 19:36
        temp.sym_rCenters{i}(1) = params.rCenters{i}(1) - tempJR_S(i-18);
        temp.sym_rCenters{i}(2) = params.rCenters{i}(2) + tempJR_S(i-18);
        temp.sym_lCenters{i}(1) = params.lCenters{i}(1) - tempJL_S(i-18);
        temp.sym_lCenters{i}(2) = params.lCenters{i}(2) + tempJL_S(i-18);
    end

    t0 = GetSecs;
    while GetSecs-t0 <= params.duration(duraIndx)
        if params.randLoc(trial) % test right visual field
            if params.randSym(trial) % test dots is symmetry dots
                for i = 1: length(tMatrix)*length(tMatrix)
                    Screen('FillOval', dp.wPtr, [1 1 1], CenterRectOnPoint(tRect{tIndex2(i)},temp.sym_rCenters{i}(1),temp.sym_rCenters{i}(2)));
                end
                for i = 1: length(sMatrix)*length(sMatrix)
                    Screen('FillOval', dp.wPtr, [1 1 1], CenterRectOnPoint(sRect{sIndex1(i)},temp.lCenters{i}(1),temp.lCenters{i}(2)));
                end
            else % test dots is random dots
                for i = 1: length(tMatrix)*length(tMatrix)
                    Screen('FillOval', dp.wPtr, [1 1 1], CenterRectOnPoint(tRect{tIndex1(i)},temp.rCenters{i}(1),temp.rCenters{i}(2)));
                end
                for i = 1: length(sMatrix)*length(sMatrix)
                    Screen('FillOval', dp.wPtr, [1 1 1], CenterRectOnPoint(sRect{sIndex1(i)},temp.lCenters{i}(1),temp.lCenters{i}(2)));
                end
            end
        else % test left
            if params.randSym(trial) %test dots is symmetry dots
                for i = 1: length(tMatrix)*length(tMatrix)
                    Screen('FillOval', dp.wPtr, [1 1 1], CenterRectOnPoint(tRect{tIndex2(i)},temp.sym_lCenters{i}(1),temp.sym_lCenters{i}(2)));
                end
                for i = 1: length(sMatrix)*length(sMatrix)
                    Screen('FillOval', dp.wPtr, [1 1 1], CenterRectOnPoint(sRect{sIndex1(i)},temp.rCenters{i}(1),temp.rCenters{i}(2)));
                end
            else % test dots is random dots
                for i = 1: length(tMatrix)*length(tMatrix)
                    Screen('FillOval', dp.wPtr, [1 1 1],CenterRectOnPoint(tRect{tIndex1(i)},temp.lCenters{i}(1),temp.lCenters{i}(2)));
                end
                for i = 1: length(sMatrix)*length(sMatrix)
                    Screen('FillOval', dp.wPtr, [1 1 1], CenterRectOnPoint(sRect{sIndex1(i)},temp.rCenters{i}(1),temp.rCenters{i}(2)));
                end
            end
        end

        make_fixation(dp,center,params.fixSize,params.fixColor)
        Screen('Flip', dp.wPtr);

        if params.debugFullDots
            KbQueueWait;
        end
        
    end
    make_fixation(dp,center,params.fixSize,params.fixColor)
    Screen('Flip', dp.wPtr);
end

