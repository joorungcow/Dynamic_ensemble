function make_fixation(dp,fixSize,color,lineWidth)

if ~exist('fixSize','var')
    fixSize = 20;
end

if ~exist('lineWidth','var')
    lineWidth = 2;
end

if ~exist('color','var')
    color = [1 1 1];
end

Screen('DrawLine', dp.wPtr, color, dp.cx-fixSize/2, dp.cy, dp.cx+fixSize/2, dp.cy, lineWidth)
Screen('DrawLine', dp.wPtr, color, dp.cx, dp.cy-fixSize/2, dp.cx, dp.cy+fixSize/2, lineWidth)
