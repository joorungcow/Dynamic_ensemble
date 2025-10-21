clear;
temp = [];
temp2 = [];
tempT = [];

sd1 = .4;
mu1 = 1.4;
g = 0.0001;
sd2 = .2;
mu2 = 1.6;
dotNum = 8; %num/2
minNum = .6;
maxNum = 2;
%%
t1 = GetSecs;

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
time = GetSecs - t1;

[sort(x1) sort(x2)]
[mean(x1) mean(x2)]
time

%%
figure; clf; hold on;
histogram(x1)
histogram(x2)
hold off;

%%
for i = 1:1000
    t1 = GetSecs;
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
    
    %[sort(x1) sort(x2)]
    %[mean(x1) mean(x2)]
    
    time = GetSecs - t1;
    temp = [temp; x1];
    temp2 = [temp2; x2];
    tempT = [tempT; time];

end

sum(tempT)/length(tempT)

temp = sort(temp);
temp2 = sort(temp2);
%%
figure; clf; hold on;
histogram(temp)
histogram(temp2)
hold off;
