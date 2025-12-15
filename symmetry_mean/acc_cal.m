clear

cTable = cbrewer('qual', 'Set1', 9);

d = dir('*_1_1_*.mat');
subjects = cell(1,length(d));
for i = 1: length(d)
    if ~strcmp(d(i).name(5),'_')
        subjects{i} = d(i).name(1:5);
    else
        subjects{i} = d(i).name(1:4);
    end
end

% Fectch raw data
for iSub = 1: length(subjects)
    d = dir(sprintf('%s*.mat',subjects{iSub}));
    d = circshift(d,9);
    for iBlock = 1: length(d)
        load(d(iBlock).name);
        %%% save the block order %%%
        if iBlock == 1
            blockOrder(:,iSub) = params.duration(3:end);
        end
        %%% raw data in each block %%%
        for iStair = 1: 4 % 1,2: random; 3,4:symmetric
            temp(:,iStair,iBlock,iSub) = result(iStair).threshold';
            ensbTempI(:,iStair,iBlock,iSub) = result(iStair).intensity';
            ensbTempR(:,iStair,iBlock,iSub) = result(iStair).response';
            ensbTempRT(:,iStair,iBlock,iSub) = result(iStair).responseTime';
        end
    end
end

% Do the real thing

temp2 = NaN * ones(5,4,2,length(subjects)); % 5 reps, 4 rand conditions (2 random, 2 symmetric)
                                              % 2 stim dur, subject
ensbI250 = NaN * ones(25,4,5,length(subjects));
ensbI500 = NaN * ones(25,4,5,length(subjects));
ensbR250 = NaN * ones(25,4,5,length(subjects));
ensbR500 = NaN * ones(25,4,5,length(subjects));
ensbRT250 = NaN * ones(25,4,5,length(subjects));
ensbRT500 = NaN * ones(25,4,5,length(subjects));


for iSub = 1: length(subjects)
    temp2(:,:,1,iSub) = squeeze(temp(25,:,blockOrder(:,iSub) == 0.25,iSub))';
    temp2(:,:,2,iSub) = squeeze(temp(25,:,blockOrder(:,iSub) == 0.50,iSub))';
    ensbI250(:,:,:,iSub) = ensbTempI(1:25,:,blockOrder(:,iSub) == 0.25, iSub);
    ensbI500(:,:,:,iSub) = ensbTempI(1:25,:,blockOrder(:,iSub) == 0.50, iSub);
    ensbR250(:,:,:,iSub) = ensbTempR(1:25,:,blockOrder(:,iSub) == 0.25, iSub);
    ensbR500(:,:,:,iSub) = ensbTempR(1:25,:,blockOrder(:,iSub) == 0.50, iSub);
    ensbRT250(:,:,:,iSub) = ensbTempRT(1:25,:,blockOrder(:,iSub) == 0.25, iSub);
    ensbRT500(:,:,:,iSub) = ensbTempRT(1:25,:,blockOrder(:,iSub) == 0.50, iSub);
end

results250 = NaN * ones(10,2,length(subjects));
results500 = NaN * ones(10,2,length(subjects));

for iSub = 1: length(subjects)
    results250(:,1,iSub) = [temp2(:,1,1,iSub); temp2(:,2,1,iSub)];
    results250(:,2,iSub) = [temp2(:,3,1,iSub); temp2(:,4,1,iSub)];
    results500(:,1,iSub) = [temp2(:,1,2,iSub); temp2(:,2,2,iSub)];
    results500(:,2,iSub) = [temp2(:,3,2,iSub); temp2(:,4,2,iSub)];
end
%% accuracy

solutions = ensbI250(1:25,:,:,:)>1.4;
sameSize = ensbI250(1:25,:,:,:)==1.4;
ensbAns250 = ensbR250;
ensbAns250(sameSize) = NaN;
corrAns = ensbAns250 == solutions;

for iSub = 1: length(subjects)
    totalCorr(1,iSub) = sum(sum(sum(corrAns(:,1:2,:,iSub)))); %random
    totalSame(1,iSub) = sum(sum(sum(sameSize(:,1:2,:,iSub))));
    totalCorr(2,iSub) = sum(sum(sum(corrAns(:,3:4,:,iSub)))); %symmetry
    totalSame(2,iSub) = sum(sum(sum(sameSize(:,3:4,:,iSub))));
end

indiAcc = totalCorr./(250-totalSame);

figure()
bar(indiAcc')

indiAccM = mean(indiAcc);
mean(indiAcc')

%% 500

solutions = ensbI500(1:25,:,:,:)>1.4;
sameSize = ensbI500(1:25,:,:,:)==1.4;
ensbAns500 = ensbR500;
ensbAns500(sameSize) = NaN;
corrAns = ensbAns500 == solutions;

for iSub = 1: length(subjects)
    totalCorr(1,iSub) = sum(sum(sum(corrAns(:,1:2,:,iSub)))); %random
    totalSame(1,iSub) = sum(sum(sum(sameSize(:,1:2,:,iSub))));
    totalCorr(2,iSub) = sum(sum(sum(corrAns(:,3:4,:,iSub)))); %symmetry
    totalSame(2,iSub) = sum(sum(sum(sameSize(:,3:4,:,iSub))));
end

indiAcc = totalCorr./(250-totalSame);

figure()
bar(indiAcc')

indiAccM = mean(indiAcc);
mean(indiAcc')

