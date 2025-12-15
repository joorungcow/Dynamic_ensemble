clear

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
        end
    end
end

% Do the real thing

temp2 = NaN * ones(5,4,2,length(subjects)); % 5 reps, 4 rand conditions (2 random, 2 symmetric)
                                              % 2 stim dur, subject
for iSub = 1: length(subjects)
    temp2(:,:,1,iSub) = squeeze(temp(25,:,blockOrder(:,iSub) == 0.25,iSub))';
    temp2(:,:,2,iSub) = squeeze(temp(25,:,blockOrder(:,iSub) == 0.50,iSub))';
end

results250 = NaN * ones(10,2,length(subjects));
results500 = NaN * ones(10,2,length(subjects));

for iSub = 1: length(subjects)
    results250(:,1,iSub) = [temp2(:,1,1,iSub); temp2(:,2,1,iSub)];
    results250(:,2,iSub) = [temp2(:,3,1,iSub); temp2(:,4,1,iSub)];
    results500(:,1,iSub) = [temp2(:,1,2,iSub); temp2(:,2,2,iSub)];
    results500(:,2,iSub) = [temp2(:,3,2,iSub); temp2(:,4,2,iSub)];
end

%% indivisual data plots
withPlots = 0;

if withPlots
    for i = 1 : length(subjects)
        figure(200+i); clf; hold on
        set(gca,'YLim',[1.2 1.6])
        scatter(1, results250(:,1,i),"blue") %random
        scatter(2, results250(:,2,i),"red") 
        scatter(3, results500(:,1,i),"blue")
        scatter(4, results500(:,2,i),"red") 
    end
end

%% boxplots

if withPlots
    for iSub = 1:length(subjects)
        figure(100+iSub); clf; hold on
        subplot(2,2,1); hold on;
        boxplot(results250(:,1,iSub),'whisker',w95)
        subplot(2,2,2); hold on;
        boxplot(results250(:,2,iSub),'whisker',w95)
        subplot(2,2,3); hold on;
        boxplot(results500(:,1,iSub),'whisker',w95)
        subplot(2,2,4); hold on;
        boxplot(results500(:,2,iSub),'whisker',w95)
    end
end

%% outlier (data point)

out250 = results250;
out500 = results500;

q3 = norminv(.75);
q95=norminv(.95);
w95=(q95-q3)/(2*q3);

for iSub = 1:length(subjects)
    for j = 1:2
        q25 = prctile(out250(:,j,iSub),25);
        q75 = prctile(out250(:,j,iSub),75);
        out25 = q25 - (w95*(q75-q25));
        out75 = q75 + (w95*(q75-q25));
        id25 = out25>out250(:,j,iSub);
        id75 = out75<out250(:,j,iSub);
        out250(id25,j,iSub) = NaN;
        out250(id75,j,iSub) = NaN;
    end
end

for iSub = 1:length(subjects)
    for j = 1:2
        q25 = prctile(out500(:,j,iSub),25);
        q75 = prctile(out500(:,j,iSub),75);
        out25 = q25 - (w95*(q75-q25));
        out75 = q75 + (w95*(q75-q25));
        id25 = out25>out500(:,j,iSub);
        id75 = out75<out500(:,j,iSub);
        out500(id25,j,iSub) = NaN;
        out500(id75,j,iSub) = NaN;
    end
end

%% mean

mean250(:,1) = mean(results250(:,1,:));
mean250(:,2) = mean(results250(:,2,:));
mean500(:,1) = mean(results500(:,1,:));
mean500(:,2) = mean(results500(:,2,:));

graphMean = [mean(mean250(:,1)) mean(mean250(:,2)); mean(mean500(:,1)) mean(mean500(:,2))];

%SEM
sem_rand_250 = std(mean250(:,1))./sqrt(length(mean250(:,1)));
sem_symm_250 = std(mean250(:,2))./sqrt(length(mean250(:,2)));
sem_rand_500 = std(mean500(:,1))./sqrt(length(mean500(:,1)));
sem_symm_500 = std(mean500(:,2))./sqrt(length(mean500(:,2)));

graphSem = [sem_rand_250 sem_symm_250; sem_rand_500 sem_symm_500];

%% mean (outlier)

mean250Ot(:,1) = nanmean(out250(:,1,:));
mean250Ot(:,2) = nanmean(out250(:,2,:));
mean500Ot(:,1) = nanmean(out500(:,1,:));
mean500Ot(:,2) = nanmean(out500(:,2,:));

graphMean2 = [mean(mean250Ot(:,1)) mean(mean250Ot(:,2)); mean(mean500Ot(:,1)) mean(mean500Ot(:,2))];

%SEM
sem_rand_250_2 = std(mean250Ot(:,1))./sqrt(length(mean250Ot(:,1)));
sem_symm_250_2 = std(mean250Ot(:,2))./sqrt(length(mean250Ot(:,2)));
sem_rand_500_2 = std(mean500Ot(:,1))./sqrt(length(mean500Ot(:,1)));
sem_symm_500_2 = std(mean500Ot(:,2))./sqrt(length(mean500Ot(:,2)));

graphSem2 = [sem_rand_250_2 sem_symm_250_2 ; sem_rand_500_2 sem_symm_500_2];

%% graph
figure(1); clf; hold on
x = categorical({'250ms','500ms'});
x = reordercats(x,{'250ms','500ms'});
y = graphMean;
b = bar(x,y);
set(gca,'YLim',[0 1.6])

[ngroups,nbars] = size(y);

xx = nan(nbars, ngroups);
for i = 1:nbars
    xx(i,:) = b(i).XEndPoints;
end

er = errorbar(xx', y, graphSem, 'k', 'linestyle', 'none');

title(sprintf('PSE of mean size judgments'))
xtips1 = b(1).XEndPoints;
ytips1 = er(1).YData;
ytips1 = ytips1 + 0.01;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
xtips2 = b(2).XEndPoints;
ytips2 = er(2).YData;
ytips2 = ytips2 + 0.01;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
xlabel('Duration');
ylabel('PSE');
legend('random','symmetric');
hold off

%% t-test
%[h,p,ci,stats]

disp("250ms")
[h,p] = ttest(mean250(:,1), mean250(:,2))
disp("500ms")
[h,p] = ttest(mean500(:,1), mean500(:,2))
disp("all")
[h,p] = ttest([mean250(:,1); mean500(:,1)], [mean250(:,2); mean500(:,2)])


%% graph
figure(2); clf; hold on
x = categorical({'250ms','500ms'});
x = reordercats(x,{'250ms','500ms'});
y = graphMean2;
b = bar(x,y);
set(gca,'YLim',[0 1.6])

[ngroups,nbars] = size(y);

xx = nan(nbars, ngroups);
for i = 1:nbars
    xx(i,:) = b(i).XEndPoints;
end

er = errorbar(xx', y, graphSem2, 'k', 'linestyle', 'none');

title(sprintf('PSE of mean size judgments -outlier'))
xtips1 = b(1).XEndPoints;
ytips1 = er(1).YData;
ytips1 = ytips1 + 0.01;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
xtips2 = b(2).XEndPoints;
ytips2 = er(2).YData;
ytips2 = ytips2 + 0.01;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
xlabel('Duration');
ylabel('PSE');
legend('random','symmetric');
hold off

%% t-test

disp("250ms")
[h,p] = ttest(mean250Ot(:,1), mean250Ot(:,2))
disp("500ms")
[h,p] = ttest(mean500Ot(:,1), mean500Ot(:,2))
disp("all")
[h,p] = ttest([mean250Ot(:,1); mean500Ot(:,1)], [mean250Ot(:,2); mean500Ot(:,2)])

%% data plots
figure(3); clf; hold on
x1 = [1 2.6];
x2 = [1.6 3.2];
y1 = graphMean(:,1)';
y2 = graphMean(:,2)';
sz = 200;
c1 = [0 0.4470 0.7410];
c2 = [0.8500 0.3250 0.0980];

set(gca,'XLim',[0.5 3.7],'FontSize',12)
xticks([1.3 2.9])
xticklabels({'250ms','500ms'})
set(gca,'YLim',[1.28 1.5])
hline = refline(0,1.4);
hline.Color = [0.5 0.5 0.5];

%indivisual data
yI = [mean250(:,1), mean250(:,2), mean500(:,1), mean500(:,2)];
cI1 = [0.5 0.6470 0.9410];
sI1 = scatter(x1,yI(:,1:2:3),[],cI1);
cI2 = [0.9500 0.5250 0.4980];
sI2 = scatter(x2,yI(:,2:2:4),[],cI2);

%linking
xl = [1 1.6];
yl = [mean250(:,1) mean250(:,2)];
line(xl,yl,'Color',[0.8 0.8 0.8]);

xl = [2.6 3.2];
yl = [mean500(:,1) mean500(:,2)];
line(xl,yl,'Color',[0.8 0.8 0.8]);

%mean
%er1 = errorbar(x1, y1, CILeng(:,1)', 'k', 'linestyle', 'none'); %95% CI, t value
%er2 = errorbar(x2, y2, CILeng(:,2)', 'k', 'linestyle', 'none');

s1 = scatter(x1,y1,sz,c1,"filled","square");
s2 = scatter(x2,y2,sz,c2,"filled","square");

er1 = errorbar(x1, y1, graphSem(:,1)', 'k', 'linestyle', 'none'); %SEM
er2 = errorbar(x2, y2, graphSem(:,2)', 'k', 'linestyle', 'none');

lgd = legend([s1,s2],{'Random','Symmetric'});
lgd.FontSize = 12;

t = title(sprintf('PSE of mean size judgments'));
t.FontSize = 17;
xt = xlabel('Duration');
yt = ylabel('PSE');
xt.FontSize = 17;
yt.FontSize = 17;

%% data plots -outlier
figure(4); clf; hold on
x1 = [1 2.6];
x2 = [1.6 3.2];
y1 = graphMean2(:,1)';
y2 = graphMean2(:,2)';
sz = 200;
c1 = [0 0.4470 0.7410];
c2 = [0.8500 0.3250 0.0980];

set(gca,'XLim',[0.5 3.7],'FontSize',12)
xticks([1.3 2.9])
xticklabels({'250ms','500ms'})
set(gca,'YLim',[1.28 1.5])
hline = refline(0,1.4);
hline.Color = [0.5 0.5 0.5];

%indivisual data
yI = [mean250Ot(:,1), mean250Ot(:,2), mean500Ot(:,1), mean500Ot(:,2)];
cI1 = [0.5 0.6470 0.9410];
sI1 = scatter(x1,yI(:,1:2:3),[],cI1);
cI2 = [0.9500 0.5250 0.4980];
sI2 = scatter(x2,yI(:,2:2:4),[],cI2);

%linking
xl = [1 1.6];
yl = [mean250Ot(:,1) mean250Ot(:,2)];
line(xl,yl,'Color',[0.8 0.8 0.8]);

xl = [2.6 3.2];
yl = [mean500Ot(:,1) mean500Ot(:,2)];
line(xl,yl,'Color',[0.8 0.8 0.8]);

s1 = scatter(x1,y1,sz,c1,"filled","square");
s2 = scatter(x2,y2,sz,c2,"filled","square");

er1 = errorbar(x1, y1, graphSem2(:,1)', 'k', 'linestyle', 'none'); %SEM
er2 = errorbar(x2, y2, graphSem2(:,2)', 'k', 'linestyle', 'none');

lgd = legend([s1,s2],{'random','symmetric'});
lgd.FontSize = 12;

t = title(sprintf('PSE of mean size judgments -outlier'));
t.FontSize = 17;
xt = xlabel('Duration');
yt = ylabel('PSE');
xt.FontSize = 17;
yt.FontSize = 17;
%%
EnsbBias250 = (1-(mean250(:,2)./mean250(:,1)))*100;
EnsbBias500 = (1-(mean500(:,2)./mean500(:,1)))*100;

EnsbBias250Ot = (1-(mean250Ot(:,2)./mean250Ot(:,1)))*100;
EnsbBias500Ot = (1-(mean500Ot(:,2)./mean500Ot(:,1)))*100;


sEnsbBias250 = std(EnsbBias250)/sqrt(length(EnsbBias250));
sEnsbBias500 = std(EnsbBias500)/sqrt(length(EnsbBias500));
sEnsbBias250Ot = std(EnsbBias250Ot)/sqrt(length(EnsbBias250Ot));
sEnsbBias500Ot = std(EnsbBias500Ot)/sqrt(length(EnsbBias500Ot));

%%
figure(5); clf; hold on;
x = categorical({'250ms','250ms -outlier', '500ms', '500ms -outlier'});
x = reordercats(x, {'250ms','250ms -outlier', '500ms', '500ms -outlier'});
y = [mean(EnsbBias250),mean(EnsbBias250Ot), mean(EnsbBias500), mean(EnsbBias500Ot)];
b = bar(x,y);
t = title(sprintf('Mean of mean size bias (%%)'));
set(gca, 'YLim', [0 1.8])

sem = [sEnsbBias250, sEnsbBias250Ot, sEnsbBias500, sEnsbBias500Ot];
er = errorbar(b.XEndPoints, y, sem, 'k', 'linestyle', 'none');

for i = 1:length(x)
    xtips = b.XEndPoints(i);
    ytips = er.YData(i);
    ytips = ytips + 0.6;
    labels = string(b.YData(i));
    text(xtips,ytips,labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
end

%%
mdl = fitlm(EnsbBias250,EnsbBias500)
figure(6); clf; hold on
p = plot(mdl);
t = title(sprintf('Comparison of the duration condition'));
t.FontSize = 17;
xt = xlabel('Bias in Perceived Mean Size (250ms, %)');
yt = ylabel('Bias in Perceived Mean Size (500ms, %)');
xt.FontSize = 15;
yt.FontSize = 15;

p(1).Marker = "o";
p(1).LineWidth = 1.3;
p(2).LineWidth = 2;
p(3).LineWidth = 2;
p(4).LineWidth = 2;

%%
mdl = fitlm(EnsbBias250Ot,EnsbBias500Ot)
figure(7); clf; hold on
p = plot(mdl);
t = title(sprintf('Comparison of the duration condition - outlier'));
t.FontSize = 17;
xt = xlabel('Bias in Perceived Mean Size (250ms, %)');
yt = ylabel('Bias in Perceived Mean Size (500ms, %)');
xt.FontSize = 15;
yt.FontSize = 15;

p(1).Marker = "o";
p(1).LineWidth = 1.3;
p(2).LineWidth = 2;
p(3).LineWidth = 2;
p(4).LineWidth = 2;

%% 2-sample t-test

calthis = 0;

if calthis
    figure()
    h1 = histogram(mean250(:,1));
    hold on
    h2 = histogram(mean250(:,2));
    
    h1.Normalization = 'probability';
    h1.BinWidth = 0.01;
    h2.Normalization = 'probability';
    h2.BinWidth = 0.01;
    
    [h,p,ci,stats] = vartest2(mean250(:,1),mean250(:,2))
    [h,p,ci,stats] = vartest2(mean250(:,1),mean250(:,2),'Tail','left') %symmetric var > random var
    
    
    figure()
    h1 = histogram(mean500(:,1));
    hold on
    h2 = histogram(mean500(:,2));
    
    h1.Normalization = 'probability';
    h1.BinWidth = 0.01;
    h2.Normalization = 'probability';
    h2.BinWidth = 0.01;
    
    [h,p,ci,stats] = vartest2(mean500(:,1),mean500(:,2))
    [h,p,ci,stats] = vartest2(mean500(:,1),mean500(:,2),'Tail','left') %symmetric var > random var
end
%% make csv

withcsv = 0;

T = readtable('subans.xlsx');
tempT = sortrows(T,2);
sortT = [tempT(6:end,:);tempT(1:5,:)];
sawPattern = logical(sortT.sawPattern);
sawPattern1 = [sawPattern;sawPattern;sawPattern;sawPattern];

if withcsv
    data1 = [mean250(:,1); mean250(:,2); mean500(:,1); mean500(:,2)];
    temp = ones(length(mean250(:,1)),1);
    sub = [1:length(mean250(:,1))]';
    pattern1 = [temp; temp*2; temp; temp*2];
    duration1 = [temp; temp; temp*2; temp*2];
    sub1 = [sub; sub; sub; sub];
    T = table(data1, sub1, pattern1, duration1, sawPattern1);
    writetable(T, 'resultcsv.csv')

    data2 = [mean250Ot(:,1); mean250Ot(:,2); mean500Ot(:,1); mean500Ot(:,2)];
    temp2 = ones(length(mean250Ot(:,1)),1);
    sub0 = [1:length(mean250Ot(:,1))]';
    pattern2 = [temp2; temp2*2; temp2; temp2*2];
    duration2 = [temp2; temp2; temp2*2; temp2*2];
    sub2 = [sub0; sub0; sub0; sub0];
    T = table(data2, sub2, pattern2, duration2, sawPattern1);
    writetable(T, 'resultcsvOT.csv')
end

