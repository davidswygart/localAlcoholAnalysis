% navigate to data folder in MATLAB "Current Folder Browser"
paths = dataAndFigDirectoryPaths(pwd());
load([paths.data, 'goodClusters.mat'])
isFemale = contains(goodClusters.sex, 'female');
%%%%%%%%%%%%%%%%%%% Units %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 2: Injection
isControl = contains(goodClusters.group,'control') | contains(goodClusters.group,'drink');
isInject = contains(goodClusters.group,'inject');

x = categorical({'control','inject'});
n = [sum(isControl & isFemale), sum(isControl & ~isFemale);
    sum(isInject & isFemale), sum(isInject & ~isFemale)];


figure(1)
bar(x,n, 'stacked')
legend('F','M')
title('Injection: units')
ylabel('# units')

%% Fig 3: Sipper
isControl = contains(goodClusters.group,'control');
isDrink = contains(goodClusters.group,'drink');
isInject = contains(goodClusters.group,'inject');

x = categorical({'control','drink','inject'});
n = [sum(isControl & isFemale), sum(isControl & ~isFemale);
    sum(isDrink & isFemale), sum(isDrink & ~isFemale);
    sum(isInject & isFemale), sum(isInject & ~isFemale)];

figure(2)
bar(x,n, 'stacked')
title('Sipper: units')
legend('F','M')
ylabel('# units')

%%
%%%%%%%%%%%%%%%%%%% Unique animals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, firstUnique, ~] = unique(goodClusters.matName);
first = goodClusters(firstUnique,:);
isFemale = contains(first.sex, 'female');
%% Fig 2: Injection
isControl = contains(first.group,'control') | contains(first.group,'drink');
isInject = contains(first.group,'inject');

x = categorical({'control','inject'});
n = [sum(isControl & isFemale), sum(isControl & ~isFemale);
    sum(isInject & isFemale), sum(isInject & ~isFemale)];


figure(3)
bar(x,n, 'stacked')
legend('F','M')
title('Injection: Unique animals')
ylabel('# Unique animals')

%% Fig 3: Sipper
isControl = contains(first.group,'control');
isDrink = contains(first.group,'drink');
isInject = contains(first.group,'inject');

x = categorical({'control','drink','inject'});
n = [sum(isControl & isFemale), sum(isControl & ~isFemale);
    sum(isDrink & isFemale), sum(isDrink & ~isFemale);
    sum(isInject & isFemale), sum(isInject & ~isFemale)];

figure(4)
bar(x,n, 'stacked')
title('Sipper: Unique animals')
legend('F','M')
ylabel('# Unique animals')