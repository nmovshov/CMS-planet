function compare_strats(degree, trueval)

close all
if nargin == 0 || isempty(degree), degree = 2; end
if nargin < 2
    data = load('truth.mat');
    if degree == 2, trueval = data.trueJ2; end
    if degree == 4, trueval = data.trueJ4; end
    if degree == 6, trueval = data.trueJ6; end
end

files = uigetfile('*.mat','multi','on');
if ~iscell(files), files = {files}; end

figure;

% first a straight plot
ah1 = axes;
ah1.Box = 'on';
hold(ah1, 'on')
set(ah1,'XTick',[128 256 512 1024 2048 4096])
set(ah1,'XScale', 'log','XMinorTick','on')
styles = {'--+', '--o', '--*','--x','--s','--d','--p','--h'};
for k=1:length(files)
    data = load(files{k});
    x = data.x;
    if degree == 2, y = data.y2; end
    if degree == 4, y = data.y4; end
    if degree == 6, y = data.y6; end
    plot(x, y, styles{k},...
        'displayname', data.name,...
        'MarkerSize',8,'LineWidth',2)
end
try
    hline(trueval, 'k--');
catch
end
ylabel(['J_',int2str(degree),'(N)']);
xlabel('N spheroids');
legend show

% then maybe a convergence plot
if exist('trueval', 'var')
    figure
    ah2 = axes;
    ah2.Box = 'on';
    hold(ah2, 'on')
    set(ah2,'XTick',[128 256 512 1024 2048 4096])
    set(ah2,'XScale', 'log','XMinorTick','on')
    set(ah2,'YScale', 'log','YMinorTick','on')
    styles = {'--+', '--o', '--*','--x','--s','--d','--p','--h'};
    for k=1:length(files)
        data = load(files{k});
        x = data.x;
        if degree == 2, y = data.y2; end
        if degree == 4, y = data.y4; end
        if degree == 6, y = data.y6; end
        plot(x, abs(y - trueval), styles{k},...
            'displayname', data.name,...
            'MarkerSize',8,'LineWidth',2)
    end
    ylabel(['\Delta J_',int2str(degree),'(N)']);
    xlabel('N spheroids');
    legend show
end
