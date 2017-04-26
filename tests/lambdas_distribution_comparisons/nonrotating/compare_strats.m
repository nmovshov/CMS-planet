function compare_strats()

files = uigetfile('*.mat','multi','on');
if ~iscell(files), files = {files}; end
close all

% NMoI plot
figure
ah1 = axes;
ah1.Box = 'on';
hold(ah1, 'on')
set(ah1,'XTick',[128 256 512 1024 2048 4096])
set(ah1,'XScale', 'log','XMinorTick','on')
styles = {'--+', '--o', '--*','--x','--s','--d','--p','--h'};
for k=1:length(files)
    data = load(files{k});
    x = data.x;
    y = data.L;
    plot(x, y, styles{k},...
        'displayname', data.name,...
        'MarkerSize',8,'LineWidth',2)
end
ylabel('NMoI');
xlabel('N spheroids');
legend show

% betanorm plot
figure
ah1 = axes;
ah1.Box = 'on';
hold(ah1, 'on')
set(ah1,'XTick',[128 256 512 1024 2048 4096])
set(ah1,'XScale', 'log','XMinorTick','on')
styles = {'--+', '--o', '--*','--x','--s','--d','--p','--h'};
for k=1:length(files)
    data = load(files{k});
    x = data.x;
    y = data.beta;
    plot(x, y, styles{k},...
        'displayname', data.name,...
        'MarkerSize',8,'LineWidth',2)
end
ylabel('betanorm');
xlabel('N spheroids');
legend show
