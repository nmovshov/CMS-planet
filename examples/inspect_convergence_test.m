%% Examine convergence test output
clear
clc
close all

%% Load results
outfile = 'convtable.mat';
out = load(outfile);
T = out.T;
disp(T)

%% Let's extract some numbers
N = unique(T.N);
nx = unique(T.nx);

%% Let's plot runtime against nx
[fh, ah] = ngraf.get_canvas('proj');
for k=1:length(N)
    ind = T.N == N(k);
    lh = plot(T.nx(ind), T.runtime(ind), '--+');
    lh.DisplayName = sprintf('N=%d',N(k));
end
xlabel('nx')
ylabel('run time [sec]')
legend('Location','nw')

%% Let's look at J2 against nx for different N
[fh, ah] = ngraf.get_canvas('proj');
for k=1:length(N)
    ind = T.N == N(k);
    lh = plot(T.nx(ind), T.J2(ind), '--+');
    lh.DisplayName = sprintf('N=%d',N(k));
end
xlabel('nx')
ylabel('$J_2$')
legend('Location','ne')

%% And let's compare everything to highest resolution run
tind = find(T.N == max(N) & T.nx == max(T.nx));
Jinf = T.J2(tind);
[fh, ah] = ngraf.get_canvas('proj');
for k=1:length(N)
    ind = T.N == N(k);
    lh = plot(T.nx(ind), abs(T.J2(ind) - Jinf)*1e6, '--+');
    lh.DisplayName = sprintf('N=%d',N(k));
end
xlabel('nx')
ylabel('$10^6(J_2-J_2^c)$')
legend('Location','ne')
hline(0.04)
