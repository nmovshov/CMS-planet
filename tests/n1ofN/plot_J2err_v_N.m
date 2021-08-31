function [fh, ah] = plot_J2err_v_N(T, kwargs)
arguments
    T (:,:) table
    kwargs.logx (1,1) logical = false
    kwargs.logy (1,1) logical = false
end
N = unique(T.N);
nx = unique(T.nx);
tind = find(T.N == max(N) & T.nx == max(nx));
Jinf = T.J2(tind);
for k=1:length(N)
    ind = T.N == N(k) & T.nx == max(nx);
    J2_err(k) = abs(T.J2(ind) - Jinf)/Jinf; %#ok<AGROW>
end
[fh, ah] = ngraf.get_canvas('proj');
plot(N, J2_err*1e6, '--+');
if kwargs.logx, ah.XScale = 'log'; end
if kwargs.logy, ah.YScale = 'log'; end
xlabel(sprintf('N (nx = %d)', max(nx)))
ylabel('$10^6(J_2-J_2^c)$')
hline(0.04)
end
