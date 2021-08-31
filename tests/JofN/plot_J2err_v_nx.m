function [fh, ah] = plot_J2err_v_nx(T, kwargs)
arguments
    T (:,:) table
    kwargs.logx (1,1) logical = false
    kwargs.logy (1,1) logical = false
end
N = unique(T.N);
nx = unique(T.nx);
tind = find(T.N == max(N) & T.nx == max(nx));
Jinf = T.J2(tind);
for k=1:length(nx)
    ind = T.nx == nx(k) & T.N == max(N);
    J2_err(k) = abs(T.J2(ind) - Jinf)/Jinf; %#ok<AGROW>
end
[fh, ah] = ngraf.get_canvas('proj');
plot(nx, J2_err*1e6, '--+');
if kwargs.logx, ah.XScale = 'log'; end
if kwargs.logy, ah.YScale = 'log'; end
xlabel(sprintf('nx (N = %d)', max(N)))
ylabel('$10^6(J_2-J_2^c)$')
hline(0.04)
end
