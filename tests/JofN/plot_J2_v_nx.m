function [fh, ah] = plot_J2_v_nx(T, kwargs)
arguments
    T (:,:) table
    kwargs.logx (1,1) logical = false
    kwargs.logy (1,1) logical = false
end
N = unique(T.N);
[fh, ah] = ngraf.get_canvas('proj');
for k=1:length(N)
    ind = T.N == N(k);
    lh = plot(T.nx(ind), T.J2(ind), '--+');
    lh.DisplayName = sprintf('N=%d',N(k));
end
if kwargs.logx, ah.XScale = 'log'; end
if kwargs.logy, ah.YScale = 'log'; end
xlabel('nx')
ylabel('$J_2$')
legend('Location','ne')

end
