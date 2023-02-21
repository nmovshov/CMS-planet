function ah = plot_runtime_v_N(T, kwargs)
arguments
    T (:,:) table
    kwargs.logx (1,1) logical = false
    kwargs.logy (1,1) logical = false
end
nx = unique(T.nx);
[~, ah] = ngraf.get_canvas('proj');
for k=1:length(nx)
    ind = T.nx == nx(k);
    lh = plot(T.N(ind), T.runtime(ind), '--+');
    lh.DisplayName = sprintf('nx=%d',nx(k));
end
if kwargs.logx, ah.XScale = 'log'; end
if kwargs.logy, ah.YScale = 'log'; end
xlabel('N')
ylabel('run time [sec]')
legend('Location','nw')

end
