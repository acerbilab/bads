function scatterplot(iter,ubest,fval,action,gpstruct,optimState,options)
%SCATTERPLOT Scatter plot of optimization iteration

if size(gpstruct.x, 1) == 0; return; end

figure(100);

nvars = size(gpstruct.x,2);
MeshSize = optimState.meshsize;

if options.Debug
    display(['GP hyper-parameters at iteration ' num2str(iter) ':']);
    gpstruct.hyp.mean
    (gpstruct.hyp.cov')
    (gpstruct.hyp.lik)
end

nrows = 2;
ncols = ceil((ceil(nvars/2) + 1)/2);

for iDim = 1:ceil(nvars/2)
    subplot(nrows,ncols,iDim);
    d1 = iDim*2 - 1;
    d2 = min(nvars,iDim*2);

    index = 1:optimState.Xmax;
    hold off;

    % Plot mesh size box
    s = MeshSize*ones(1,nvars);
    sstar = ubest;
    plot([sstar(d1) + s(d1),sstar(d1) + s(d1),sstar(d1) - s(d1),sstar(d1) - s(d1),sstar(d1) + s(d1)], ...
        [sstar(d2) - s(d2),sstar(d2) + s(d2),sstar(d2) + s(d2),sstar(d2) - s(d2),sstar(d2) - s(d2)], ...
        '-','LineWidth',1,'Color',[0.25 0.5 0]); hold on;

    % Plot all points
    idx = isfinite(optimState.Y(index(1:end-1)));            
    plot(optimState.U(idx,d1),optimState.U(idx,d2),'.','MarkerSize',5,'MarkerFaceColor',[0 0 0]);
    idx = ~isfinite(optimState.Y(index(1:end-1)));
    plot(optimState.U(idx,d1),optimState.U(idx,d2),'*','MarkerSize',5,'MarkerFaceColor',[0 0 0]);

    % Plot training set
    if ~isfield(gpstruct,'erry') || isempty(gpstruct.erry); erry = false(size(gpstruct.y)); else erry = gpstruct.erry; end
    idx = ~erry;
    plot(gpstruct.x(idx,d1),gpstruct.x(idx,d2),'o','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','none');
    idx = erry;
    plot(gpstruct.x(idx,d1),gpstruct.x(idx,d2),'*','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5]);

    miny = min(optimState.Y(index));
    maxy = max(optimState.Y(index));
    deg = (optimState.Y(index) - miny)./(maxy - miny);
    %for iX = 1:optimState.Xmax
    %    col = bsxfun(@plus, deg*0.8, 0.1*[1 1 1]);
    %    plot(optimState.X(iX,d1),optimState.X(iX,d2),'o','MarkerFaceColor',col(iX,:)); hold on;      
    %end

    % Plot true minimum
    if ~isempty(options.TrueMinX)
        plot(options.TrueMinX(d1),options.TrueMinX(d2),'o','MarkerFaceColor',[0.25 0 1], 'MarkerEdgeColor', 'none');
    end

    % Plot last evaluated point
    plot(optimState.U(index(end),d1),optimState.U(index(end),d2),'o','MarkerFaceColor',[1 0.4 0.4]);

    % Plot incumbent
    plot(sstar(d1),sstar(d2),'d','MarkerFaceColor',[0.25 1 0], 'MarkerEdgeColor', 'none');

    box off;
    set(gca,'TickDir','out');
    
    % Zoom level
    s = 4*MeshSize*ones(1,nvars);
    abox = [sstar(d1),sstar(d1),sstar(d2),sstar(d2)] + [-s(d1),s(d1),-s(d2),s(d2)];    
    axis(abox);
end

subplot(nrows,ncols,nrows*ncols);
cla(gca);
box off;
axis off;
text(0.2,0.9,['Iteration: ' num2str(iter)]);
text(0.2,0.8,['f-count: ' num2str(optimState.funccount)]);
text(0.2,0.7,['f(x): ' num2str(fval,'%.6g')]);
text(0.2,0.6,['MeshScale: ' num2str(MeshSize,'%.6g')]);
text(0.2,0.5,['Method: ']);
text(0.2,0.4,['Action: ' action]);

set(gcf,'Color','w');
drawnow;
    
end
