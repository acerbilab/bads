function u = gridunits(x,optimState)
%GRIDUNITS Convert vector(s) coordinates to grid-normalized units

if isfield(optimState,'trinfo')
    if size(x,1) == 1
        u = transvars(x,'dir',optimState.trinfo);
    else
        u = zeros(size(x));
        for i = 1:size(x,1)
            u(i,:) = transvars(x(i,:),'dir',optimState.trinfo);
        end
    end
else
    if size(x,1) == 1
        u = (x-optimState.x0)./optimState.scale;
    else
        u = bsxfun(@rdivide,bsxfun(@minus,x,optimState.x0),optimState.scale);
    end
end