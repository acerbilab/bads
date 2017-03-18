function x = origunits(u,optimState)
%ORIGUNITS Convert grid-normalized vector(s) coordinate to original units

if isfield(optimState,'trinfo')
    if size(u,1) == 1
        x = transvars(u,'inv',optimState.trinfo);
    else
        x = zeros(size(u));
        for i = 1:size(u,1)
            x(i,:) = transvars(u(i,:),'inv',optimState.trinfo);
        end
    end
else
    if size(u,1) == 1
        x = (u.*optimState.scale) + optimState.x0;
    else
        x = bsxfun(@plus,bsxfun(@times,u,optimState.scale),optimState.x0);
    end
end