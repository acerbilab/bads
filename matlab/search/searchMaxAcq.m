function z = searchMaxAcq(data,xi)
%SEARCHMAXACQ Wrapper of acquisition function for SEARCHMAX

if isstruct(xi) && isnumeric(data) % Arguments in wrong order
    temp = data;
    data = xi;
    xi = temp;
end

z = feval(data.acq{:},xi(:)',data.target,data.gpstruct,data.optimState,0);

end
