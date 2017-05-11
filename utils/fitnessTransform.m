function [z,nonlinf,nonlinmu,deltay] = fitnessTransform(y)
%FITNESSTRANSFORM Nonlinear transformation of objective function

% y-value order
[~,yord] = sort(y,'ascend');
%y = log(numel(yord)/2+1)-log(yord);
%y = -numel(yord)*(gpstruct.y./sum(gpstruct.y));

%hold off;
%histogram(gpstruct.y,20); drawnow;
%pp = [prctile1(gpstruct.y,0),prctile1(gpstruct.y,10),prctile1(gpstruct.y,25),prctile1(gpstruct.y,50),prctile1(gpstruct.y,75),prctile1(gpstruct.y,90),prctile1(gpstruct.y,100)]
%diff(pp)

% Compute percentiles
iq10 = max(1,round(0.10*numel(yord)));
iq25 = max(1,round(0.25*numel(yord)));
iq50 = min(numel(yord),max(1,round(0.5*numel(yord))));
iq75 = min(numel(yord),round(0.75*numel(yord)));

z = y;

if 0
    nonlinf = @(x_,mu_,delta_) 1*sqrt((x_ - mu_)/1 +1) + mu_ - 1;
elseif 0    
    nonlinf = @(x_,mu_,delta_) mu_ + delta_*(2./(1+exp(-2*(x_ - mu_)/delta_)) - 1);
elseif 1
    nonlinf = @(x_,mu_,delta_) x_.*(x_ < mu_) + 2*delta_*(sqrt((x_-mu_)/delta_+1) + mu_/2/delta_ -1).*(x_ >= mu_);
else
    nonlinf = @(x_,mu_,delta_) x_.*(x_ < mu_) + (delta_*(log((x_-mu_)/delta_ +1)+mu_/delta_)).*(x_ >= mu_);
end

if 0
    idx = yord(iq75:numel(yord));
    gamma = 100;
    nonlinmu = y(yord(iq75));
    deltay = gamma*max(1,y(yord(iq75)) - y(yord(iq25)));
    z(idx) = nonlinf(y(idx),nonlinmu,deltay);
else
    idx = yord(1:numel(yord));
    gamma = 40;
    nonlinmu = y(yord(1));
    deltay = gamma*max(1,y(yord(iq50)) - y(yord(iq25)));
    z(idx) = nonlinf(y(idx),nonlinmu,deltay);
end

% Plot trasnformed vs non-transformed values
if 0
    figure(999)
    hold off;                
    scatter(y,y); hold on; scatter(y,z);
    set(gcf,'Color','w'); set(gca,'TickDir','out');
    drawnow;
end
% y = (yord-1)./(numel(yord)-1) + 5*(yord-1).^2./(numel(yord)-1);


% y(yord(ceil(ntrain/2):end)) = y(yord(ceil(ntrain/2)));

end