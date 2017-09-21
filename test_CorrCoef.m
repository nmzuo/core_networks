
function [zout, p]=test_CorrCoef(z1, n1, z2, n2)
% compare the Corr Coef r1 and r2, from different sample size group: n1!=n2
% Ref: Hinkle DE, Wiersma W, Jurs SG. Applied Statistics for the Behavioral Sciences. 
%      5th ed. Boston: Houghton Mifflin; 2003.
% http://www.statisticssolutions.com/comparing-correlation-coefficients/

% Question
% http://stats.stackexchange.com/questions/12558/comparing-correlation-coefficients/249247#249247
% Citation
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4383486/
% Revised: Dec 2, 2016

%    z1=fisherz(r1);
%    z2=fisherz(r2);
    zout = (z1-z2)./sqrt(1./ (n1-3)+1./(n2-3));
    p=2*(1-normcdf(abs(zout))); % two tailed
end

function z = fisherz(r)
    z=0.5 * log((1+r)./(1-r));
end
