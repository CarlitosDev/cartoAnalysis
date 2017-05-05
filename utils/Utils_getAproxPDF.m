function [vals_rate_norm,vals_rate,vals,zflag] = Utils_getAproxPDF(imaVOL,plot_it)

%Use:
%[vals_rate_norm,vals_rate,vals,zflag] = mypdf2(imaVOL,plot_it)
%

data = round(double(imaVOL(:)));  %Just integers!
ldata =length(data);
zflag=0; %obsolete

% PDF aproximation
old_min=0;
minval = min(data);
if (minval<0)
    old_min = minval;
    data = data - minval;
    minval = 0;
else
    minval = 0;
end

vals = minval:max(data);

vals_rate = zeros (size(vals));

for p=1:ldata
    if (data(p)>0),vals_rate (data(p)) = vals_rate(data(p))+1;end
end
vals_rate_norm = vals_rate./sum(vals_rate);

%Recover original values
if (old_min)
    vals = vals + old_min - zflag;
end

%Plot

if(plot_it)
figure,
plot(vals,vals_rate_norm);
title('Probability Density Function Approximation')
xlabel('Values');
ylabel('Probability');

end