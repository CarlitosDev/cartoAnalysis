%% Find the trend of a simple signal

% Sin-based waveform with fd = 20Hz

load(fullfile('testSplines', 'dataForSplines.mat'))
idx    = 2;
yTrend = dataForICA(idx, :);

%% Find the trend using a spline - JL

siz = numel(yTrend);
Fe  = 1; %0.1 ? %10 ??
L_w = 200;

t = (0:siz-1)'/Fe;
n_t = 1:L_w/2:siz-L_w;
s_m=zeros(length(n_t),1);
t_m = t(n_t+L_w/2-1);
for k = 1:length(n_t)
    s_m(k)=mean(yTrend(n_t(k):n_t(k)+L_w-1));
end
pp=csaps(t_m,s_m);
s_pp=ppval(pp,t)';
yTrendHat = yTrend - s_pp;

figure,
subplot(211)
plot(yTrend)
hold on
plot(s_pp, 'r')
hold on
plot(n_t, s_m, 'r*')
subplot(212)
plot(yTrend-mean(yTrend)),
hold on
plot(yTrendHat, 'r')



%%
inputData = yTrend-mean(yTrend);

numberOfSamples = numel(inputData);
numSmpFactor    = factor(numberOfSamples);
[uniqueFactors, ~, ic] = unique(numSmpFactor);



% add a snippet to do it automatically.
% So far, let's pick 5^3

stepSize = 5^3;
xControlPoints    = 0:stepSize:numberOfSamples;
xControlPoints(1) = 1;

% let's smooth the signal before getting the values for the control points
filterKernel = ones(1, stepSize)./stepSize;
filteredInputData = conv(inputData, filterKernel, 'same');
yControlPoints    = filteredInputData(xControlPoints);


pp        = csaps(xControlPoints, yControlPoints);
s_pp      = ppval(pp,1:numberOfSamples);
yTrendHat = yTrend - s_pp;

figure,
plot(inputData)
hold on
plot(filteredInputData, 'r')
hold on,
plot(yTrendHat, 'g')






%% Find the trend using a spline - JL

siz = numel(yTrend);
Fe  = 1; %0.1 ? %10 ??
L_w = 200;

% convol
convWin = ones(1, L_w);
s_m     = conv(convWin, yTrend, 'same');



t = (0:siz-1)'/Fe;
n_t = 1:L_w/2:siz-L_w;
s_m=zeros(length(n_t),1);
t_m = t(n_t+L_w/2-1);
yTrend = yTrend - mean(yTrend);
for k = 1:length(n_t)
    s_m(k)=mean(yTrend(n_t(k):n_t(k)+L_w-1));
end
pp=csaps(t_m,s_m);
s_pp=ppval(pp,t)';
yTrendHat = yTrend - s_pp;

figure,
subplot(211)
plot(trend)
hold on
plot(s_pp, 'r')
hold on
plot(n_t, s_m, 'r*')
subplot(212)
plot(y-mean(y)),
hold on
plot(yTrendHat, 'r')