%% Find the trend of a simple signal

% Sin-based waveform with fd = 20Hz

f  = 20;
Fs = 1e3;
t  = 1:2000;
A  = 8;
y  = (A/2)*sin(2*pi*f/Fs*t) + A;
y(t(100):t(200))  = 0;
y(t(600):t(700))  = 0;
y(t(1500):t(1600))= 0;

figure,
subplot(311)
plot(t, y, 'LineWidth', 1.5)
title('Clean signal')

% Add the 'trend' as a 1.5Hz sin plus some noise
f  = 1.5;
Fs = 1e3;
t  = 1:2000;
A  = 1;
trend = (A/2)*cos(2*pi*f/Fs*t) + 1.5*rand(1, length(t));

hold on,
subplot(312)
plot(t, trend, 'LineWidth', 1.5)
title('trend')

yTrend = y + trend;

hold on,
subplot(313)
plot(t, yTrend, 'LineWidth', 1.5)
title('''trendy''')

%% Find the trend using a spline - JL

siz = numel(yTrend);
Fe  = 1; %0.1 ? %10 ??
L_w = 200;

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