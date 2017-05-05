function filt_x = bpfilter(x_input,fs)

% Filtra paso banda con Chebyshev orden 8 y banda de paso plana
% y rizado en la banda atenuada. 

%Frecuencias de paso 0.5 y 60 Hz
Wn = [0.2 60]*2/fs;% frecuencia de corte normalizada

rizado = 40;
orden = 8; %si es paso banda, se pone N/2

[b,a]=cheby2(0.5*orden,rizado,Wn);

[fl,cl]=size(x_input);

filt_x=zeros(fl,cl);

for k=1:fl
    filt_x(k,:)=filtfilt(b,a,x_input(k,:));    
end



% %For data sampled at 1000 Hz, design a ninth-order lowpass Chebyshev Type
% %II filter with stopband attenuation 20 dB down from the passband and a
% %stopband edge frequency of 300 Hz, which corresponds to a normalized value
% %of 0.6:
% 
% orderFilter  = 12; % f Wst is a two-element vector cheby2(n,R,Wst,'s') returns an order 2*n bandpass analog filter 
% rippleFilter = 20; % attenuation
% Wn = [0.2 60]*2/fs;
% ftype = 'bandpass';
% 
% [z,p,k] = cheby2(0.5*orderFilter,rippleFilter,Wn,ftype);
% 
% [sos,g] = zp2sos(z,p,k);   % Convert to SOS form
% Hd = dfilt.df2tsos(sos,g);  % Create a dfilt object
% h = fvtool(Hd);      % Plot magnitude response
% set(h,'Analysis','freq')   % Display frequency respon