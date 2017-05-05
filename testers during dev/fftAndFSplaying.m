%
% initializations
% ===============
 T = 6 ; % period (secs/cycle)
 m = 11 ; %
 N = (2^m) ; % number of discrete data
 ffreq = 2*pi/T ; % fundamental frequency
 
 %
 fs = N/T;
 maxFreq = fs/2;
  
 t = linspace(0,T,N+1) ;
 f = 5 -2*cos( ffreq*t )...
 +3*cos( 7*ffreq*t )...
 +8*sin( 80*ffreq*t );

% think this is just for mathematical consistency...
 fhat = f ;
 fhat(1) = (f(1)+f(N+1))/2 ;
 fhat(N+1) = [] ;
 F = fft(fhat,N) ;
%
% use only the first half of the DFTs
% ===================================
%
 F=F(1:N/2) ;
 k=0:(N/2-1) ;
 omega=k*ffreq ; % in units of rads/sec
 %deltaFreq = maxFreq/(N/2);
 deltaFreq = fs/N;
 
 % extracting the coefficients
 % ---------------------------
 A = 2*real(F)/N ;
 A(1)= A(1)/2 ;
 B =-2*imag(F)/N ;
 
 frequencyAxes = omega/(2*pi*deltaFreq);
 
powerCoeffs =sqrt(A.^2 + B.^2);
%figure, plot(omega/(2*pi), powerCoeffs) 
numCoeffs=100;
figure, plot(frequencyAxes(1:numCoeffs), powerCoeffs(1:numCoeffs)) 
 
 %
 L = 20;
 fapprox = A(1)*ones(size(t)) ;
 for k=1:L
 fapprox = fapprox + A(k+1)*cos(omega(k+1)*t)...
 + B(k+1)*sin(omega(k+1)*t);
 end
 
 figure,
 plot(t, f), hold on, 
 plot(t, fapprox, 'r')