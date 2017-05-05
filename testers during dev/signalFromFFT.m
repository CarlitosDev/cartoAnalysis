idx = 2

idx = 7
data = dataIn.signalData(idx, :);


nmbOfCoeffs  = 80;

N = length ( data ); % data is assumed to be 1D

if ( nmbOfCoeffs >= N )
    nmbOfCoeffs = N-1;
end
mInd = 1:nmbOfCoeffs;


%   FFT
dataFFT = fft ( data ) / N;


%   FS Coefficients
C0 = real (  dataFFT(1) ) ;
Ck = 2*dataFFT( mInd + 1 );


figure,
stem(abs(Ck))



%   II - Reconstruct from Ck

tSteps = N;
tR = 0:N/(N-1):N; % tReconstruction 
FS_data = C0 * ones ( tSteps,1 );
for k = 1:nmbOfCoeffs
    for n = 1:tSteps
        FS_data(n) = FS_data(n) + real ( Ck(k) * exp(j*2*pi*k*tR(n)/N) );
    end
end


figure,
xaxis = 1:N;
plot(xaxis,data,tR,FS_data,'r')

% % correct drift
% alignedData = data - FS_data;

residuals = data - FS_data';
figure,
plot(xaxis, residuals);