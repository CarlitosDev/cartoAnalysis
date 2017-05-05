aux = cathererTable{12,'signal'};
% aux = filterChebyShev(aux, 1e3);
n = length(aux);

figure, subplot(311),
plot(aux), axis tight;

optcor = 'coeff';
maxlag = 400;
[acor,lag]=xcorr(aux,aux,maxlag,optcor);
subplot(312), plot(lag,acor); axis tight

k=0;
lon = 500;
sal = 250;
acorlittle = zeros(size(xcorr(aux(1:lon),aux(1:lon),maxlag)));
for i=1:sal:length(aux)
    i,
    miend = i+lon-1;
    if miend<=length(aux)
        auxmini = aux(i:miend);
        [acormini,lagmini]=xcorr(auxmini,auxmini,maxlag,optcor);
        acorlittle = acorlittle+acormini;
        k = k+1;
    end
end

subplot(313)
plot(lagmini,acormini), axis tight