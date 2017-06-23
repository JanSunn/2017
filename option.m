close all;
clear all;
N = 10000;%#samples
S=[];
r = [0.06];  %mean interest rate as initial 
C = 1350;%Strike price
r0 = 0.06; %mean interest rate
sigma = 0.1;%volatility
T = 5.0;%maturity time
tau = 1.0; %ralaxation time
gamma = 0.3; %measure of volatility
dt = 5/1000;%lentgh time step

for n=1:N
    currS = 1000;
    currR = 0.06;
    t = 0;
    while t<T
        currS = currS+currR*currS*dt+sigma*currS*sqrt(dt)*randn();
        currR = currR+(r0-currR)*dt/tau+gamma*sqrt(abs(currR))*sqrt(dt)*randn(); 
        t = t+dt;
    end
S = [S, currS];
end

 edges=[0:20:10000];

tProb = sum(S>1350)/N
cSample = exp(-r0*T)*max(S-1350, 0);
histogram(S,edges);
meanHat = sum(cSample)/N
varHat = (sum((cSample-meanHat).^2))/(N-1)