%%% simulating signals
amp1 = 0.1;                     %intensity of s1 in x2
amp2 = 0.5;                     %intensity of s2 in x1
N = 600;                        %signal size in samples

%time 
t = (1:N)/N;                    %time vector (1 second)

%original signals
s1 = cos(2*pi*2.5*t);
s2 = cos(2*pi*10*t);

%simulating mixed signals
shat1 = s1 + amp2*s2;
shat2 = amp1*s1 + s2;

figure(1);
plot(t,s1,t,shat1)
figure(2);
plot(t,s2,t,shat2)

%estimating SIR, SDR, SAR with bss_eval_sources
[SDR,SIR,SAR,perm]=bss_eval_sources([shat1;shat2],[s1;s2])
