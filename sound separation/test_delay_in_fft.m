%%% simulating signals
N = 200;
d = 5;
doa1 = 20  *pi/180;
t = (1:N)/N;
c = 343;
fs = N;

s1 = cos(2*pi*2*t);
s1_d = delay_f(s1,(d/c)*sin(doa1),N);

s1_f = fft(s1);
s1_d_f = fft(s1_d);
s1_f = s1_f(1:N/2);
s1_d_f = s1_d_f(1:N/2);

figure(1);
plot([abs(s1_f); abs(s1_d_f)]')
figure(2);
plot([angle(s1_f); angle(s1_d_f)]')

