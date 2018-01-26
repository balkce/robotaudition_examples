%

doa1 = 20  *pi/180;     %direction of arrival of first signal
doa2 = -40 *pi/180;     %direction of arrival of second signal

amp_out1 = 0.35;        %pre-amplification of first signal
amp_out2 = 0.7;         %pre-amplification of second signal

d = 5;                  %distance between microphones in meters
M = 8;                  %number of microphones

N = 200;                %signal size in samples
%%% simulating signals
t = (1:N)/N;            %time vector (1 second)
c = 343;                %speed of sound
fs = N;                 %sampling frequency same as signal size (1 second)

%original signals
s1 = cos(2*pi*2.5*t);
s2 = trianglewave(10,N)*0.5;

%microphones (input signals)
X_t = zeros(M,N);
X_t(1,:) = s1+s2;
for m = 2:M
	X_t(m,:) = delay_f(s1,(m*d/c)*sin(doa1),N)+delay_f(s2,(m*d/c)*sin(doa2),N);
end

figure(1);
plot(t,s1,t,s2)


%%% doing PCA from FFT

%fft'ing data
X_f = zeros(M,N);
for n = 1:M
  X_f(n,:) = fft(X_t(n,:));
end

R = X_f*X_f';

[V,D] = eig(R);
[D,I] = sort(diag(D),1,'descend');
V = V(:,I);

o_f=V.'*X_f;

o = zeros(M,N);
for n = 1:M
  o(n,:) = real(ifft(o_f(n,:)));
end
%o = real(ifft(o_f));

figure(2);
plot(t,s1,t,s2,t,o(1,:)*amp_out1,t,o(3,:)*amp_out2)
axis([0 1 -1 1])

figure(3);
plot(t,s2,t,o(2,:)*0.5)
axis([0 1 -1 1])
