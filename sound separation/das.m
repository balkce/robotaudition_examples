%

doa1 = 20  *pi/180;             %direction of arrival of first signal
doa2 = -40 *pi/180;             %direction of arrival of second signal

doa_steer = doa1;               %direction to steer the beamformer (original: doa1)

d = 10;                         %distance between microphones in meters (original: 10)

M = 8;                          %number of microphones (original: 8)

N = 200;                        %signal size in samples

%%% simulating signals
t = (1:N)/N;                    %time vector (1 second)
c = 343;                        %speed of sound
fs = N;                         %sampling frequency same as signal size (1 second)

%original signals
s1 = cos(2*pi*2.5*t);
s2 = trianglewave(10,N)*0.5;

figure(1);
plot(t,s1,t,s2)

%microphones (input signals)
X = zeros(M,N);
X(1,:) = s1+s2;
for m = 2:M
	X(m,:) = delay_f(s1,(m*d/c)*sin(doa1),N)+delay_f(s2,(m*d/c)*sin(doa2),N);
end

figure(2);
plot(t,X(1,:))


%%% doing Delay-and-Sum

%calculating the steering vector
w_c = zeros(M,N);

w = ((1:N)/N)*fs;
w_c(1,:) = ones(1,N);
for m = 1:M-1
	for f = 1:round(N/2)
	    w_c(m+1,f)=exp(-i*(2*pi*w(f)*m*d/c)*sin(doa_steer));    % steering vector for this frequency

	    w_c(m+1,end-f+1)=exp(i*(2*pi*w(f+1)*m*d/c)*sin(doa_steer));    % negative steering vector for the mirror frequency
	end
end
w_c = w_c/M;

%fft
for m=1:M
	X(m,:) = fft(X(m,:));
end

%applying beamformer
o_f = zeros(1,N);
for f = 1:N
	o_f(f) = w_c(:,f)'*X(:,f);
end

o = real(ifft(o_f));

figure(3);
plot(t,s1,t,s2,t,o)
axis([0 1 -1 1])
