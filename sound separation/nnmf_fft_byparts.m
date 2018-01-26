%

doa1 = 20  *pi/180;             %direction of arrival of first signal
doa2 = -40 *pi/180;             %direction of arrival of second signal

amp_out1 = 0.8;                 %pre-amplification of first signal
amp_out2 = 1.2;                 %pre-amplification of second signal

d = 5;                          %distance between microphones in meters
M = 5;                          %number of microphones

N = 200;                        %signal size in samples

do_not_delay = 1                %flag to not delay the original signals as part of the microphone model

%%% simulating signals
t = (1:N)/N;                    %time vector (1 second)
c = 343;                        %speed of sound
fs = N;                         %sampling frequency same as signal size (1 second)

%original signals
s1 = cos(2*pi*2.5*t);
s2 = trianglewave(10,N)*0.5;

%microphones (input signals)
X_t = zeros(M,N);
X_t(1,:) = s1+s2;
for m = 2:M
	if(do_not_delay == 1)
		X_t(m,:) = s1*amp_out1+s2*amp_out2;
	else
		X_t(m,:) = delay_f(s1,(m*d/c)*sin(doa1),N)+delay_f(s2,(m*d/c)*sin(doa2),N);
	end
end

figure(1);
plot(t,s1,t,s2)


%%% doing NNMF

%fft'ing data
X_f = zeros(M,N);
for n = 1:M
  X_f(n,:) = fft(X_t(n,:));
end

% X should only have positive values
% no need to substract minimum values

% Initialize loop variables
o_f = (2 * rand(2,N/2))+(2j * rand(2,N/2));
o_f(abs(o_f) < 0) = 0;
W = 2 * rand(M,2);
W(W < 0) = 0;

err = Inf;

its = 0;
maxIters = 10000;

X_f2 = abs(X_f(:,1:N/2));
while ((err > eps) && (its < maxIters))
    its = its + 1;
    
    o_f_abs = abs(o_f);
    o_f_pha = angle(o_f);
    
    W = (X_f2 * o_f_abs') * pinv(o_f_abs * o_f_abs');
    W(W < 0) = 0;

    o_f_abs = pinv(W' * W) * W' * X_f2;
    o_f_abs(o_f_abs < 0) = 0;
    
    o_f = o_f_abs.*exp(j*o_f_pha);
    
    X_est = W*o_f;
    
    err = mean(sqrt(sum((X_f2-abs(X_est)).^2)));
end

o_f2 = [o_f fliplr(o_f)];

o = zeros(2,N);
for n = 1:2
  o(n,:) = real(ifft(o_f2(n,:)));
end

figure(2);
plot(t,o')
%axis([0 1 0 2])

figure(3);
plot(abs(X_f2)')
%axis([0 1 0 2])

figure(4);
plot(abs(o_f)')
%axis([0 1 0 2])
