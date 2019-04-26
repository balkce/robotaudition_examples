%

doa1 = 20  *pi/180;             %direction of arrival of first signal
doa2 = -40 *pi/180;             %direction of arrival of second signal

amp_out1 = 0.8;                 %pre-amplification of first signal
amp_out2 = 1.2;                 %pre-amplification of second signal

d = 4;                          %distance between microphones in meters
M = 8;                          %number of microphones

N = 200;                        %signal size in samples

do_not_delay = 1                %flag to not delay the original signals as part of the microphone model

%%% simulating signals
t = (1:N)/N;                    %time vector (1 second)
c = 343;                        %speed of sound
fs = N;                         %sampling frequency same as signal size (1 second)

%original signals
s1 = cos(2*pi*2*t);
s2 = trianglewave(10,N)*0.5;

%microphones (input signals)
X = zeros(M,N);
X(1,:) = s1+s2;
for m = 1:M-1
	if(do_not_delay == 1)
		X(m+1,:) = s1*amp_out1+s2*amp_out2;
	else
		X(m+1,:) = delay_f(s1,(m*d/c)*sin(doa1),N)+delay_f(s2,(m*d/c)*sin(doa2),N);
	end
end

figure(1);
plot(t,s1,t,s2)


%%% doing NNMF

% X should only have positive values

X = X - min(min(X));

% Initialize loop variables
rand('seed',0.0001)
o = 2 * rand(2,N);
W = 2 * rand(2,M);

err = Inf;

its = 0;
maxIters = 5000;

X_inv = pinv(X);
while ((err > eps) && (its < maxIters))
    its = its + 1;
    
    W = (X * o') * pinv(o * o');
    W(W < 0) = 0;

    o = pinv(W' * W) * W' * X;
    o(o < 0) = 0;
    
    X_est = W*o;
    
    err = mean(sqrt(sum((X-X_est).^2)));
end

figure(2);
plot(t,o')
%axis([0 1 0 2])
