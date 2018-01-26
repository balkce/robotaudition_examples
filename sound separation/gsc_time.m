%

doa1 = 20  *pi/180;
doa2 = -40 *pi/180;

doa_steer = doa1;

d = 5;

M = 5;

amp_out = 0.1;
mu0 = 0.005;

N = 1000;
Nw = 10;
%%% simulating signals
t = (1:N)/N;
KL = Nw/2;
KR = Nw/2;
c = 343;
fs = N;

s1 = cos(2*pi*2.5*t);
s2 = trianglewave(10,N)*0.5;

figure(1);
plot(t,s1,t,s2)

X = zeros(M,N);

X(1,:) = s1+s2;
for m = 2:M
	X(m,:) = delay_f(s1,(m*d/c)*sin(doa1),N)+delay_f(s2,(m*d/c)*sin(doa2),N);
end

figure(2);
plot(t,X(1,:))


%%% doing GSC
o = zeros(1,N);
o_fbf = zeros(1,N);
g = zeros(Nw+1,M-1);

delay = floor(sin(doa_steer)*d/c*fs);

for k = KL+(abs(delay)*M-1)+1+Nw:N-KR-(abs(delay)*M-1)-1
    %fixed beamforming
    o_fbf(k) = 0;
    for m=1:M
	    o_fbf(k) = o_fbf(k) + X(m,k-(delay*(m-1)));
    end
    
    %estimated noise
    u = zeros(Nw+1,M-1);
    for m=2:M
    	u(:,m-1) = X(m,k-KL-delay:k+KR-delay)' - X(m-1,k-KL:k+KR)';
    end
    
    %obtaining output sample
    o(k) = o_fbf(k) - sum(sum(g.*u));
    
    %updating filter
    for m=1:M-1
	    p_u = sum(u(:,m).^2);
	    g(:,m) = g(:,m) + (mu0/p_u)*o(k)*u(:,m);
    end
end

o = o*amp_out;

figure(3);
plot(t,s1,t,s2,t,o)
axis([0 1 -1 1])
