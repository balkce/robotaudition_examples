%

doa1 = 20  *pi/180;             %direction of arrival of first signal
doa2 = -40 *pi/180;             %direction of arrival of second signal

doa_steer = doa1;               %direction to steer the beamformer (original: doa1)

d = 4;                          %distance between microphones in meters (original: 4)

M = 8;                          %number of microphones (original: 8)

mu0 = 0.01;                     %base adaptation rate
mu_max = 0.0001;                %maximum adaptation rate

N = 200;                        %signal size in samples
%%% simulating signals
t = (1:N)/N;                    %time vector (1 second)
c = 343;                        %speed of sound
fs = N;                         %sampling frequency same as signal size (1 second)

%original signals
s1 = cos(2*pi*2*t);
s2 = trianglewave(10,N)*0.5;

figure(1);
plot(t,s1,t,s2)

%microphones (input signals)
X = zeros(M,N);
X(1,:) = s1+s2;
for m = 1:M-1
    X(m+1,:) = delay_f(s1,(m*d/c)*sin(doa1),N)+delay_f(s2,(m*d/c)*sin(doa2),N);
end

figure(2);
plot(t,X(1,:))


%%% doing GSC throughout time

%applying the appropriate delays to input signals
for m=1:M-1
    X(m+1,:) = delay_f(X(m+1,:),-(m*d/c)*sin(doa_steer),N);
end

%calculating the upper part of GSC
y_u = sum(X)/M;

%calculating the lower part
x_n = zeros(M-1,N);
for m=1:M-1
    x_n(m,:) = X(m+1,:)-X(m,:);
end
y_n = sum(x_n);

%applying beamformer
Nw = 16;
o = zeros(1,N);
g = zeros(M-1,Nw);
updater = zeros(M-1,Nw);
mu = 0;
for k = Nw:round(N)
    this_y_u = y_u(k);

    this_x_n = x_n(:,k-Nw+1:k);

    this_y_n = sum(sum(g.*this_x_n));

    o(k) = this_y_u - this_y_n;

    this_o = o(k-Nw+1:k);

    %calculating output power
    p_o = sum(this_o.^2);

    %updating filter (one mu per input)
    for m = 1:M-1
      %calculating noise power
      p_x_n = sum(this_x_n(m,:).^2);
      
      if(mu0*p_x_n/p_o < mu_max)
        mu = mu0/p_o;
      else
        mu = mu0/p_x_n;
      end
      
      g(m,:) = g(m,:) + mu * o(k) * this_x_n(m,:);
    end
end

figure(3);
plot(t,s1,t,s2,t,o)
axis([0 1 -1 1])
