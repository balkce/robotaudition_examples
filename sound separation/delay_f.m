function y = delay_f(x,time,fs)

%disp(['Delaying ' num2str(time) ' sec.'])

%          theta = arcsin(c*t/d)
%     sin(theta) = c*t/d
% sin(theta)*d/c = t

x_f = fft(x);
N = length(x_f);

%w = ((1:N2)/N2)*(fs/2);
w = [0 1:N/2 (-N/2+1):-1]/N*fs;

y_f = zeros(1,length(x_f));
for f = 1:N
    e=exp(-i*2*pi*w(f)*time);    % steering vector for this frequency
    y_f(f) = x_f(f)*e;
end

%w = [0:N2 -N2+1:-1]*N2*2/fs;
%y_f = x_f.*exp(-1j *2*pi*w*time);

y = real(ifft(y_f));

%figure(2); plot(1:length(x),x, 1:length(x),y)
