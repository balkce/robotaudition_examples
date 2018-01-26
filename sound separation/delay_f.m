function y = delay_f(x,time,fs)

%disp(['Delaying ' num2str(time) ' sec.'])

%          theta = arcsin(c*t/d)
%     sin(theta) = c*t/d
% sin(theta)*d/c = t

x_f = fft(x);
N2 = length(x_f)/2;

w = ((0:N2)/N2)*(fs/2);
y_f = zeros(1,length(x_f));
for f = 1:N2
    e=exp(-i*2*pi*w(f)*time);    % steering vector for this frequency
    y_f(f) = x_f(f)*e;

    e=exp(i*2*pi*w(f+1)*time);    % negative steering vector for the mirror frequency
    y_f(end-f+1) = x_f(end-f+1)*e;
end

%w = [0:N2 -N2+1:-1]*N2*2/fs;
%y_f = x_f.*exp(-1j *2*pi*w*time);

y = real(ifft(y_f));

%figure(2); plot(1:length(x),x, 1:length(x),y)
