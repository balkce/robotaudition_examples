%

max_delay = 15;                 %maximum delay possible in samples

delay_x1 = 0;                   %delay of signal 1 in mic 1
delay_x2 = 10;                  %delay of signal 2 in mic 1

delay_y1 = 5;                   %delay of signal 1 in mic 2
delay_y2 = 0;                   %delay of signal 2 in mic 2

noise_w = 0;                    %noise presence (between 0 and 1)
reverb_w = 0.35;                %reverb presence (between 0 and 1)

K = 100;                        %signal size in samples
%%%%%%%%

t = 1:K;                        %time vector

%original signals
s1 = e.^(-pi*((t-6).^2)/10);
s2 = 0.75*e.^(-pi*((t-20).^2)/K);

%microphones (input signals)
x = [zeros(1,delay_x1) s1(1:end-delay_x1)] + [zeros(1,delay_x2) s2(1:end-delay_x2)];
y = [zeros(1,delay_y1) s1(1:end-delay_y1)] + [zeros(1,delay_y2) s2(1:end-delay_y2)];

%adding reverberation
x = add_reverb(x,reverb_w);
y = add_reverb(y,reverb_w);

%adding noise
x = x + randn(1,K)*noise_w/10;
y = y + randn(1,K)*noise_w/10;

figure(1); plot(t,x,t,y); axis([1 K 0 1]); title('Senales de entrada')

%%% GCC

%centering
x_c = x - mean(x);
y_c = y - mean(y);

%fft'ing input signals
x_f = fft(x_c);
y_f = fft(y_c);

%cross-correlation via Pearson
ccv_f = x_f.*conj(y_f);
ccv = real(ifft(ccv_f))/(norm(x_c)*norm(y_c));
ccv = fftshift(ccv);
mid = (size(ccv,2)/2)+1;
figure(2); plot(-max_delay:max_delay,ccv(mid-max_delay:mid+max_delay)); axis([-max_delay max_delay -0.2 1]); title('Correlacion Cruzada - Pearson') 

%cross-correlation via GCC-PHAT
ccv_fp = x_f.*conj(y_f)./abs(x_f.*conj(y_f));
ccvp = real(ifft(ccv_fp));
ccvp = fftshift(ccvp);
mid = (size(ccvp,2)/2)+1;
figure(3); plot(-max_delay:max_delay,ccvp(mid-max_delay:mid+max_delay)); axis([-max_delay max_delay -0.2 1]); title('GCC-PHAT')

