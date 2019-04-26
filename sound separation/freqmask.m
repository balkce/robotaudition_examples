%

doa1 = 20  *pi/180;             %direction of arrival of first signal
doa2 = -40 *pi/180;             %direction of arrival of second signal

doa_steer = doa1;               %direction to steer the beamformer (original: doa1)

phase_diff_threshold = 0.01     *pi/180;   %mask threshold in degrees (original:0.01)

d = 4;                         %distance between microphones in meters (original: 4)

M = 2;                          %number of microphones (only works with 2)

N = 200;                        %signal size in samples

%%% simulating signals
t = (1:N)/N;                    %time vector (1 second)
c = 343;                        %speed of sound
fs = N;                         %sampling frequency same as signal size (1 second)
wp = [0 1:N/2 (N/2-1):-1:1]/N*fs;

%original signals
s1 = cos(2*pi*2*t);
s2 = trianglewave(10,N)*0.5;

figure(1);
subplot(2,1,1)
plot(t,s1,t,s2)
subplot(2,1,2)
plot(wp,abs(fft(s1)),wp,abs(fft(s2)))

%microphones (input signals)
X = zeros(M,N);
X(1,:) = s1+s2;
X(2,:) = delay_f(s1,(d/c)*sin(doa1),N)+delay_f(s2,(d/c)*sin(doa2),N);

figure(2);
subplot(3,1,1)
plot(t,X(1,:))
subplot(3,1,2)
plot(wp,abs(fft(X(1,:))))
subplot(3,1,3)
plot(wp,arg(fft(X(1,:))),wp,arg(fft(X(2,:))))


%%% doing frequency masking

%calculating the steering vector
w_c = zeros(M,N);

w = [0 1:N/2 (-N/2+1):-1]/N*fs;
w_c(1,:) = ones(1,N);
for f = 1:N
	w_c(2,f)=exp(-i*(2*pi*w(f)*d/c)*sin(doa_steer));    % steering vector for this frequency
end

%fft
for m=1:M
	X(m,:) = fft(X(m,:));
end

o_f = zeros(1,N);
freq_mask = zeros(1,N);

o_f(1) = X(1,1);
for f = 2:N
  %aligning the other microphone given the steering vector
  align_m2 = w_c(2,f)'*X(2,f);
  
  this_m1_phase = arg(X(1,f));
  this_m2_phase = arg(align_m2);
  
  %calculating phase difference
  phase_diff = this_m1_phase - this_m2_phase;
  
  if(abs(phase_diff) < phase_diff_threshold)
    freq_mask(f) = 1;
  else
    freq_mask(f) = 0;
  end
  
  o_f(f) =freq_mask(f)*X(1,f);
end

o = real(ifft(o_f));

figure(3);
subplot(3,1,1)
plot(t,s1,t,s2,t,o)
axis([0 1 -1 1])
subplot(3,1,2)
plot(wp,abs(X(1,:)),wp,freq_mask*100,'r')
subplot(3,1,3)
plot(wp,arg(X(1,:)),wp,arg(w_c(2,:)'.*X(2,:).'))
