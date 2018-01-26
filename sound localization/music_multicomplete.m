%

doas = [40 -30];    %directional of arrival of both signals in degrees

d = 20;             %distance between microphones in meters

noise_w = 0.5;      %noise presence (between 0 and 1)

K = 100;            %signal size in samples
%%%%%%%%

freq = 2.5;         %base frequency for signal 1

c = 343;            %speed of sound
fs = K;             %sampling frequency same as signal size (1 second)
t = (1:K)/K;        %time vector (1 second)

r = 2;              %number of signals in signal sub-space

%original signals
s1 = exp(j*(2*pi*freq*t));
s2 = exp(j*(2*pi*0.75*freq*t));

N = 3;              %number of microphones
steer = zeros(N,r);
steer(1,:) = [1 1];    %first microphones is reference, no delay
steer(2,:) = exp(-i*2*pi*freq*(d/c)*sin(doas*pi/180));  % second mic, delayed one distance
steer(3,:) = exp(-i*2*pi*freq*(2*d/c)*sin(doas*pi/180));% third mic, delayed double distance

%data matrix with noise
X = steer*[s1; s2] + randn(N,K)*noise_w/10;

figure(1); plot(t,X); axis([min(t) max(t) -2 2]); title('Senales de entrada')

%covariance matrix
R = X*X'/K;

%eigendecomposicion of covariance matrix
% Q: vectors
% D: values
[Q,D] = eig(R);

%sorting eigenvalues
[D,I] = sort(diag(D),1,'descend');

%sorting eigenvectors
Q = Q(:,I);

%getting signal eigenvectors
Qs = Q(:,1:r);   %this could be done without knowing r

%getting noise eigenvectors
Qn = Q(:,r+1:N);


%%% MUSIC
%define angles to look for orthogonality
angles = -90:0.1:90;

%compute steering vectors corresponding to values in angles
a1 = zeros(N,length(angles));
a1(1,:) = ones(1,length(angles));        %first microphones is reference, no delay
a1(2,:) = exp(-i*2*pi*freq*(d/c)*sin(angles*pi/180));  % second mic, delayed one distance
a1(3,:) = exp(-i*2*pi*freq*(2*d/c)*sin(angles*pi/180));% third mic, delayed double distance

%compute MUSIC spectrum
for k=1:length(angles)
	music_spectrum(k)=(a1(:,k)'*a1(:,k))/(a1(:,k)'*Qn*Qn'*a1(:,k));
end

figure(2)
plot(angles,abs(music_spectrum))
