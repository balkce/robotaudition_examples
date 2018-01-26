%
doas = [50 20]; %directional of arrival of both signals in degrees

c = 343;        %speed of sound
N = 3;          %number of microphones
K = 100;        %signal size in samples
d = 20;         %distance between microphones in meters
noise_w = 0;    %noise presence (between 0 and 1)
freq = 2.5;     %base frequency for signal 1

t = (1:K)/K;    %time vector (1 second)
r = length(doas); %number of signals in signal sub-space

%steering vectors
steer = exp(-i*2*pi*(d/c)*(0:N-1)'*sin(doas*pi/180));

%original signals
s1 = exp(j*(2*pi*freq*t));
s2 = exp(j*(2*pi*0.75*freq*t));

sig = [s1; s2];

%noise (uncorrelated)
noise = sqrt(noise_w/2)*(randn(N,K)+i*randn(N,K));

%data matrix
X = steer*sig+noise;

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
a1 = exp(-i*2*pi*(d/c)*(0:N-1)'*sin(angles*pi/180));

%compute MUSIC spectrum
for k=1:length(angles)
	music_spectrum(k)=(a1(:,k)'*a1(:,k))/(a1(:,k)'*Qn*Qn'*a1(:,k));
end

figure(1)
plot(t,X)

figure(2)
plot(angles,abs(music_spectrum))
