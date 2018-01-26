%
angle = 30;     %directional of arrival of the signal in degrees
d = 20;         %distance between microphones in meters

noise_w = 0;    %noise presence (between 0 and 1)

K = 200;        %signal size in samples
%%%%%%%%

freq = 2.5;     %base frequency of signal

c = 343;        %speed of sound
t = (1:K)/K;    %time vector (1 second)

r = 1;          %number of signals in signal sub-space

s1 = exp(j*(2*pi*freq*t));  %defining the original signal

N = 2;          %number of microphones

x = s1; %first mic, steering vector equal to 1, no delay
y = s1*exp(-i*2*pi*freq*(d/c)*sin(angle*pi/180));   % second mic, delayed one distance

%adding noise
x = x + randn(1,K)*noise_w/10;
y = y + randn(1,K)*noise_w/10;


figure(1); plot(t,x,t,y); axis([min(t) max(t) -1 1]); title('Senales de entrada')

%data matrix
X = [x; y];

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
a1(1,:) = ones(1,length(angles)); %first microphones is reference, no delay
a1(2,:) = exp(-i*2*pi*freq*(d/c)*sin(angles*pi/180));   % second mic, delayed one distance

%compute MUSIC spectrum
for k=1:length(angles)
	music_spectrum(k)=(a1(:,k)'*a1(:,k))/(a1(:,k)'*Qn*Qn'*a1(:,k));
end

figure(2)
plot(angles,abs(music_spectrum));  title('MUSIC')
