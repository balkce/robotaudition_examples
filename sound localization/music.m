%

doas = [50]*pi/180;     %directional of arrival of both signals in degrees

N = 3;                  %number of microphones
K = 20;                 %signal size in samples
d = 0.5;                %distance between microphones in meters
noise_var = 0.1;        %noise presence (between 0 and 1)
r = length(doas);       %number of signals in signal sub-space

%steering vectors
A = exp(-i*2*pi*d*(0:N-1)'*sin(doas));

%signals
sig = round(rand(r,K))*2-1;

%noise (uncorrelated)
noise = sqrt(noise_var/2)*(randn(N,K)+i*randn(N,K));

%data matrix
X = A*sig+noise;

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
a1 = exp(-i*2*pi*d*(0:N-1)'*sin(angles*pi/180));

%compute MUSIC spectrum
for k=1:length(angles)
	music_spectrum(k)=(a1(:,k)'*a1(:,k))/(a1(:,k)'*Qn*Qn'*a1(:,k));
end

figure(1)
plot(1:K,X(1,:),1:K,X(2,:),1:K,X(3,:))

figure(2)
plot(angles,abs(music_spectrum))
