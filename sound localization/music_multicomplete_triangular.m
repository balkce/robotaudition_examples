%

doas = [30 20];

d = 20;

noise_w = 0;   %entre 0 y 1

K = 100;
%%%%%%%%

freq = 2.5;

c = 343;
fs = K;
t = (1:K)/K;

r = 2;
s1 = exp(j*(2*pi*freq*t));
s2 = exp(j*(2*pi*0.75*freq*t));

N = 3;
steer = zeros(N,r);
steer(1,:) = [1 1];    %first microphones is reference, no delay
steer(2,:) = exp(-i*2*pi*(d/c)*sin(doas*pi/180));  % second mic, delayed as usual
steer(3,:) = exp(-i*2*pi*(d/c)*sin((doas-60)*pi/180));% third mic, delays with an angular offset

%ruido
noise = sqrt(noise_w/2)*(randn(N,K)+i*randn(N,K));

%data matrix
X = steer*[s1; s2]+noise;

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
a1(2,:) = exp(-i*2*pi*(d/c)*sin(angles*pi/180));  % second mic, delayed as usual
a1(3,:) = exp(-i*2*pi*(d/c)*sin((angles-60)*pi/180));% third mic, delays with an angular offset

%compute MUSIC spectrum
for k=1:length(angles)
	music_spectrum(k)=(a1(:,k)'*a1(:,k))/(a1(:,k)'*Qn*Qn'*a1(:,k));
end

figure(2)
plot(angles,abs(music_spectrum))
