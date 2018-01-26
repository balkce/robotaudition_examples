%

doa1 = 20  *pi/180;             %direction of arrival of first signal
doa2 = -40 *pi/180;             %direction of arrival of second signal

amp_out1 = 0.8;                 %pre-amplification of first signal
amp_out2 = 1.2;                 %pre-amplification of second signal

d = 5;                          %distance between microphones in meters
M = 2;                          %number of microphones

N = 200;                        %signal size in samples

do_not_delay = 1                %flag to not delay the original signals as part of the microphone model

%%% simulating signals
t = (1:N)/N;                    %time vector (1 second)
c = 343;                        %speed of sound
fs = N;                         %sampling frequency same as signal size (1 second)

%original signals
s1 = cos(2*pi*2.5*t);
s2 = trianglewave(10,N)*0.5;

%microphones (input signals)
X = zeros(M,N);
X(1,:) = s1+s2;
for m = 2:M
	if(do_not_delay == 1)
		X(m,:) = s1*amp_out1+s2*amp_out2;
	else
		X(m,:) = delay_f(s1,(m*d/c)*sin(doa1),N)+delay_f(s2,(m*d/c)*sin(doa2),N);
	end
end

figure(1);
plot(t,s1,t,s2)


%%% doing ICA

%centering data
for m = 1:M
	X(m,:) = X(m,:) - mean(X(m,:));
end

%whitening data
R = X*X';

[U D ~] = svd(R,'econ');
DIM  = size(D,1);
T = zeros(m);
for i = 1:DIM
    T = T + 1/sqrt(D(i,i)) * U(:,i) * U(:,i)';
end
X = T * X;

rand('seed',0.0005);
w = rand(M,M);
for i = 1:M
    w(i,:) = w(i,:) / norm(w(i,:));
end

% Initialize loop variables
err = ones(M,1);
its = 0;
maxIters = 10000;
while ((max(err) > eps) && (its < maxIters))
    % Increment iteration counter
    its = its + 1;
    
    % Save last weight matrix
    w_old = w;
    
    % for each weight vector
    for i = 1:M
        % Last independent components
        si = w_old(i,:) * X;
        
        % Compute negentropy scores
        % negentropy function : f(u) = -exp(-u^2/2)
        g = si .* exp(-0.5 * (si.^2));
        gp = -1.0 * ((si.^2) .* exp(-0.5 * (si.^2)));
        
        % Update weights in the direction of maximum negentropy
        w(i,:) = mean(X .* repmat(g,m,1),2)' - mean(gp) * w_old(i,:);
        
        % Normalize weight vector
        w(i,:) = w(i,:) / norm(w(i,:));
    end
    
    % Decorrelate weight vectors
    [U,S,~] = svd(w,'econ');
    Sinv = diag(1./diag(S));
    w = U * Sinv * U' * w;
    
    % Compute innovation
    for i = 1:M
        err(i) = 1 - w(i,:) * w_old(i,:)';
    end
    
    % Display innovation
    %disp(['Iteration ' num2str(its) ': max(1 - <w' num2str(its) ',w' num2str(its-1) '>) = ' num2str(max(err))])
end

o=w'*X;

figure(2);
plot(t,o'*5)
axis([0 1 -1 1])
