%THIS IS NOT FINISHED

doa1 = 20  *pi/180;
doa2 = -40 *pi/180;

doa_steer = doa2;

d = 15;

M = 8;

amp_out = 1;
mu = 0.005;

Nw = 512;
N = 6400;
%%% simulating signals
fs = 3200;
t = (1:N)/fs;
c = 343;

s1 = cos(2*pi*2.5*t);
s2 = trianglewave(10,N)*0.5;

figure(1);
plot(t,s1,t,s2)

X = zeros(M,N);

X(1,:) = s1+s2;
for m = 2:M
	X(m,:) = delay_f(s1,(m*d/c)*sin(doa1),N)+delay_f(s2,(m*d/c)*sin(doa2),N);
end

figure(2);
plot(t,X(1,:))


%%% doing Geometric Source Separation

%calculating the steering vector
w_c = zeros(M,Nw);
w_gss = zeros(M,Nw);

w = ((1:Nw)/Nw)*fs;
w_c(1,:) = ones(1,Nw);
for m = 1:M-1
	for f = 1:round(Nw/2)
	    w_c(m+1,f)=exp(-i*2*pi*w(f)*m*d/c*sin(doa_steer));    % steering vector for this frequency

	    w_c(m+1,end-f+1)=exp(i*2*pi*w(f+1)*m*d/c*sin(doa_steer));    % negative steering vector for the mirror frequency
	end
end
w_c = w_c/M;
w_gss = w_c;

%applying beamformer
o = zeros(1,N);
g = zeros(1,M-1);
this_o = zeros(1,Nw);
for k = Nw:round(Nw/2):N-round(Nw/2)
	last_o = this_o;
	this_frame = X(:,k-round(Nw/2)+1:k+round(Nw/2));
	
	this_frame_f = zeros(M,Nw);
	for m = 1:M
		this_frame_f(m,:) = fft(this_frame(m,:).*hamming(Nw)');
	end
	
	this_o_f = zeros(1,Nw);
	for f = 1:Nw
		this_o_f(f) = w_gss(:,f)'*this_frame_f(:,f);
	end
	
	%updating filters
	for f = 1:Nw
		alpha = norm(this_frame);
		Rx = this_frame(:,f)*this_frame(:,f)';
		Ry =
		E = 
		
		w_gss(:,f) = w_gss(:,f) - mu*(alpha*J1 + J2);
	end
	
	this_o = real(ifft(this_o_f));

	o(k-round(Nw/2)+1:k) = last_o(round(Nw/2)+1:end) + this_o(1:round(Nw/2));
	%o(k-Nw+1:k) = this_o;
	%o(k-round(Nw/2)+1:k+round(Nw/2)) = this_o;

	%figure(2)
	%this_f_test = fft(this_frame(2,:));
	%this_o_f = w_c(2,:).*this_f_test*M;
	%this_f_test = real(ifft(this_o_f));
	%plot(1:Nw,this_frame(2,:),1:Nw,this_f_test)
	
	%figure(3)
	%plot(1:Nw,last_o,Nw-round(Nw/2)+1:Nw+round(Nw/2),this_o,Nw-round(Nw/2)+1:Nw,o(k-round(Nw/2)+1:k))
	%pause
	
end

%o = real(ifft(o_f));
o = o*amp_out;

figure(3);
plot(t,s1,t,s2,t,o)
axis([0 1 -1.5 1.5])
