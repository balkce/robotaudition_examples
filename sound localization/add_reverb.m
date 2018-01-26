function reverb_data = add_reverb(data,reflection_coef)

if (size(data,1) > size(data,2))
	data = data';
end

fs = 1000;               %sample rate
room_dim=[4 4 2.5];     %dimensions of room (meters)
mic_pos=[1 1 1.6];      %location of mic in room
src_pos=[2 2 1];        %location of sound source in room
%reflection_coef=0.9;    %reflection coefficient
virtual_sources_num=12; %number of virtual sources: (2*virtual_sources_num+1)^3

%creating room response
room_response=rir(fs, mic_pos, virtual_sources_num, reflection_coef, room_dim, src_pos)';
plot(room_response)

%figure(4);plot(room_response)

%Adding reverb to data
reverb_data = fconv(data, room_response,1);   %adding reverb
reverb_data = reverb_data(1:size(data,2));


function [h]=rir(fs, mic, n, r, rm, src);
%RIR   Room Impulse Response.
%   [h] = RIR(FS, MIC, N, R, RM, SRC) performs a simple room impulse
%         response calculation.
%
%      FS =  sample rate.
%      MIC = row vector giving the x,y,z coordinates of
%            the microphone.  
%      N =   The program will account for (2*N+1)^3 virtual sources 
%      R =   reflection coefficient for the walls, in general -1<R<1.
%      RM =  row vector giving the dimensions of the room.  
%      SRC = row vector giving the x,y,z coordinates of 
%            the sound source.
%
%   EXAMPLE:
%
%      >>fs=44100;
%      >>rm=[20 19 21];
%      >>mic=[19 18 1.6];
%      >>src=[5 2 1];
%      >>r=0.3;
%      >>n=12;
%      >>h=rir(fs, mic, n, r, rm, src);
%
%   NOTES:
%
%   1) To implement this filter, you will need to do a fast 
%      convolution.  The program FCONV.m will do this. It is available
%      at http://www.2pi.us/code/fconv.m
%   2) All distances are in meters.
%   3) A paper has been written on this model.  It is available at:
%      http://www.2pi.us/rir.html
%      
%
%Version 3.2
%Copyright © 2003 Stephen G. McGovern

%Some of the following comments are references to equations in my paper.

nn=[-n:1:n];                          % Index for the sequence
rms=nn+0.5-0.5*(-1).^nn;              % Part of equations 2,3,& 4
srcs=(-1).^(nn);                      % part of equations 2,3,& 4
xi=[srcs*src(1)+rms*rm(1)-mic(1)];    % Equation 2 
yj=[srcs*src(2)+rms*rm(2)-mic(2)];    % Equation 3 
zk=[srcs*src(3)+rms*rm(3)-mic(3)];    % Equation 4 

[i,j,k]=meshgrid(xi,yj,zk);           % convert vectors to 3D matrices
d=sqrt(i.^2+j.^2+k.^2);               % Equation 5
time=round(fs*d/343)+1;               % Similar to Equation 6
              
[e,f,g]=meshgrid(nn, nn, nn);         % convert vectors to 3D matrices
c=r.^(abs(e)+abs(f)+abs(g));          % Equation 9
e=c./d;                               % Equivalent to Equation 10

h=full(sparse(time(:),1,e(:)));       % Equivalent to equation 11

function [y]=fconv(x, h, normalize)
%FCONV Fast Convolution
%   [y] = FCONV(x, h) convolves x and h, and normalizes the output  
%         to +-1.
%
%      x = input vector
%      h = input vector
% 
%      See also CONV
%
%   NOTES:
%
%   1) I have a short article explaining what a convolution is.  It
%      is available at http://stevem.us/fconv.html.
%
%
%Version 1.0
%Coded by: Stephen G. McGovern, 2003-2004.

if ~exist('normalize','var'); normalize = 0; end

Ly=length(x)+length(h)-1;  % 
Ly2=pow2(nextpow2(Ly));    % Find smallest power of 2 that is > Ly
X=fft(x, Ly2);		   % Fast Fourier transform
H=fft(h, Ly2);	           % Fast Fourier transform
Y=X.*H;        	           % 
y=real(ifft(Y, Ly2));      % Inverse fast Fourier transform
y=y(1:1:Ly);               % Take just the first N elements
if normalize; y=y/max(abs(y)); end;           % Normalize the output

