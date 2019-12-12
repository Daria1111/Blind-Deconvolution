clear;clc;
N = 64; % array length

L = 3; % dimension of the PSF subspace

B = randn(N,L);
%B = randn(N,L); %random PSF subspace
alpha = randn(L,1);  % PSF coefficient              
h = B*alpha; % generate calibration vector/PSF vector

nspikes = 5; % number of spikes

% generating the spike locations satisfying a separation of 1/N
t1spikes = rand(1,nspikes);  

dmin = min(pdist(t1spikes'));
while (dmin<1/N)
fprintf('too close!\n');
t1spikes = rand(1,nspikes);        
dmin = min(pdist(t1spikes'));
end


% generate complex amplitudes of the spikes
dynamic_range=10; % in dB

x = zeros(nspikes,1);

x = exp(-1i*2*pi*rand(nspikes,1)).*(1 + 10.^(rand(nspikes,1).*(dynamic_range/20)));

% generate the frequency domain signal
y =  exp(1i*2*pi*([0:N-1]'*t1spikes))*x;



y_obs = diag(h) * y + wgn(64,1,0); % uncalirated observations


%%
addpath 'C:\Program Files\Mosek\9.1\toolbox\R2015a'
cvx_solver mosek
cvx_begin sdp quiet
    variable u(N-1,1) complex
    variable M(L,L) hermitian
    variable Z(N,L) complex
    variable t
    variable aux(N,1) complex
    dual variable dual_var;
    minimize  1/2*(t*N+trace(M) )
    subject to
          [toeplitz([t; u]), Z;
            Z' , M ]>=0
         dual_var: y_obs ==  diag(Z*B.');

cvx_end

%%


t = linspace(0,1,1e4);
% 

dualpoly = @(t)  exp(-1i*2*pi*t'*[0:N-1]) * diag(dual_var')* B/sqrt(N);
% %dualpol = zeros(length(t),1);
dualpol = sum(abs(dualpoly(t)).^2,2);
%figure, plot(dualpol)
%hold on,
%grid on;


% Isolate roots on the unit circle

for i = 1:L
    aux_u(:,i) = - conv((diag(dual_var)* conj(B(:,i))),flipud(diag(dual_var')* (B(:,i))))/N;
end

aux = sum(aux_u,2);

aux(N) = aux(N) +1 ;

roots_pol = roots((aux)); %all roots

 

% Isolate roots on the unit circle
roots_detected = roots_pol(abs(1-abs(roots_pol)) < 1e-4);
[auxsort,ind_sort]=sort(real(roots_detected));
roots_detected = roots_detected(ind_sort);

% Roots are double so take 1 and out of 2 and compute argument
t_rec = angle(roots_detected(1:2:end))/2/pi;
% Argument is between -1/2 and 1/2 so convert angles
t_rec(t_rec < 0)= t_rec(t_rec < 0) + 1;  

figure, plot(cos(2*pi*t_rec),sin(2*pi*t_rec),'*')
hold on,
plot(cos(2*pi*t1spikes), sin(2*pi*t1spikes),'ro')

legend('detected support using dualpoly','true support')
t = linspace(0,1,1e4);

plot(cos(2*pi*t), sin(2*pi*t),'b-')
grid on;