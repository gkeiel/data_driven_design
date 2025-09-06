%% initialization
clc; clear; close all; format shortg;

% data
T   = 1;
z   = tf('z',T);
G   = 0.5/( z-0.9 );
u   = [1; 1];
N   = 2;
y   = zeros(1,N+1);
r   = zeros(1,N+1);
for k = 1:N
    y(k+1) = 0.9*y(k) + 0.5*u(k);
end

% reference model T_d(z)
a   = 0.6;
T_d = ( 1-a )/( z-a );


%% VRFT design
% controller C = rho'*Cb
Cb    = [tf([1 0],[1 -1],T); tf([0 1],[1 -1],T)]; % PI structure

% calculate regressions
for k = 1:N
    r(k) = (1/(1-a))*(y(k+1) -a*y(k));
end
e           = r(1:end-1)-y(1:end-1);
L           = tf(1,1,T);
H           = L*Cb;
[b_H1,a_H1] = tfdata(H(1),'v');
[b_H2,a_H2] = tfdata(H(2),'v');
phi_L       = [filter(b_H1,a_H1,e)', filter(b_H2,a_H2,e)']
[b_L,a_L]   = tfdata(L,'v');
u_L         = filter(b_L,a_L,u);

% estimate parameters
rho = (phi_L'*phi_L)\(phi_L'*u_L)
C   = minreal( rho(1)*Cb(1) +rho(2)*Cb(2) )