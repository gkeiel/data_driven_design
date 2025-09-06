%% initialization
clc; clear; close all; format shortg;

% data
T   = 0.01;
z   = tf('z',T);
G   = 0.5*( z-0.5 )/( z-0.9 )^2;
N   = 500;
t   = 0:T:N*T;
u   = ones(length(t),1);
y   = lsim(G,u,t)';

% reference model T_d(z)
a   = 0.3;
T_d = ( 1-a )/( z-a );


%% VRFT design
% controller C = rho'*Cb
Cb    = tf([1 0],[1 -1],T); % I structure

% calculate regressions
[b_T,a_T] = tfdata(T_d^-1,'v');
r         = zeros(1,N+1);
for k = 1:N
    r(k) = (1/(1-a))*(y(k+1) -a*y(k));
end
e           = r(1:end-1)-y(1:end-1);
L           = tf(1,1,T);
%L           = T_d*(1 -T_d);
H           = L*Cb;
[b_H,a_H] = tfdata(H,'v');
phi_L       = filter(b_H,a_H,e)';
[b_L,a_L]   = tfdata(L,'v');
u_L         = filter(b_L,a_L,u(1:end-1));

% estimate parameters
rho = (phi_L'*phi_L)\(phi_L'*u_L);
C   = rho*Cb
T_r = minreal( (C*G)/(1 +C*G),1e-5 );
zpk_Tr = zpk(T_r);