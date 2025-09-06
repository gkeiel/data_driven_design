%% initialization
clc; clear; close all; format shortg;

% data
T   = 1;
z   = tf('z',T);
G   = ( -z+1.1 )*( z-0.8 )/( z-0.9 )^3;
N   = 600;
t   = 0:T:N*T;
u   = 0.5 +0.5*square(2*pi*t/200);
y   = lsim(G,u,t)';

% reference model T_d(z)
T_d = ( 0.125*z^2 )/( z-0.5 )^3;


%% VRFT design
% controller C = rho'*Cb
Cb = [1;
      tf([1 0],[1 -1],T);
      tf([1 -1],[1 0],T)];
  
% calculate e_v(k)
r_v = zeros(1,N+1);
for k = 3:N
    r_v(k) = 8*( y(k+1) -1.5*y(k) +0.75*y(k-1) -0.125*y(k-2) );
end
e_v = r_v-y;

% calculate regressions
L         = tf(1,1,T);
L         = minreal( T_d*(1 -T_d) );
[b_L,a_L] = tfdata(L,'v');
[b_1,a_1] = tfdata(L*Cb(1),'v');
[b_2,a_2] = tfdata(L*Cb(2),'v');
[b_3,a_3] = tfdata(L*Cb(3),'v');
phi_l     = [filter(b_1,a_1,e_v)', filter(b_2,a_2,e_v)', filter(b_3,a_3,e_v)'];
u_l       = filter(b_L,a_L,u');

% estimate parameters
k_0   = 1;
phi_L = phi_l(k_0:end,:);
u_L   = u_l(k_0:end);
rho   = (phi_L'*phi_L)\(phi_L'*u_L);

% results
C      = rho(1)*Cb(1)+rho(2)*Cb(2)+rho(3)*Cb(3);
T_r    = feedback(C*G,1);
zpk_C  = zpk(C)
zpk_Tr = zpk(T_r);