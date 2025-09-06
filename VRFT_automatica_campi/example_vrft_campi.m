%% VRFT design (C_0 not belonging to C)
clc; clear; close all; format shortg;

% data
T   = 0.05;
num = [0 1];
den = [1 -0.6];
G   = tf(num,den,T);
N   = 50;
t   = 0:T:N*T;
u   = idinput(length(t),'prbs',[0 0.5],[-1 1]); u(1:2) = 0;
y   = lsim(G,u,t)';

% reference model M(z)
z   = tf('z',T);
M   = 0.6/( z-0.4 );

% controller Ch = rho'*C0
C0  = tf([1 0],[1 -1],T);

% ideal controller
C   = tf([0.6 -0.36],[1 -1],T);

% calculate regressions
r = zeros(1,N+1);
for k = 1:N
    r(k) = (1/0.6)*(y(k+1) -0.4*y(k));
end
e = r(1:end-1)-y(1:end-1);
L = tf(1,1,T);
%L = M*(1 -M);
[b_H,a_H] = tfdata(L*C0,'v');
phi_l     = filter(b_H,a_H,e)';
[b_L,a_L] = tfdata(L,'v');
u_l       = filter(b_L,a_L,u);

% estimate parameters
k0    = 2;
phi_L = phi_l(k0+1:end);
u_L   = u_l(k0+1:end-1); 
rho   = (phi_L'*phi_L)\(phi_L'*u_L)
Ch    = rho*C0;

% results and plots
figure(1)
T_r = Ch*G/(1+Ch*G);
step(M);
hold on;
step(T_r)
legend('M','T_r');
grid on;


%% VRFT design (C_0 belonging to C)
clc; clear; close all; format shortg;

% data
T   = 0.05;
num = [0 1];
den = [1 -0.6];
G   = tf(num,den,T);
N   = 500;
t   = 0:T:N*T;
u   = idinput(length(t),'prbs',[0 0.5],[-1 1]); u(1:2) = 0;
y   = lsim(G,u,t)';

% reference model M(z)
z   = tf('z',T);
M   = 0.6/( z-0.4 );

% controller Ch = rho'*C0
C0  = [tf([1 0],[1 -1],T); tf([0.26 -0.36],[1 -1],T)];

% ideal controller
C   = tf([0.6 -0.36],[1 -1],T);

% calculate regressions
r = zeros(1,N+1);
for k = 1:N
    r(k) = (1/0.6)*(y(k+1) -0.4*y(k));
end
e = r(1:end-1)-y(1:end-1);
L = tf(1,1,T);
H = L*C0;
[b_H1,a_H1] = tfdata(H(1),'v');
[b_H2,a_H2] = tfdata(H(2),'v');
phi_l       = [filter(b_H1,a_H1,e)', filter(b_H2,a_H2,e)'];
[b_L,a_L]   = tfdata(L,'v');
u_l         = filter(b_L,a_L,u);
phi_L       = phi_l;
u_L         = u_l(1:end-1);

% estimate parameters
rho = (phi_L'*phi_L)\(phi_L'*u_L)
Ch  = rho(1)*C0(1) +rho(2)*C0(2);

% results and plots
figure(1)
T_r = Ch*G/(1+Ch*G);
step(M);
hold on;
step(T_r)
legend('M','T_r');
grid on;