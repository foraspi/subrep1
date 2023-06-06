% System matrices to study linear-qadratic dynamics in absolute values
% clear all

%--------------------------- System Parameters ----------------------------
% Run SMIB_2states_gen ???? or set manually
%------------------------ End of System Parameters ------------------------


%---------------------------- System Equations ----------------------------
syms x [2 1]
syms x_op [2 1]
syms u
syms u_op
f1(x) = -K_D / (2 * H) * x(1) -  abs(E_htc) * E_Bmod * sin(x(2))/(2*X_T*H);
f2(x) = 2 * pi * sys_freq * x(1);
g1(u) = u / (2 * H);
g2(u) = u * 0;
f = {f1;f2};
g = {g1;g2};
%------------------------ End of System Equations -------------------------


%----------------------- Taylor Polynom Method ----------------------------
tay_f = taylor(f,x,x_op,'Order',3);
x_kr = kron(x,x);
sx = string(x);
sx_kr = string(x_kr);
A_lin = zeros(size(tay_f,1),size(x,1));
Asqr = zeros(size(tay_f,1),size(x_kr,1));
for i = 1:size(tay_f,1)
    [c, t] = coeffs(tay_f(i), x);
    c = double(subs(c,x_op,x0));
    st = string(t);
    if ~isempty(c)
        A_lin (i,:) = get_lin(c,st,sx);
        Asqr (i,:) = get_sqr(c,st,sx_kr);
    end
end
tay_g = taylor(g,u,u_op,'Order',3);
su = string(u);
B_lin = zeros(size(tay_g,1),size(u,1));
for i = 1:size(tay_g,1)
    [c, t] = coeffs(tay_g(i), u);
    c = double(subs(c,u_op,u0));
    st = string(t);
    if ~isempty(c)
        B_lin (i,:) = get_lin(c,st,su);
    end 
end
%-------------------- End of Taylor Polynom Method ------------------------


%----------------------- Construction of Matrices -------------------------
n = size(A_lin,1);
n2 = n + n * n;
I = eye(n, n);
A21 = kron(A_lin, I) + kron(I, A_lin);
A_lq = zeros(n + n^2, n + n^2);
A_lq(1:n, 1:n) = A_lin;
A_lq(1:n, n+1:n2) = Asqr;
A_lq(n+1:n2, n+1:n2) = A21;
B_lq = zeros(n + n^2, size(u0,1));
B_lq(1:n,:) = B_lin;
N1 = zeros(size(A_lq));
N1(n+1:n2, 1:n) = kron(B_lin, I) + kron(I, B_lin);
Cx_lin = double(subs(subs(tay_f,x_op,x0),x,zeros(size(x))));
Cu_lin = double(subs(subs(tay_g,u_op,u0),u,zeros(size(u))));
Cx_lq = zeros(n2,1);
Cu_lq = zeros(n2,1);
Cx_lq(1:n) = Cx_lin;
Cu_lq(1:n) = Cu_lin;
%-------------------- End of Construction of Matrices ---------------------

%%
%----------------------- Time Domain Simulation ---------------------------

dt = 1e-4; % time step
t = 0:dt:10.0; % simulation times
K_u = 0.5; % magnitude of the input perturbation
freq_u = 0.5; % frequency of the input perturbation

% начальные услови€
x_init = [0; delta_0r];
u_init = Tm;
x0 = zeros(size(A_lin,1),1);
x0(2) = 0.0; % задаЄм возмущение второй переменной состо€ни€
u0 = 0;


t_size = size(t, 2);

%%                       ---- Linear System ----
x = x_init + x0;
u = u_init + u0;

y_lin = zeros(t_size, 2);
for i = 1:t_size
    y_lin(i, :) = x;
    dx = A_lin * x + B_lin * u + Cx_lin + Cu_lin;
    x = x + dx * dt;
end
% figure()
% plot(t, y_lin)
% title ('ѕереходной процесс линейной системы(u = 0.5 * sin(t))');
% legend('omega_r', 'delta', 'psi_{fd}', 'psi_{1d}', 'psi_{1q}', 'psi_{2q}');

%%                       ---- Qadratic System ----
a_size2 = size(A_lq, 1);
y_bilin2 = zeros(t_size, 2);
x = [x_init + x0; kron(x_init + x0,x_init + x0)];
u = u_init + u0;

% возмущение через управление
%u_1 = 0;
%u_2 = K_u * sin(freq_u*t);


for i = 1:t_size
    y_bilin2(i, :) = x(1:2);
    dx = A_lq * x + N1 * x * u+ B_lq * u + Cx_lq + Cu_lq;
    x = x + dx * dt;
end

% figure()
% plot(t, y_bilin2)
% title ({['ѕереходной процесс билинейной системы: ']; ['u = ',num2str(K_u) ' * sin(t)']});
% legend('omega_r', 'delta', 'psi_{fd}', 'psi_{1d}', 'psi_{1q}', 'psi_{2q}');

%                       ---- Nonlinear System ----
y_nolin = zeros(t_size, 2);
x_init = [0; delta_0r];
x = x_init;
u = 0; %вектор возмущени€

% возмущение через управление
%u_3 = E_fd0 + K_u * sin(freq_u*t);
%u_1 = P;

% возмущение ненулевыми начальными услови€ми
x = x + x0;
u_1 = P;

for i = 1:t_size
    y_nolin(i, :) = x;
    u(1) = u_1;
    dx = fsolvseq(sys_eq,x,u);
    x = x + dx * dt;
end
% figure()
% plot(t, y_nolin)
% title ('ѕереходной процесс нелинейной системы(u = E_{fd0} + 0.1 * sin(t))');
% legend('omega_r', 'delta', 'psi_{fd}', 'psi_{1d}', 'psi_{1q}', 'psi_{2q}');

figure()
%plot(t,y_lin, t,y_bilin2, t, y_nolin - [domega_0 delta_0r])
plot(t,y_lin, t,y_bilin2, t, y_nolin)
title (num2str(x0(2)))
% figure()
% plot(t,y_lin, t,y_bilin2, t, y_nolin)


%============================== FUNCTIONS =================================
function c_lin = get_lin(c,t,x)
    % get linear coefficients
    c_lin = zeros (1,size(x,1));
    for i_x = 1 : size(x,1)
        id = find(strcmp(x(i_x),t));
        if ~isempty(id)
            c_lin(i_x) = c(id);
        end    
    end
end

function c_sqr = get_sqr (c,t,x_kr)
    % get quadratic coefficients
    exception_x_kr = string;
    c_sqr = zeros (1,size(x_kr,1));
    for i_x = 1 : size(x_kr,1)
        if sum(x_kr(i_x) == exception_x_kr) == 1 % if the element from x_kr has already been
            continue %  then we move on to the next step
        end
        tmp_x = x_kr(i_x);
        tf = strcmp(tmp_x,t);
        id = find(tf);
        if ~isempty(id)
            c_sqr(i_x) = c(id);
        end
        exception_x_kr(i_x) = tmp_x;
    end
end

%-----------------------SOLVE the DYNAMIC EQUATIONS------------------------
function [fx] = fsolvseq(eq,x,u)
    
    % Solving EPS dynamic equations
    domega = eq{1,1}(x,u);
    ddelta = eq{2,1}(x,u);
    % Constructing of the state variables increments vector
    fx = [domega; ddelta];
    
end
%--------------------END of SOLVE the DINAMIC EQUATIONS--------------------