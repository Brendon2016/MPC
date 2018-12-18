%% Inverted Pendulum with Animation
% 
% This example shows how to model an inverted pendulum.  The animation is created using MATLAB(R) Handle Graphics(R).  The animation block is a masked
% S-function.  For more information, use the context menu to look under the Animation block's mask
% and open the S-function for editing.

%% model
M = .5;                 % mass of the cart
m = 0.2;                % mass of the pendulum
b = 0.1;                % coefficient of friction for cart 
I = 0.006;              % mass moment of inertia of the pendulum
g = 9.8;                % gravity acceleration
l = 0.3;                % length to pendulum center of mass

dt = 0.01;              % control cycle

p = I*(M+m)+M*m*l^2;    %denominator for the A and B matrices

A = [0      1              0           0;
     0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
     0      0              0           1;
     0 -(m*l*b)/p       m*g*l*(M+m)/p  0];
A1 = A * dt + eye(4);   % state transition matrix of discrete model
B = [     0;
     (I+m*l^2)/p;
          0;
        m*l/p];
B1 = B * dt;            % control matrix of discrete model
C = [1 0 0 0;           % measure matrix
     0 0 1 0];

D = [0;
     0];
C1 = eye(length(C));    % for simulation only
C2 = [0 0 1 0];         % for simulation only
D1 = zeros(size(B));    % for simulation only
%% control
Q = [5000 0 0 0;        % state weight matrix
     0 0.1 0 0;
     0 0 5000 0;
     0 0 0 0.1];
R = 0.1;                % control weight
[K,S,e] = lqr(A,B,Q,R); % sovle LQR


%% path
steps = 1000;
xr = zeros(4 , steps);  % reference state
ur = zeros(1,steps);    % reference control input
for i =1:steps          % generate reference state and control input
    if i <= 100
        v = i*0.001;
    elseif i> steps -100
        v = (steps - i) * 0.001;
    end
    ur(i) = v;
    xr(1,i+1) = xr(1,i) + v*dt;
end

%% simulation
steps_stabilize = 200;
x = [0 0 0 0]';              % init state
x_h = zeros(4,steps+steps_stabilize+1);         % records of state
x_h(:,1) = x;

for i=1:steps
    e = x - xr(:,i+1);          % error state of trajectory tracking control
    %e = x - [0 0 0 0]';        % state regulator
    u = -K*e;                   % calculate control input
    %e = A1*e + B1*u;
    x = A1*x + B1*(u+ur(i));    % get next state
    x_h(:,i+1) = x;
end

%% for stabilize the error of inverted pendulum caused by tracking
for i = 1:steps_stabilize
    e = x - xr(:,end);          % error state of trajectory tracking control
    u = -K*e;                   % calculate control input
    x = A1*x + B1*(u+ur(end));    % get next state
    x_h(:,steps+i+1) = x;
end
figure()
plot(x_h(3,:));
title('angle of inverted pendulum/rad');

figure()
plot(xr(1,:));
hold on
plot(x_h(1,:));
title('position of cart/m');
%open_system('penddemo');
%set_param('penddemo', 'StopTime', '10');
%sim('penddemo');