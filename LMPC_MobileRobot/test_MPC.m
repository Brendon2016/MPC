clear
clc
%% creat a car object
dt = 0.01;
car = Bicycle( 'L' , 1.378 ,...                         % wheel distance
                'x0' , [0 0.1 0] ,...                     % init state
                'dt' , dt ,...                        % control period
                'covar', diag([0.1 0.01].^2 ) ,...      % 
                'steermax', 0.8 ,...
                'speedmax', 1);
 car_sim = Bicycle( 'L' , 1.378 ,...                    % wheel distance
                'x0' , [0 0 0] ,...                     % init state
                'dt' , dt ,...                        % control period
                'covar', diag([0 0] ) ,...              % 
                'steermax', 0.8 ,...
                'speedmax', 1);

 %% path init
 n = 5;
 steps = 1000;
 xr = zeros(3 ,steps+1);
 ur = zeros(2 ,steps);
 xr(:,1) = car_sim.x;
 steer = 0;
 speed = 0;
 for i =1:steps
     if i<=10
         speed = i*0.1;
     elseif i>(steps-10) 
         speed = (steps-i)*0.1;
     end
     if i>220 && i<240
         steer = steer + 0.05;
     elseif i>700 && i< 720
         steer = steer - 0.05;
     end
     v = speed*cos(steer);
     w = speed * sin(steer)/car_sim.L;
     ur(:,i) = [v w]';
     car_sim.step(speed , steer);
     xr(:,i+1) = car_sim.x;
 end
 
%% control parameters 
Q = [1 0 0;
     0 1 0;
     0 0 0.005];
R = 0.0000001*[0.1 0;
           0 1];
nx = 3;
nz = 2;

x_h = car.x;
u_h = [];
for i=1:steps-5
     %% calculate matrixs H F
     QQ = zeros(n*nx , n*nx);
     RR = zeros(n*nz , n*nz);
     SS = zeros(n*nx , n*nz);
     TT = zeros(n*nx , nx);
     
     Ak = zeros(nx , nx , n);                                           %k=1:n represent Ak|k to Ak+n-1|k
     Bk = zeros(nx , nz , n);
     
     for j = 1 : n
         QQ((j-1)*nx+1:j*nx , (j-1)*nx+1:j*nx) = Q;                     % calculate matrix Q¡¢R
         RR((j-1)*nz+1:j*nz , (j-1)*nz+1:j*nz) = R;
         
         Ak(:,:,j) = [1 0 -ur(1,j+i-1)*sin(xr(3,j+i-1))*dt;             % calculate Ak+j|k and Bk+j|k  
                      0 1  ur(1,j+i-1)*cos(xr(3,j+i-1))*dt;
                      0 0                                1];
         Bk(:,:,j) = [cos(xr(3,j+i-1))*dt 0;
                      sin(xr(3,j+i-1))*dt 0;
                                        0 dt];
         if j == 1                                                      % calculate matrix T
             TT(1:3,:) = Ak(:,:,j);
         else
             TT((j-1)*nx+1:j*nx,:) = TT((j-2)*nx+1:(j-1)*nx,:)*Ak(:,:,j);
         end

         t = j;                                                         % matrix SS
         while t>0
             SS((j-1)*nx+1:j*nx , (t-1)*nz+1:t*nz) = Bk(:,:,t);
             for k=(j-1):(-1):t
                 SS((j-1)*nx+1:j*nx , (t-1)*nz+1:t*nz) = Ak(:,:,k+1)*SS((j-1)*nx+1:j*nx , (t-1)*nz+1:t*nz);
             end
             t = t-1;
         end
     end
     H = 2*(RR + SS'*QQ*SS);
     F = 2*SS'*QQ*TT;
     %Y = 2*(Q + TT'*QQ*TT);
     
    %% control loop
     x = car.x - xr(:,i);
     z_c = -inv(H)*F*x;                                 % solve QP
     u = z_c(1:nz , 1)+ur(:,i);
     u_h = [u_h u];
     steer = atan(u(2)*car.L/u(1));
     speed = u(1)/cos(steer);
     odo = car.step(speed , steer);
     x_h = [x_h car.x];
end

%% plot
figure(1)
plot(x_h(1,:),x_h(2,:));
hold on
plot(xr(1,:),xr(2,:));
title("trajectory");

figure(2)
plot(x_h(1,:)-xr(1,1:end-5));
title("ex");

figure(3)
plot(x_h(2,:)-xr(2,1:end-5));
title("ey");
%hold on
% plot(xr(1,1:100),ur(1,:));
%stairs(u_h); 
