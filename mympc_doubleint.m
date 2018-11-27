%% Define Plant Model
A = [1 1;1 0];
B = [0 1]';
C = [1 0];
X = zeros(2 , 1);
Z = zeros(1 , 1);
%% Select Q¡¢R matrixs and predict step
Q = C' * C;
R = 1/10;
P = [1 0 ; 0 1];
mpcn = 2;
nx = length(X);
nz = length(Z);

%% Calculate matrixs H¡¢F¡¢Y
QQ = zeros(mpcn*nx , mpcn*nx);
RR = zeros(mpcn*nz , mpcn*nz);
SS = zeros(mpcn*size(B));
TT = [];
for i = 1:mpcn
    QQ((i-1)*nx+1:i*nx , (i-1)*nx+1:i*nx) = Q;
    RR((i-1)*nz+1:i*nz , (i-1)*nz+1:i*nz) = R;
    TT = [TT ; A^i];
    t = i;
    while t>0
        SS((i-1)*nx+1:i*nx , (t-1)+1:t) = A^(i-t) * B;
        t = t-1;
    end
end
QQ((mpcn-1)*nx+1:mpcn*nx , (mpcn-1)*nx+1:mpcn*nx) = P;

H = 2*(RR + SS'*QQ*SS);
FT = 2*TT'*QQ*SS;
Y = 2*(Q + TT'*QQ*TT);

x0 = [-10 0]';
y_h = C*x0;
u_h = [];
for j = 1:100
    z_c = -inv(H)*FT'*x0;
    u = z_c(1:nz , 1);
    x0 = A*x0 + B*u;
    y_h = [y_h , C*x0];
    u_h = [u_h , u];
end

plot(y_h);
hold on
stairs(u_h);



