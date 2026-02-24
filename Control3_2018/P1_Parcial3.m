clear;clc;

%% Datos
G = 6.627E-11;
R = 6.371E6;
M = 5.972E24;

h=6E5;
m=5E2;

mu = 0.15;
v0 = 9045.5;
theta = -pi/10;


dt= 1;
tf=1E5;

%% INICIALIZANDO
t= 0: dt: tf;

vx = zeros(1,length(t));
vy = zeros(1,length(t));
vz = zeros(1,length(t));

x = zeros(1,length(t));
y = zeros(1,length(t));
z = zeros(1,length(t));

ax = zeros(1,length(t));
ay = zeros(1,length(t));
az = zeros(1,length(t));

F_grav = zeros(3,length(t));
F_fric = zeros(3,length(t));


%%
x(1) = 0;
y(1) = 1.2*R;
z(1) = 0;
vx(1) = v0 * cos(theta);
vy(1) = v0 * sin(theta);
vz(1) = 0;

for i= 1:length(t)-1

    r= [x(i) y(i) z(i)];
    v= [vx(i) vy(i) vz(i)];
    Fg_i = Fg(G, M, m, r);
    Fr_i = Fr(R, h, mu, r, v);
    a = Get_a(Fg_i, Fr_i, m);
    ax(i) = a(1);
    ay(i) = a(2);
    az(i) = a(3);

    % Update position and velocity
    x(i+1) = x(i) + vx(i) * dt;
    y(i+1) = y(i) + vy(i) * dt;
    z(i+1) = z(i) + vz(i) * dt;

    vx(i+1) = vx(i) + ax(i) * dt;
    vy(i+1) = vy(i) + ay(i) * dt;
    vz(i+1) = vz(i) + az(i) * dt;

    F_grav(:,i) = transpose(Fg_i);
    F_fric(:,i) = transpose(Fr_i);
    
    if sqrt(x(i)^2 + y(i)^2 + z(i)^2) <= R
        break
    end
end

t_new= t(sqrt(x.^2 + y.^2) > R);

x_new= x(sqrt(x.^2 + y.^2) > R);
y_new= y(sqrt(x.^2 + y.^2) > R);
z_new= z(sqrt(x.^2 + y.^2) > R);


vx_new= vx(sqrt(x.^2 + y.^2) > R);
vy_new= vy(sqrt(x.^2 + y.^2) > R);
vz_new= vz(sqrt(x.^2 + y.^2) > R);

F = F_grav(:,sqrt(x.^2 + y.^2) > R);
F_grav = F;

F = F_fric(:,sqrt(x.^2 + y.^2) > R);
F_fric = F;

%%

x_t1 = -R : 1E2 : R ;
y_t1 = sqrt(R^2 - x_t1.^2);
y_t2 = -sqrt(R^2 - x_t1.^2);

x_t2 = -(R+h) : 1E2 : (R+h) ;
y_t3 = sqrt((R+h)^2 - x_t2.^2);
y_t4 = -sqrt((R+h)^2 - x_t2.^2);

figure(1)
hold on
plot(x_t1,y_t1,'b', x_t1, y_t2,'b')
plot(x_t2,y_t3,'b--', x_t2, y_t4,'b--')
plot(x_new,y_new,'g')

xlabel('X (m)') ;
ylabel('Y (m)') ;
title('Trayectoria del Satelite') ;
leg=legend('Tierra (R)','' ,'Atmosfera (h)','','Trayectoria','Location','southwest') ; 

axis equal
%%

t_impacto = t(length(x_new));
t_i_h = (t_impacto/3600) + 10;
t_i_m = ( t_i_h - fix(t_i_h) )*60;
t_i_s = ( t_i_m - fix(t_i_m) )*60;


t_i_h = fix(t_i_h);
t_i_m = fix(t_i_m);
t_i_s = fix(t_i_s);

fprintf('El tiempo de impacto es a las %.f horas %.f minutos %.f segundos +- 1 segundo\n',t_i_h,t_i_m,t_i_s) ;

x_impacto = x(t_impacto);
y_impacto = y(t_impacto);
latitud = atan(y_impacto/x_impacto)*(180/pi);

fprintf('La latitud del impacto es %.2f grados\n',latitud) ;

%%

v_new= [vx_new;vy_new;vz_new];


W_i=dot(F_fric,v_new);
[~, W ] = Int(t_new,W_i, 0);


K_inicial = 0.5*m*norm([vx_new(1); vy_new(1);vz_new(1)])^2;
K_final = 0.5*m*norm([vx_new(end); vy_new(end);vz_new(end)])^2;
U_inicial = -(G*M*m)/(norm([x_new(1); y_new(1);z_new(1)]));
U_final = -(G*M*m)/(norm([x_new(end); y_new(end);z_new(end)]));

dE =K_final - K_inicial + U_final -U_inicial;

fprintf('El trabajo que realiza el rozamiento es %.2f (J)\n', W) ;
fprintf('La Energia que se pierde en el recorrido es %.2f (J)\n',-dE) ;


%%
function [F] = Fr(R, h, mu, r,v)

    if norm(r) <= R +h
        miu= (1 -( (norm(r)-R) / h ))*mu;
    else
        miu=0;
    end

    F=-miu.*v;
end


function [F] = Fg(G, M, m, r)
    F = -G*M*m.*(r./norm(r)^3);
end

function [a] = Get_a(Fg,Fr,m)
    a= (Fg+Fr)./m;
end

%INTEGRAL
function [Int_Indef, Int_def ] = Int(x,y, y0)

Int_Indef =zeros(size(y)) ;
areas=diff(x).*(y(1:end-1)+y(2:end))/2 ;
Int_def=sum(areas) ;
Int_Indef=[0, cumsum(areas)];

Int_def=y0 + Int_def ;
Int_Indef=y0 + Int_Indef;

end