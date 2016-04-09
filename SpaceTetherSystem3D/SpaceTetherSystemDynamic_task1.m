clc; clear;

global m;
global c;
global gravity;
global tension;
global Fk;
global Fe;
global r0;

G=6.6740831e-11;
Mearth=5.9742e24;
Rearth=6378100;
omega = [0;0;pi/240/180];

m=1;
l0=500;
c=3;
N=5;

gravity=@(r) -G*Mearth*m*r/norm(r)^3;
tension=@(dr) -c*(1-l0/norm(dr))*dr;
Fk=@(v) -2*m*cross(omega,v);
Fe=@(r) -m*cross(omega, cross(omega,r));
r0 = @(t) [0;0;Rearth + 500];
%устойчивое вертикальное положение, форму тросса, точки отличаются по
%массе, последние точки отличатся по массе на 2 порядка
tspan = [1:0.2:125];
q0=cell(2*N,1);
rRand=rand(3,1)*2-1;
q0{1}=r0(0)+(l0+rand()*l0/3)*rRand/norm(rRand);
vRand=rand(3,1)*2-1;
q0{2}=0.05*l0*vRand/norm(vRand);
for i=3:2:2*N-1
    rRand=rand(3,1)*2-1;
    vRand=rand(3,1)*2-1;
    
    q0{i} = q0{i-2}+(l0+rand()*l0/10)*rRand/norm(rRand);
    q0{i+1} = 0.05*l0*vRand/norm(vRand);
end
q0=cell2mat(q0);
[t,sol]=ode45(@ODEfunc, tspan, q0);
%счетчик цикла - индекс точки, для которой строятся графики
figure(); %годографы радиусов векторов
hold on;
tmp=r0(tspan(1):tspan(2));
plot3(tmp(1),tmp(2),tmp(3), 'LineStyle', '-.', 'Marker', '.', 'Color', 'g', 'MarkerSize', 20);
for i=0:0 %N
    plot3(sol(:,i*3+1),sol(:,i*3+2),sol(:,i*3+3), 'LineStyle', '-.', 'Marker', '.', 'Color', 'k', 'MarkerSize', 10);
end

figure(); %годографы векторов скорости
hold on;
for i=0:0 %N
    plot3(sol(:,i*3+4),sol(:,i*3+5),sol(:,i*3+6),  'LineStyle', '-.', 'Marker', '.', 'Color', 'r', 'MarkerSize', 10);
end

figure('Position', [50, 50, 600, 600]);
hold on;
shft=r0(tspan);
for k=1:length(sol(:,1))
    i=0:6:6*N-1;
    plot3([shft(1),sol(k,i+1)]-shft(1),[shft(2),sol(k,i+2)]-shft(2),[shft(3),sol(k,i+3)]-shft(3), 'LineStyle', '-', 'Marker', '.', 'Color', 'b', 'MarkerSize', 10);
    axis([-N*l0 N*l0 -N*l0 N*l0 -N*l0 N*l0]/2);
    M(k) = getframe;
    clf();
end
close();

%figure('Position', [50, 50, 600, 600]);
%movie(M,5)
