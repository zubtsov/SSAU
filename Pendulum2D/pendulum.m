clc; clear;

global invm;
global P;
global lengths;

N=7;
g = 9.8;

masses=zeros(N,1);
lengths=zeros(N,1);

for i=1:N
    masses(i)=rand()+1;
    lengths(i)=rand()+0.5;
end

m=zeros(3*N, 3*N);
for i=1:N
    m(3*i-2,3*i-2)=masses(i);
    m(3*i-1,3*i-1)=masses(i);
    m(3*i,3*i)=masses(i)*lengths(i)^2/3;
end
invm = inv(m);

P=zeros(3*N, 1);
for i=1:N
    P(i*3-1)=-masses(i)*g;
    %вращение рассматривается относительно центра масс, поэтому момент силы
    %тяжести равен нулю
end

q0 = zeros(6*N, 1);
q0(6)=pi*(rand()-1/2);
q0(4)=sin(q0(6))*lengths(1);
q0(5)=-cos(q0(6))*lengths(1);
for i=2:N
    q0(i*6)=pi*(rand()-1/2);
    q0(i*6-2)=q0((i-1)*6-2)+sin(q0((i-1)*6))*lengths(i-1)+sin(q0(i*6))*lengths(i);
    q0(i*6-1)=q0((i-1)*6-1)-cos(q0((i-1)*6))*lengths(i-1)-cos(q0(i*6))*lengths(i);
end

tspan = [0:0.04:20];
[t,sol]=ode45(@pendulumODE, tspan, q0);

figure('Position', [50, 50, 700, 600]);
coords=zeros(N+1,2);
coords(1,:)=[0,0];
hold on;
for k=1:length(tspan)
    for j=2:N+1
        coords(j, 1)=coords(j-1, 1)+2*lengths(j-1)*sin(sol(k, (j-1)*6));
        coords(j, 2)=coords(j-1, 2)-2*lengths(j-1)*cos(sol(k, (j-1)*6));
    end
    
    plot(coords(:,1), coords(:,2), 'Marker', '.', 'Color', 'r', 'MarkerSize', 10);
    axis([-N*3.5 N*3.5 -N*5 N*2]/2);
    M(k) = getframe;
    clf();
end
%figure();
%movie(M);
movie2avi(M, 'pendulum.avi');