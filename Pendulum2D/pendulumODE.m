function [dq] = pendulumODE(t, q)
    global lengths;
    global invm;
    global P;
    %q=[... xi' yi' phii' xi yi phii ...]'
    N=length(q)/6; %число звеньев маятника
    
    Q=zeros(2*N, 3*N);
    Q(1,1)=1; Q(1,3)=-lengths(1)*cos(q(6));
    Q(2,2)=1; Q(2,3)=-lengths(1)*sin(q(6));
    for i=2:N
        Q(2*i-1,3*(i-2)+1)=-1; 
        Q(2*i-1,3*(i-2)+3)=-lengths(i-1)*cos(q((i-1)*6));
        Q(2*i-1,3*(i-2)+4)=1; 
        Q(2*i-1,3*(i-2)+6)=-lengths(i)*cos(q(i*6));
       
        Q(2*i,3*(i-2)+2)=-1; 
        Q(2*i,3*(i-2)+3)=-lengths(i-1)*sin(q((i-1)*6));
        Q(2*i,3*(i-2)+5)=1; 
        Q(2*i,3*(i-2)+6)=-lengths(i)*sin(q(i*6)); 
    end
    
    b=zeros(2*N,1);
    b(1)=-lengths(1)*q(3)^2*sin(q(6));
    b(2)=lengths(1)*q(3)^2*cos(q(6));
    for i=2:N
        b(2*i-1)=-lengths(i)*q(i*6-3)^2*sin(q(i*6)) - lengths(i-1)*q((i-1)*6-3)^2*sin(q((i-1)*6));
        b(2*i)=lengths(i)*q(i*6-3)^2*cos(q(i*6)) + lengths(i-1)*q((i-1)*6-3)^2*cos(q((i-1)*6));
    end
    
    lambda = linsolve(Q*invm*Q', Q*invm*P-b);
    w=invm*(P-Q'*lambda);
    
    dq=zeros(6*N,1);
    for i=1:N
        dq(6*i-5:6*i-3)=w(i*3-2:i*3);
        dq(6*i-2)=q(6*i-5);
        dq(6*i-1)=q(6*i-4);
        dq(6*i)=q(6*i-3);
    end
end

