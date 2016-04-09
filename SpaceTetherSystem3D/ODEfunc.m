function [dq] = ODEfunc(t, q)
    global m;
    global gravity;
    global tension;
    global Fk;
    global Fe;
    global r0;
    
    l=length(q);
    q=reshape(q,[3, l/3]);
    q=num2cell(q,1);
    dq = cell(1,l/3);
    for i=1:2:l/3-1
        dq{i}=q{i+1};
    end
    %силу натяжения рассчитать заранее
    dq{2}=1/m*(gravity(q{1})  + tension(q{1}-r0(t)) + tension(q{1}-q{3})  + Fk(q{2}) + Fe(q{1}));
    for i=4:2:l/3-1
        dq{i}=1/m*(gravity(q{i-1})  + tension(q{i-1}-q{i-3}) + tension(q{i-1}-q{i+1})  + Fk(q{i}) + Fe(q{i-1}));
    end
    last = uint8(l/3);
    dq{last}=1/m*(gravity(q{last-1}) + tension(q{last-1}-q{last-3})   + Fk(q{last}) + Fe(q{last-1}));
    dq=reshape(cell2mat(dq),l,1);
end

