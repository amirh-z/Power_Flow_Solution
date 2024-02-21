ns = 0; 
ng = 0; 
Vm = 0; 
delta = 0; 
yload = 0; 
deltad = 0;
nbus = length(bus(:,1));
for k = 1:nbus
    n = bus(k,1);
    kb(n) = bus(k,2); 
    Vm(n) = bus(k,3);
    delta(n) = bus(k, 4);
    Pd(n) = bus(k,5); 
    Qd(n) = bus(k,6); 
    Pg(n) = bus(k,7); 
    Qg(n) = bus(k,8);
    Qmin(n) = bus(k, 9); 
    Qmax(n) = bus(k, 10);
    if Vm(n) <= 0   
        Vm(n) = 1.0;    
        V(n) = 1 + 1j*0;
    else delta(n) = delta(n)*pi/180;
        V(n) = Vm(n)*(cos(delta(n))+1j*sin(delta(n)));
        P(n) = (Pg(n)-Pd(n))/basemva;
        Q(n) = (Qg(n)-Qd(n))/basemva;
        S(n) = P(n)+1j*Q(n);
    end
end
for k = 1:nbus
    if kb(k) == 1  
        ns = ns+1;  
    else
    end
    if kb(k) == 2 
        ng = ng+1; 
    else
    end
    ngs(k) = ng;
    nss(k) = ns;
end
Ym = abs(Y);  
t = angle(Y);
m = 2*nbus-ng-2*ns;
maxerror = 1; 
converge = 1;
iter = 0;
% Start of iterations
clear A  DC   J  DX
while maxerror >= accuracy && iter <= maxiter 
    for i = 1:m
        for k = 1:m
            A(i,k) = 0;      
        end
    end
    iter = iter+1;
    for n = 1:nbus
        nn = n-nss(n);
        lm = nbus+n-ngs(n)-nss(n)-ns;
        J1 = 0;
        J2 = 0; 
        J3 = 0;
        J4 = 0;
    for i = 1:nbr
        if nl(i) == n || nr(i) == n
            if nl(i) == n 
                l = nr(i); 
            end
            if nr(i) == n 
                l = nl(i); 
            end
            J1 = J1+ Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
            J3 = J3+ Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
        if kb(n)~= 1
            J2 = J2+ Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
            J4 = J4+ Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
        else
        end
        if kb(n) ~= 1  && kb(l) ~=1
            lk = nbus+l-ngs(l)-nss(l)-ns;
            ll = l -nss(l);
      
            A(nn, ll) = -Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
            if kb(l) == 0  
                A(nn, lk) = Vm(n)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
            end
            if kb(n) == 0  
                A(lm, ll) = -Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n)+delta(l)); 
            end
            if kb(n) == 0 && kb(l) == 0  
                A(lm, lk) = -Vm(n)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
            end
        else
        end
        else 
        end
    end
    Pk = Vm(n)^2*Ym(n,n)*cos(t(n,n))+J3;
    Qk = -Vm(n)^2*Ym(n,n)*sin(t(n,n))-J1;
    if kb(n) == 1 
        P(n)= Pk;
        Q(n) = Qk; 
    end   
    if kb(n) == 2  
        Q(n)= Qk;
        if Qmax(n) ~= 0
            Qgc = Q(n)*basemva + Qd(n);
            if iter <= 7                  
                if iter > 2                
                    if Qgc  < Qmin(n)      
                        Vm(n) = Vm(n) + 0.01;   
                    elseif Qgc  > Qmax(n)   
                        Vm(n) = Vm(n) - 0.01;
                    end 
                else
                end
            else
            end
        else
        end
    end
    if kb(n) ~= 1
        A(nn,nn) = J1;  
        DC(nn) = P(n)-Pk;
    end
    if kb(n) == 0
        A(nn,lm) = 2*Vm(n)*Ym(n,n)*cos(t(n,n))+J2;  
        A(lm,nn)= J3;        
        A(lm,lm) = -2*Vm(n)*Ym(n,n)*sin(t(n,n))-J4;  
        DC(lm) = Q(n)-Qk;
    end
    end
    DX = A\DC';
    for n = 1:nbus
        nn = n-nss(n);
        lm = nbus+n-ngs(n)-nss(n)-ns;
        if kb(n) ~= 1
            delta(n) = delta(n)+DX(nn);
        end
        if kb(n) == 0
            Vm(n) = Vm(n)+DX(lm);
        end
    end
    maxerror = max(abs(DC));
    if iter == maxiter && maxerror > accuracy 
        fprintf('\nWARNING: Iterative solution did not converged after ')
        fprintf('%g', iter)
        fprintf(' iterations.\n\n')
        fprintf('Press Enter to terminate the iterations and print the results \n')
        converge = 0; 
        pause
    else
    end
   
end

if converge ~= 1
    tech = ('                      ITERATIVE SOLUTION DID NOT CONVERGE'); 
else 
    tech = ('                   Power Flow Solution by Newton-Raphson Method');
end   
V = Vm.*cos(delta)+1j*Vm.*sin(delta);
deltad = delta*180/pi;
k=0;
for n = 1:nbus
    if kb(n) == 1
        k = k+1;
        S(n)= P(n)+1j*Q(n);
        Pg(n) = P(n)*basemva + Pd(n);
        Qg(n) = Q(n)*basemva + Qd(n);
        Pgg(k) = Pg(n);
        Qgg(k) = Qg(n);     
    elseif  kb(n) ==2
        k = k+1;
        S(n) = P(n)+1j*Q(n);
        Qg(n) = Q(n)*basemva + Qd(n);
        Pgg(k) = Pg(n);
        Qgg(k) = Qg(n);  
    end
    yload(n) = (Pd(n)- 1j*Qd(n))/(basemva*Vm(n)^2);
end
busdata(:,3)=Vm';
busdata(:,4)=deltad';
Pgt = sum(Pg); 
Qgt = sum(Qg); 
Pdt = sum(Pd); 
Qdt = sum(Qd); 