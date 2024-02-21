ns = 0; 
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
    else delta(n) = pi/180*delta(n);
        V(n) = Vm(n)*(cos(delta(n))+1j*sin(delta(n)));
        P(n) = (Pg(n)-Pd(n))/basemva;
        Q(n) = (Qg(n)-Qd(n))/basemva;
        S(n) = P(n) + 1j*Q(n);
    end
    if kb(n) == 1
        ns = ns+1;
    else
    end
    nss(n) = ns;
end
Ym = abs(Y); 
t = angle(Y);
ii = 0;
for ib = 1:nbus
     if kb(ib) == 0 || kb(ib) == 2
        ii = ii+1;
        jj = 0;
        for jb = 1:nbus
            if kb(jb) == 0 || kb(jb) == 2
            jj = jj+1;
            B1(ii,jj) = imag(Y(ib,jb));
            else   
            end
        end
     else
     end
end

ii = 0;
for ib = 1:nbus
     if kb(ib) == 0
        ii = ii+1;
        jj = 0;
        for jb = 1:nbus
            if kb(jb) == 0
                jj = jj+1;
                B2(ii,jj) = imag(Y(ib,jb));
            else
            end
        end
     else 
     end
end
B1inv = inv(B1); 
B2inv = inv(B2);
maxerror = 1; 
converge = 1; 
iter = 0;
% Start of iterations
while maxerror >= accuracy && iter <= maxiter 
    iter = iter+1;
    id = 0; 
    iv = 0;
    for n = 1:nbus
        nn = n-nss(n);
        J11 = 0;
        J33 = 0;
        for i = 1:nbr
            if nl(i) == n || nr(i) == n
                if nl(i) == n 
                    l = nr(i);
                end
                if nr(i) == n
                    l = nl(i);
                end
                J11 = J11+ Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
                J33 = J33+ Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
            else  
            end
        end
        Pk = Vm(n)^2*Ym(n,n)*cos(t(n,n))+J33;
        Qk = -Vm(n)^2*Ym(n,n)*sin(t(n,n))-J11;
        if kb(n) == 1 
            P(n) = Pk; 
            Q(n) = Qk; 
        end  
        if kb(n) == 2  
            Q(n) = Qk;
            Qgc = Q(n)*basemva + Qd(n);
            if Qmax(n) ~= 0
                if iter <= 20                 
                    if iter >= 10              
                        if Qgc  < Qmin(n)       
                            Vm(n) = Vm(n) + 0.005;   
                        elseif Qgc > Qmax(n)    
                            Vm(n) = Vm(n) - 0.005;
                        end 
                    else
                    end
                else
                end
            else
            end
        end
        if kb(n) ~= 1
            id = id+1;
            DP(id) = P(n)-Pk;
            DPV(id) = (P(n)-Pk)/Vm(n);
        end
        if kb(n) == 0
            iv = iv+1;
            DQ(iv) = Q(n)-Qk;
            DQV(iv) = (Q(n)-Qk)/Vm(n);
        end
    end
    Dd = -B1\DPV';
    DV = -B2\DQV';
    id = 0;
    iv = 0;
    for n = 1:nbus
        if kb(n) ~= 1
            id = id+1;
            delta(n) = delta(n)+Dd(id);
        end
    if kb(n) == 0
        iv = iv+1;
        Vm(n) = Vm(n)+DV(iv); 
    end
    end
    maxerror = max(max(abs(DP)),max(abs(DQ)));
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
    tech= ('                      ITERATIVE SOLUTION DID NOT CONVERGE'); 
else 
    tech=('                   Power Flow Solution by Fast Decoupled Method');
end   
k = 0;
V = Vm.*cos(delta)+1j*Vm.*sin(delta);
deltad = 180/pi*delta;

clear A DC DX

for n = 1:nbus
     if kb(n) == 1
        S(n)=P(n)+1j*Q(n);
        Pg(n) = P(n)*basemva + Pd(n);
        Qg(n) = Q(n)*basemva + Qd(n);
        k = k+1;
        Pgg(k) = Pg(n);
     elseif  kb(n) == 2
        S(n) = P(n)+1j*Q(n);
        Qg(n) = Q(n)*basemva + Qd(n);
        k = k+1;
        Pgg(k) = Pg(n);
     end
     yload(n) = (Pd(n)- 1j*Qd(n))/(basemva*Vm(n)^2);
end
bus(:,3) = Vm'; 
bus(:,4) = deltad';
Pgt = sum(Pg); 
Qgt = sum(Qg); 
Pdt = sum(Pd); 
Qdt = sum(Qd); 
clear Pk Qk  DP DQ J11 J33 B1 B1inv B2 B2inv DPV  DQV Dd delta ib id ii iv jb jj
