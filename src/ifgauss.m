Vm = 0;
delta = 0;
yload = 0;
deltad = 0;
nbus = length(bus(:,1));
for k = 1:nbus
    n = bus(k,1);
    kb(n) = bus(k,2);
    Vm(n) = bus(k,3);
    delta(n) = bus(k,4);
    Pd(n) = bus(k,5); 
    Qd(n) = bus(k,6); 
    Pg(n) = bus(k,7); 
    Qg(n) = bus(k,8);
    Qmin(n) = bus(k,9); 
    Qmax(n) = bus(k,10);
    if Vm(n) <= 0  
        Vm(n) = 1.0; 
        V(n) = 1 + 1j*0;
    else
        delta(n) = delta(n)*pi/180;
        V(n) = Vm(n)*(cos(delta(n)) + 1j*sin(delta(n)));
        P(n) = (Pg(n)-Pd(n))/basemva;
        Q(n) = (Qg(n)-Qd(n))/basemva;
        S(n) = P(n) + 1j*Q(n);
    end
    DV(n) = 0;
end
num = 0;
AcurBus = 0;
converge = 1;
Vc = zeros(nbus,1)+1j*zeros(nbus,1);
Sc = zeros(nbus,1)+1j*zeros(nbus,1);

iter = 0;
maxerror = 10;
while maxerror >= accuracy && iter <= maxiter
    iter = iter+1;
    for n = 1:nbus
        YV = 0+1j*0;
    for L = 1:nbr
        if nl(L) == n
            k = nr(L);
            YV = YV + Y(n,k)*V(k);
        elseif nr(L) == n
            k = nl(L);
            YV = YV + Y(n,k)*V(k);
        end
    end
    Sc = conj(V(n))*(Y(n,n)*V(n) + YV) ;
    Sc = conj(Sc);
    DP(n) = P(n) - real(Sc);
    DQ(n) = Q(n) - imag(Sc);
    if kb(n) == 1
        S(n) = Sc; 
        P(n) = real(Sc); 
        Q(n) = imag(Sc); 
        DP(n) = 0; 
        DQ(n) = 0;
        Vc(n) = V(n);
    elseif kb(n) == 2
        Q(n) = imag(Sc); 
        S(n) = P(n) + 1j*Q(n);
        if Qmax(n) ~= 0
            Qgc = Q(n)*basemva + Qd(n);
                if abs(DQ(n)) <= 0.005 && iter >= 10 
                    if DV(n) <= 0.045                
                        if Qgc < Qmin(n)           
                            Vm(n) = Vm(n) + 0.005;       
                            DV(n) = DV(n)+0.005;          
                        elseif Qgc > Qmax(n)         
                            Vm(n) = Vm(n) - 0.005;        
                            DV(n) = DV(n)+0.005; 
                        end
                    else
                    end
                else
                end
        else
        end
    end
    if kb(n) ~= 1
        Vc(n) = (conj(S(n))/conj(V(n)) - YV )/ Y(n,n);
    else   
    end
    if kb(n) == 0
        V(n) = V(n) + accel*(Vc(n)-V(n));
    elseif kb(n) == 2
        VcI = imag(Vc(n));
        VcR = sqrt(Vm(n)^2 - VcI^2);
        Vc(n) = VcR + 1j*VcI;
        V(n) = V(n) + accel*(Vc(n) -V(n));
    end
    end
    maxerror = max( max(abs(real(DP))), max(abs(imag(DQ))) );
    if iter == maxiter && maxerror > accuracy
        fprintf('\nWARNING: Iterative solution did not converged after ')
        fprintf('%g', iter), fprintf(' iterations.\n\n')
        fprintf('Press Enter to terminate the iterations and print the results \n')
        converge = 0;   
        pause  
    else
    end
end
if converge ~= 1
   tech= ('                      ITERATIVE SOLUTION DID NOT CONVERGE'); 
else 
   tech=('                   Power Flow Solution by Gauss-Seidel Method');
end   
k=0;
for n = 1:nbus
     Vm(n) = abs(V(n));
     deltad(n) = angle(V(n))*180/pi;
     if kb(n) == 1
        S(n) = P(n)+1j*Q(n);
        Pg(n) = P(n)*basemva + Pd(n);
        Qg(n) = Q(n)*basemva + Qd(n);
        k = k+1;
        Pgg(k) = Pg(n);
     elseif  kb(n) == 2
        k = k+1;
        Pgg(k) = Pg(n);
        S(n) = P(n)+1j*Q(n);
        Qg(n) = Q(n)*basemva + Qd(n);
     end
     yload(n) = (Pd(n)- 1j*Qd(n))/(basemva*Vm(n)^2);
end 
Pgt = sum(Pg);  
Qgt = sum(Qg);
Pdt = sum(Pd); 
Qdt = sum(Qd); 
bus(:,3) = Vm';
bus(:,4) = deltad';
clear  AcurBus  DP  DQ  DV  L Sc Vc VcI VcR YV converge delta