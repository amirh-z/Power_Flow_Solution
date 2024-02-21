SLT = 0;
fprintf('\n\n\n')
fprintf('                           Line Flow and Losses \n')
fprintf('     --Line--  Power at bus & line flow    --Line loss--\n')
fprintf('     from  to    MW      Mvar     MVA       MW      Mvar\n')

for n = 1:nbus
busprt = 0;
    for L = 1:nbr
        if busprt == 0
        fprintf('   \n')   
        fprintf('%6g', n) 
        fprintf('      %9.3f', P(n)*basemva)
        fprintf('%9.3f', Q(n)*basemva)
        fprintf('%9.3f\n', abs(S(n)*basemva))

        busprt = 1;
        else
        end
        if nl(L)==n      
            k = nr(L);
            In = (V(n) - V(k))*y(L);
            Ik = (V(k) - V(n))*y(L);
            Snk = V(n)*conj(In)*basemva;
            Skn = V(k)*conj(Ik)*basemva;
            SL  = Snk + Skn;
            SLT = SLT + SL;
        elseif nr(L) == n  
            k = nl(L);
            In = (V(n) - V(k))*y(L);
            Ik = (V(k) - V(n))*y(L);
            Snk = V(n)*conj(In)*basemva;
            Skn = V(k)*conj(Ik)*basemva;
            SL  = Snk + Skn;
            SLT = SLT + SL;
            else 
        end
        if nl(L)== n || nr(L)== n
            fprintf('%12g', k)
            fprintf('%9.3f', real(Snk)) 
            fprintf('%9.3f', imag(Snk))
            fprintf('%9.3f', abs(Snk))
            fprintf('%9.3f', real(SL))
        if nl(L) ==n 
            fprintf('%9.3f\n', imag(SL))
        else 
            fprintf('%9.3f\n', imag(SL))
        end
        else
        end
    end
end
SLT = SLT/2;
fprintf('   \n')
fprintf('    Total loss                         ')
fprintf('%9.3f', real(SLT))
fprintf('%9.3f\n', imag(SLT))
clear Ik In SL SLT Skn Snk