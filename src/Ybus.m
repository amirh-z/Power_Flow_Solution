nl = branch (:,1);
nr = branch(:,2); 
R = branch(:,3);
X = branch(:,4);

nbr = length( branch(:,1)); % number of branches
nbus =  max(max(nl), max(nr));
Z = R + 1j * X ; % branch impedance

y = ones(nbr, 1)./Z; % branch admittance
Y = zeros(nbus, nbus); % initialize Y to zero
for k = 1 : nbr % formation of the off diagonal elements
    if nl(k) > 0 && nr(k) > 0
        Y(nl(k), nr(k)) = Y(nl(k), nr(k)) - y( k );
        Y(nr(k), nl(k)) = Y(nl(k), nr(k));
    end
end
for n = 1 : nbus % formation of the diagonal elements
    for k = 1 : nbr
        if nl(k) == n || nr(k) == n
            Y(n, n) = Y(n, n) + y(k);
        else
        end
    end
end