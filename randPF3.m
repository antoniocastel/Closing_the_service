function park = randPF3(n)

%% code from CLOSING THE SERVICE: Contrasting Activity-Based and Time-Based Systematic Closure Policies
% % By Antonio Castellanos, Andrew Daw, Amy Ward and Galit B. Yom-Tov



%% Originally from 2024 “Conditional uniformity and Hawkes processes”. Mathematics of Operations Research 49(1):40–57 https: //doi.org/10.1287/moor.2022.1348.
    %% by Andrew Daw
    % step 1: sample uniformly at random from the factor group on 1 : n + 1
    park = randi(n+1, 1, n);
    
    parkSort = sort(park);
    spaces = zeros(1, 2*n+1);
    
    pref = parkSort(1);
    spaces(pref) = 1;
    spaces(n+1 + mod(pref, n+1)) = 1;
    
    for i = 2 : n
        
        pref = max(parkSort(i), pref + 1);
        
        while spaces(n+1 + mod(pref, n+1)) == 1
            pref = pref + 1;
        end
       

            
        spaces(pref) = 1;
        if pref < n+1
            spaces(pref + n+1) = 1;
        elseif pref > n+1
            spaces(mod(pref, n+1)) = 1;
        end

        
    end
    
    empty = sum((1:(n+1)).*(ones(1,n+1) - spaces(1:(n+1))));
    
    park = mod(park - empty, n+1); 
    
end

