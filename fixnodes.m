function fix_nod = fixnodes(n_supports, dofs)

k = 1;
for i = 1:size(n_supports)
    for j = 1:dofs
        fix_nod(k, 1) = n_supports(i, 1);
        fix_nod(k, 2) = j;
        fix_nod(k, 3) = 0;
        k = k + 1;
    end 
end

