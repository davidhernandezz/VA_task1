function fix_nod = fixnodes_p2_1(n_supports, dofs, displ_supp, y_displ_support)
% fa la matriu de nodes, dof i desplaçament, pero imposa desplaçament en y
% per a un suport

k = 1;
for i = 1:size(n_supports)
    for j = 1:dofs
        fix_nod(k, 1) = n_supports(i, 1);
        fix_nod(k, 2) = j;
        if i == displ_supp && j == 2
            fix_nod(k, 3) = y_displ_support;
        else
            fix_nod(k, 3) = 0;
        end

        k = k + 1;
    end
end