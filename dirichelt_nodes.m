function in_d = dirichelt_nodes(fix_nod, dofs)

in_d = (fix_nod(:, 1) - 1) * dofs + fix_nod(:, 2);


    