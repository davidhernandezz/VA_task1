function in_d = dirichelt_nodes(fix_nod, dofs)

k = size(fix_nod, 1);
in_d = zeros(k, 1);

in_d = (fix_nod(:, 1) - 1) * 6 + fix_nod(:, 2);


    