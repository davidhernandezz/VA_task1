function fix_nod = fixnodes_p2(n_supports, ref_node, dofs, x_rot, y_rot)

% per fer la matriu de moviments coneguts, Dirichelt, es fa igual que per
% al problema 1, pero el dof 2, y, ja que es desonegut, es la peça que s'ha
% de ficar. amés s'afagueix la rotació en x i y, dofs 4 i 5, del reference
% node, la qual es coneguda ja que es el displacement per en signe
% contrari, per a que aquest sigui 0

k = 1;
for i = 1:size(n_supports)
    for j = 1:dofs
        if j == 2

        else
            fix_nod(k, 1) = n_supports(i, 1);
            fix_nod(k, 2) = j;
            fix_nod(k, 3) = 0;
            k = k + 1;
        end
    end
end

ref_node_fix = [ref_node 4 -x_rot; ref_node 5 -y_rot];

% aixo junta les dues matrius en una sola
fix_nod = cat(1, fix_nod, ref_node_fix);
