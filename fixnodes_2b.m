function fix_nod = fixnodes_2b(n_supports, dofs, new_dimensions)
% fa la matriu de nodes, dof i desplaçament, pero imposa desplaçament en y
% per a un suport

fix_nod = zeros(dofs * length(n_supports),3);

for i = 1:length(n_supports)

    fix_nod((dofs*i-(dofs-1):dofs*i),1) = n_supports(i);
    fix_nod((dofs*i-(dofs-1):dofs*i),2) = 1:dofs;
    fix_nod((dofs*i-(dofs-1):dofs*i),3) = 0;

end

fix_nod([2;8;14;20;26;32], 3) = new_dimensions;