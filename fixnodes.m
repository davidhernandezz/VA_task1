function fix_nod = fixnodes(n_supports, dofs)

fix_nod = zeros(dofs * length(n_supports),3);

for i = 1:length(n_supports)

    fix_nod((dofs*i-(dofs-1):dofs*i),1) = n_supports(i);
    fix_nod((dofs*i-(dofs-1):dofs*i),2) = 1:dofs;
    fix_nod((dofs*i-(dofs-1):dofs*i),3) = 0;

end
