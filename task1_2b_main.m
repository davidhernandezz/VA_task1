clear
close all
clc

%% DATA INPUT 
% load K & M matrix
load('fe_model.mat');

% define n_dof 
dofs = 6;

% define the nodes where supports are located
n_supports = [10735; 13699; 16620; 19625; 22511; 4747];
% define the referecne node
ref_node = 1305;
% displacement of the reference node
x_rot = -50*10^-3;
y_rot = 20*10^-3;

%% CALCULATIONS

fix_nod = fixnodes_2b(n_supports, ref_node, dofs, x_rot, y_rot);

% Dirichelt index vector
in_d = (fix_nod(:, 1) - 1) * dofs + fix_nod(:, 2);
% Dirichelt displacements vetor
u_d = fix_nod(:, 3);

A = transpose(1:size(K, 1));
in_n = setdiff(A, in_d);

% gravity acceleration vector
g = [0; 9.81; 0; 0; 0; 0];

g_vect = repmat(g, 152340/dofs, 1);

F = M * g_vect;

F_n = F(in_n);



K_nn = K(in_n, in_n);
K_dd = K(in_d, in_d);
K_nd = K(in_n, in_d);
K_dn = K(in_d, in_n);

% displacements and forces vectors
u_n = K_nn\(F_n - K_nd * u_d);

u(in_n, 1) = u_n;
u(in_d, 1) = u_d;

u = transpose(reshape(u, [dofs, length(K)/dofs]));

% new shims dimensions
new_dimensions = 1000 - u(n_supports(:, 1), 2);


% shauria de revisar les unitats dels imputs, que nose si s'ha de ficar en
% mrad o rad o en que



