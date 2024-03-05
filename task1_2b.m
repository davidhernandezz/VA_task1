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
% displacement of the reference node [rad]
x_rot = -500*10^-6;
z_rot = 200*10^-6;

%% CALCULATIONS

resultu = transpose(function_2a(dofs));

B = [0; 0; 0; -x_rot; 0; -z_rot];

% new shims dimensions
new_dimensions = resultu\B;

fix_nod = fixnodes_2b(n_supports, dofs, new_dimensions);

% Dirichelt index vector
in_d = (fix_nod(:, 1) - 1) * dofs + fix_nod(:, 2);
% Dirichelt displacements vetor
u_d = fix_nod(:, 3);

in_n = setdiff(transpose(1:length(K)), in_d);

% gravity acceleration vector
g = [0; 0; 0; 0; 0; 0];

g_vect = repmat(g, length(K)/dofs, 1);

Fext = M * g_vect;

F_n_ext = Fext(in_n);
F_d_ext = Fext(in_d);

K_nn = K(in_n, in_n);
K_dd = K(in_d, in_d);
K_nd = K(in_n, in_d);
K_dn = K(in_d, in_n);

% displacements and forces vectors
u_n = K_nn\(F_n_ext - K_nd * u_d);

u = zeros(length(K),1);

u(in_n, 1) = u_n;
u(in_d, 1) = u_d;

u = transpose(reshape(u, [dofs, length(K)/dofs]));

F_d = K_dd*u_d + K_dn*u_n;

F = zeros(length(K),1);

F(in_n) = F_n_ext;
F(in_d) = F_d + F_d_ext;
F = transpose(reshape(F, [dofs, length(K)/dofs]));

resultu = u(ref_node, :);
resultF = F(ref_node, :);

disp(resultu);
disp(resultF);

new_dimensions = 1 + u(n_supports, 2);
