clear
close all
clc

%% DATA INPUT
% Load K & M matrix
load('fe_model.mat');

% Define number of degrees of freedom at each node 
dofs = 6;

% Define the nodes where supports are located
n_supports = [10735; 13699; 16620; 19625; 22511; 4747];

% Define the reference node
ref_node = 1305;

% Unit displacement in y direction for the supports
y_displ_support = 1; 

resultu = zeros(length(n_supports),dofs);
resultF = zeros(length(n_supports),dofs);

% For loop that computes the reference point displacement for a given displacement in the y axis of each of the supports
for displ_supp = 1:length(n_supports)

    %% CALCULATIONS
    
fix_nod = zeros(dofs * length(n_supports),3);

    % Fix node matrix definition
    for i = 1:length(n_supports)

        fix_nod((dofs*i-(dofs-1):dofs*i),1) = n_supports(i);
        fix_nod((dofs*i-(dofs-1):dofs*i),2) = 1:dofs;
        fix_nod((dofs*i-(dofs-1):dofs*i),3) = 0;
    end

    fix_nod(dofs * displ_supp - 4, 3) = y_displ_support;

    % Dirichelt index vector
    in_d = (fix_nod(:, 1) - 1) * dofs + fix_nod(:, 2);
    % Dirichelt displacements vetor
    u_d = fix_nod(:, 3);

    % Neumann index matrix
    in_n = setdiff(transpose(1:length(K)), in_d);

    % Gravity acceleration vector
    g = [0; 0; 0; 0; 0; 0];

    % Vector that repeats the g for all nodes
    g_vect = repmat(g, length(K)/dofs, 1);

    % Gravity force calaculation
    Fext = M * g_vect;

    % The force is divided in both index
    F_n_ext = Fext(in_n);
    F_d_ext = Fext(in_d);

    % Stiffness matrix division according to index
    K_nn = K(in_n, in_n);
    K_dd = K(in_d, in_d);
    K_nd = K(in_n, in_d);
    K_dn = K(in_d, in_n);

    % Displacements and forces vectors
    u_n = K_nn\(F_n_ext - K_nd * u_d);

    u = zeros(length(K),1);

    u(in_n, 1) = u_n;
    u(in_d, 1) = u_d;

    % Displacements reshape for better visualitation
    u = transpose(reshape(u, [dofs, length(K)/dofs]));

    % Supports' forces computation
    F_d = K_dd*u_d + K_dn*u_n;

    F = zeros(length(K),1);

    % Force vector definition
    F(in_n) = F_n_ext;
    F(in_d) = F_d + F_d_ext;
    F = transpose(reshape(F, [dofs, length(K)/dofs]));


    % Displacements of the reference node
    resultu(displ_supp,:) = u(ref_node, :);

    % Reaction forces of the reference node
    resultF(displ_supp,:) = F(ref_node, :);

end
