clear
close all
clc

%% DATA INPUT
% Load K & M matrix
load('fe_model.mat');

% Gravity acceleration vector
g = [0; 9.81*10^3; 0; 0; 0; 0];

% Define number of degrees of freedom at each node  
dofs = 6;

% define the nodes where supports are located
n_supports = [10735; 13699; 16620; 19625; 22511; 4747];


%% COMPUTATIONS
% fix_nodes matrix
for i = 1:length(n_supports)
    fix_nod((dofs*i-(dofs-1):dofs*i),1) = n_supports(i);
    fix_nod((dofs*i-(dofs-1):dofs*i),2) = 1:dofs;
    fix_nod((dofs*i-(dofs-1):dofs*i),3) = 0;
end

% Dirichelt index vector
in_d = (fix_nod(:, 1) - 1) * dofs + fix_nod(:, 2);

% Dirichelt displacements vetor
u_d = fix_nod(:, 3);

% Neumann index matrix
in_n = setdiff(transpose(1:length(K)), in_d);

% Stiffness matrix division according to index
K_nn = K(in_n, in_n);
K_dd = K(in_d, in_d);
K_nd = K(in_n, in_d);
K_dn = K(in_d, in_n);

% Eigenvectors and eigenvalues calculation
[V, D] = eigs(K_nn,M(in_n, in_n),5,'smallestabs');

% Frequency calculation
lamda = diag(D);
w = sqrt(lamda);
freq = (w/(2*pi));

