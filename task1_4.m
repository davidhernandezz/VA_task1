clear
close all
clc


%% input data definition
% load K & M matrix
load('fe_model.mat');

% gravity acceleration vector
g = [9.81*10^3; 0; 0; 0; 0; 0];

% define n_dof 
dofs = 6;

% define the nodes where supports are located
n_supports = [10735; 13699; 16620; 19625; 22511; 4747];


%% Preliminary computations
% fix_nodes matrix
fix_nod = zeros(dofs * length(n_supports),3);

for i = 1:length(n_supports)
    fix_nod((dofs*i-(dofs-1):dofs*i),1) = n_supports(i);
    fix_nod((dofs*i-(dofs-1):dofs*i),2) = 1:dofs;
    fix_nod((dofs*i-(dofs-1):dofs*i),3) = 0;
end

% Dirichelt index vector
in_d = (fix_nod(:, 1) - 1) * dofs + fix_nod(:, 2);
% Dirichelt displacements vetor
u_d = fix_nod(:, 3);

in_n = setdiff(transpose(1:length(K)), in_d);

g_vect = repmat(g, length(M)/dofs, 1);

Fext = M * g_vect;

F_n_ext = Fext(in_n);
F_d_ext = Fext(in_d);

K_nn = K(in_n, in_n);
K_dd = K(in_d, in_d);
K_nd = K(in_n, in_d);
K_dn = K(in_d, in_n);

% displacements and forces vectors
u_n = K_nn\(F_n_ext - K_nd * u_d);

u(in_n, 1) = u_n;
u(in_d, 1) = u_d;

u = transpose(reshape(u, [dofs, length(K)/dofs]));

F_d = K_dd*u_d + K_dn*u_n;

F(in_n) = F_n_ext;
F(in_d) = F_d + F_d_ext;

[V, D] = eigs(K_nn,M(in_n, in_n),11,'smallestabs');
lamda = diag(D);
w = sqrt(lamda);
freq = (w/(2*pi));
%%%%%%%%
F_n = F_n_ext;
phi = V;
freq_mfr = (0:10:2000);
b = 0.02;
for i = 1:size(freq_mfr, 2)
    w_mfr = 2 * pi * freq_mfr(1, i);
    for j = 1:size(V, 2)
        F_mfr = transpose(phi(:, j)) * F_n;
        mass_term = transpose(phi(:, j)) * M(in_n, in_n) * phi(:, j);
        stiffness_term = transpose(phi(:, j)) * K_nn * phi(:, j);
        dump_term = 2 * mass_term * b * 2 * pi * freq(j, 1);

        xi(j, 1) = (F_mfr) / (-w_mfr^2*mass_term + stiffness_term + (1i*dump_term*w_mfr));
    end
    x(:, i) = phi * xi;
end

ref = (1305 - 1) * 6 + 1;
figure(1)
plot(freq_mfr, abs(x(ref, :)))
xlabel('Frequency [Hz]')
ylabel('x displacement')
grid on 
grid minor


%%%%%%%%
% phi = V;
% F_mfr = transpose(phi) * F_n_ext;
% 
% mass_term = transpose(phi) * M(in_n, in_n) * phi;
% stiff_term = transpose(phi) * K_nn * phi;
% b = 0.02;
% 
% freq_mfr = (0:2:2000);
% xi = zeros(1, size(V, 2));
% for j = 1:size(freq_mfr, 2)
%     w_mfr(1, j) = freq_mfr(1, j) * 2*pi;
% 
%     mass_term_final = -mass_term * w_mfr(1, j)^2;
%     dumping = 1i* b * w_mfr(1, j) * 2 *mass_term_final * ;
% 
%     for k = 1:size(V, 2)
%         xi(k, j) = (F_mfr(k, 1) / mass_term_final(k, k) + stiff_term(k, k) + dumping);
%     end
% end
% 
% x = phi * xi;
% 
% ref = (1305 - 1) * 6 + 1;
% plot(freq_mfr, x(ref, :))
