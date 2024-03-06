clear
close all
clc
% This code takes some time to fully compile

%% DATA INPUT
% load K & M matrix
load('fe_model.mat');

% gravity acceleration vector
g = [9.81*10^3; 0; 0; 0; 0; 0];

% define n_dof 
dofs = 6;

% define the nodes where supports are located
n_supports = [10735; 13699; 16620; 19625; 22511; 4747];


%% CALCULATIONS
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

% Neumann index matrix
in_n = setdiff(transpose(1:length(K)), in_d);

% Vector that repeats the g for all nodes
g_vect = repmat(g, length(M)/dofs, 1);

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

% displacements and forces vectors
u_n = K_nn\(F_n_ext - K_nd * u_d);

u(in_n, 1) = u_n;
u(in_d, 1) = u_d;

u = transpose(reshape(u, [dofs, length(K)/dofs]));

% Supports' forces computation
F_d = K_dd*u_d + K_dn*u_n;

% Force vector definition
F(in_n) = F_n_ext;
F(in_d) = F_d + F_d_ext;


% Eigenmodes & eigenvalues calc, for modal shapes and frequencies
[V, D] = eigs(K_nn,M(in_n, in_n),11,'smallestabs');

% Frequencies calculation from eigenvalues
lamda = diag(D);
w = sqrt(lamda);
freq = (w/(2*pi));

F_n = F_n_ext;

% Transform matrix is modal shape matrix
phi = V;

% Freq. range to study
freq_mfr = (0:2:2000);

% Dumping ratio
b = 0.02;

% First loop that runs through the freq. range to study
for i = 1:size(freq_mfr, 2)
    % Hz to rad/s
    w_mfr = 2 * pi * freq_mfr(1, i);
    % Second loop to run through all the eigenvalues, the system's
    % frequencies that are considered (11)
    for j = 1:size(V, 2)
        % Force transformation into modal coord
        F_mfr = transpose(phi(:, j)) * F_n;
        % Mass matrix transformation into modal coord
        mass_term = transpose(phi(:, j)) * M(in_n, in_n) * phi(:, j);
        % Stiffness matrix transformation into modal coord
        stiffness_term = transpose(phi(:, j)) * K_nn * phi(:, j);
        % dumping ratio definition
        dump_term = 2 * mass_term * b * 2 * pi * freq(j, 1);
        % displacement in modal coord calculation
        xi(j, 1) = (F_mfr) / (-w_mfr^2*mass_term + stiffness_term + (1i*dump_term*w_mfr));
    end
    % transform modal displ. to physical displ.
    x(:, i) = phi * xi;
end

%% PLOT
ref = (1305 - 1) * 6 + 1;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

fig1 = figure(1);
plot(freq_mfr, abs(x(ref, :)))
xlabel('Frequency [Hz]', 'FontSize', 12)
ylabel('x displacement', 'FontSize', 12)
ylim([0 4.7*10^-3])
grid on 
grid minor
box on
set(fig1,'Units','points');
pos = get(fig1,'Position');
set(fig1,'PaperPositionMode','Auto','PaperUnits','Points','PaperSize',[pos(3), pos(4)])
title_1 = "modal_freq_resp";
print(fig1, title_1,'-dpdf','-r0')

