% Particle in a 1D Infinite Potential Well
% ---------------------------------------
% This script numerically solves the time-independent Schrödinger equation
% for a particle confined in a 1D infinite potential well using a finite-
% difference method. Energy eigenvalues and eigenstates are computed and
% compared to analytic solutions.

%Set Up
L = 1e-9;
N = 1000;
n_show = 5;

% Constants
hbar = 1.054e-34; % joules
m    = 9.109e-31; % kilograms

% Spatial Grid
% Discretize the domain [0, L] into N evenly spaced points
 x = linspace(0,L,N).';
dx = x(2) - x(1);

% Second Derivative Operator
% Construct finite-difference approximation of d^2/dx^2
% using a tridiagonal matrix
main = -2 * ones(N,1);
off  = 1 * ones (N-1,1);
D2   = (diag(main)+diag(off, +1) + diag(off, -1)) /dx^2;

%Apply infinite well boundary conditions
% Enforce psi(0) = psi(L) = 0 by solving only on interior points
D2i = D2 (2:end-1, 2:end-1);
xi  = x(2:end-1);

% H = -(hbar^2 / 2m) * d^2/dx^2 for V(x) = 0 inside the well
H = -(hbar^2)/(2*m) * D2i;

% Solve H psi = E psi
[vecs, vals] = eig(full(H));
[E_all, idx] = sort(diag(vals), 'ascend');
psi_i        = vecs(:, idx);


E = E_all(1:n_show);
psi = psi_i(:, 1:n_show);

% Normalization
% Ensure integral |psi|^2 dx = 1
psi_full = zeros(N, n_show);
for k = 1:n_show
    normk = sqrt(trapz(xi, abs(psi(:,k)).^2));
    phik  = psi(:,k) / normk;
    psi_full(2:end-1, k) = phik;
end
eV  = 1.602176634e-19;
E_eV = E / eV;


n_vec = (1:n_show).';
E_analytic = (n_vec.^2) * (pi^2 * hbar^2 / (2*m*L^2));
E_analytic_eV = E_analytic / eV;


fprintf('n     Numeric E (eV)    Analytic E (eV)   Rel. error\n');
for k = 1:n_show
    relerr = abs(E_eV(k) - E_analytic_eV(k)) / E_analytic_eV(k);
    fprintf('%-2d    %10.6f        %10.6f       %8.3g\n', ...
        k, E_eV(k), E_analytic_eV(k), relerr);
end
figure('Color','w');
subplot(2,1,1);
plot(x*1e9, zeros(size(x)), 'k', 'LineWidth', 1.2);
grid on; xlabel('x (nm)'); ylabel('V(x) (eV)');
title('Infinite Box: V(x)=0 inside (walls enforced by Dirichlet BC)');

subplot(2,1,2); hold on; grid on;
colors = lines(n_show);
for k = 1:n_show
    plot(x*1e9, abs(psi_full(:,k)).^2 + E_eV(k), 'LineWidth', 1.6, ...
         'Color', colors(k,:));
end
xlabel('x (nm)'); ylabel('|{\psi}|^2 (offset) & Energy (eV)');
title('Lowest eigenstates (|{\psi}|^2 curves are vertically offset by their E_n)');
legend(arrayfun(@(k) sprintf('n=%d (%.3f eV)', k, E_eV(k)), 1:n_show, 'uni', 0), ...
       'Location','bestoutside');

