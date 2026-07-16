% Aliakbar Hoveydapour 
% =========================================================
%  Heat Transfer 2D Finite Difference - Gauss-Seidel Project
% =========================================================
clear; clc;

% --- Material & grid ---
k     = 180;        % W/(m·K)  - thermal conductivity
dx    = 0.01;       % m  (= dy, square grid)
dy    = dx;
sigma = 5.67e-8;    % W/(m²·K⁴)  Stefan-Boltzmann

% --- Boundary conditions ---
% BC1 – Bottom face: constant temperature
T_bottom = 60 + 273.15;       % K

% BC2 – Left face: convection + radiation
T_inf1 = 200 + 273.15;        % K  (hot fluid)
T_sur  = 150 + 273.15;        % K  (radiation surroundings)
eps    = 0.8;                  % emissivity
h1     = 50;                   % W/(m²·K)

% BC3 – Right face of thin base: convection
T_inf2 = 25 + 273.15;         % K
h2     = 20;                   % W/(m²·K)

% BC4 – Top surface of thin base: INSULATED  (q'' = 0)

% BC5 – All other outer surfaces: convection
T_inf3 = 25 + 273.15;         % K
h3     = 25;                   % W/(m²·K)

% --- Initialisation ---
N = 41;
T = ones(N,1) * 350;          % initial guess [K]
T(1:7) = T_bottom;             % fix constant-T nodes

% --- Solver settings ---
max_iter = 50000;
tol      = 1e-6;               % tighter tolerance (K)
err      = Inf;
iter     = 0;

fprintf('Running Gauss-Seidel with linearised radiation \n');

% =========================================================
%  Gauss-Seidel loop
% =========================================================
while err > tol && iter < max_iter

    T_old = T;

    % ----------------------------------------------------------
    %  Compute LINEARISED radiation coefficient at left-face nodes
    %  h_r(T) = eps*sigma*(T_sur^2 + T^2)*(T_sur + T)
    %  so  q_rad = h_r*(T_sur - T)   [fully linear in T now]
    % ----------------------------------------------------------
    hr8  = eps * sigma * (T_sur^2 + T(8)^2)  * (T_sur + T(8));
    hr15 = eps * sigma * (T_sur^2 + T(15)^2) * (T_sur + T(15));
    hr22 = eps * sigma * (T_sur^2 + T(22)^2) * (T_sur + T(22));
    hr26 = eps * sigma * (T_sur^2 + T(26)^2) * (T_sur + T(26));

    % ----------------------------------------------------------
    %  Left boundary – convection + LINEARISED radiation
    %  General interior-edge energy balance :
    % ----------------------------------------------------------
   
    % Node 8  (left edge, interior row)
    T(8)  = ( k*(T(1) + T(15) + 2*T(9))  + (2*dx*(h1*T_inf1 + hr8 *T_sur)) ) ...
            / ( 4*k + 2*dx*(h1 + hr8)  );

    % Node 15 (left edge, interior row)
    T(15) = ( k*(T(8) + T(22) + 2*T(16)) + 2*dx*(h1*T_inf1 + hr15*T_sur) ) ...
            / ( 4*k + 2*dx*(h1 + hr15) );

    % Node 22 (left edge, interior row)
    T(22) = ( k*(T(15)+ T(26) + 2*T(23)) + 2*dx*(h1*T_inf1 + hr22*T_sur) ) ...
            / ( 4*k + 2*dx*(h1 + hr22) );

    % Node 26 (outer corner – left Conv+Rad, top Conv3)
    %   corner balance: k*(Ts + Te) + dx*(h1*T_inf1 + hr*T_sur) + dx*h3*T_inf3
    %   -----------------------------------------------------------------------
    T(26) = ( k*(T(22) + T(27))           + dx*(2*h1*T_inf1 + 2*hr26*T_sur + 2*h3*T_inf3) ) ...
            / ( 2*k + dx*(2*h1 + 2*hr26 + 2*h3) );

    % ----------------------------------------------------------
    %  9,10,16,17,23,24 are Internal nodes 
    % ----------------------------------------------------------
    T(9)  = (T(2)  + T(16) + T(8)  + T(10)) / 4;
    T(10) = (T(3)  + T(17) + T(9)  + T(11)) / 4;
    T(16) = (T(9)  + T(23) + T(15) + T(17)) / 4;
    T(17) = (T(10) + T(24) + T(16) + T(18)) / 4;
    T(23) = (T(16) + T(27) + T(22) + T(24)) / 4;
    T(24) = (T(17) + T(28) + T(23) + T(25)) / 4; 
    % ----------------------------------------------------------
    %  Top-of-base-block boundary nodes (Conv3)
    % ----------------------------------------------------------

    % Node 27 (top edge, interior)
    T(27) = ( k*(2*T(23) + T(26) + 2*T(28) + T(30)) + 2*dx*h3*T_inf3 ) ...
            / ( 6*k + 2*dx*h3 );

    % Node 28 (interface – internal)
    T(28) = (T(24) + T(31) + T(27) + T(29)) / 4;

    % Node 29 (top edge, interior – Conv3)
    T(29) = ( k*(2*T(28) + T(32) + T(25)) + 2*dx*h3*T_inf3 ) ...
            / ( 4*k + 2*dx*h3 );

    % Node 25 (right edge of upper-base region – Conv3)
    T(25) = ( k*(2*T(24) + T(29) + T(18)) + 2*dx*h3*T_inf3 ) ...
            / ( 4*k + 2*dx*h3 );

    % ----------------------------------------------------------
    %  Lower-base internal nodes
    % ----------------------------------------------------------
    T(11) = (T(4)  + T(18) + T(10) + T(12)) / 4;
    T(12) = (T(5)  + T(19) + T(11) + T(13)) / 4;
    T(13) = (T(6)  + T(20) + T(12) + T(14)) / 4;

    % Node 18 – inner corner (top-base / base junction, one insulated side)
    T(18) = ( k*(2*T(17) + T(25) + 2*T(11) + T(19)) + dx*h3*T_inf3 ) ...
            / ( 6*k + dx*h3 );

    % Insulated top surface (nodes 19, 20)
    T(19) = (2*T(12) + T(18) + T(20)) / 4;
    T(20) = (2*T(13) + T(19) + T(21)) / 4;

    % Node 14 – right edge, lower base (Conv2)
    T(14) = ( k*(2*T(13) + T(7) + T(21)) + 2*dx*h2*T_inf2 ) ...
            / ( 4*k + 2*dx*h2 );

    % Node 21 – top-right corner (insulated top, Conv2 right)
    T(21) = ( k*(T(14) + T(20)) + dx*h2*T_inf2 ) ...
            / ( 2*k + dx*h2 );

    % ----------------------------------------------------------
    %  Top block (nodes 30–41, all outer surfaces Conv3)
    % ----------------------------------------------------------

    % Left-edge nodes (Conv3)
    T(30) = ( k*(2*T(31) + T(27) + T(33)) + 2*dx*h3*T_inf3 ) / ( 4*k + 2*dx*h3 );
    T(33) = ( k*(2*T(34) + T(30) + T(36)) + 2*dx*h3*T_inf3 ) / ( 4*k + 2*dx*h3 );
    T(36) = ( k*(2*T(37) + T(33) + T(39)) + 2*dx*h3*T_inf3 ) / ( 4*k + 2*dx*h3 );

    % Right-edge nodes (Conv3)
    T(32) = ( k*(2*T(31) + T(29) + T(35)) + 2*dx*h3*T_inf3 ) / ( 4*k + 2*dx*h3 );
    T(35) = ( k*(2*T(34) + T(32) + T(38)) + 2*dx*h3*T_inf3 ) / ( 4*k + 2*dx*h3 );
    T(38) = ( k*(2*T(37) + T(35) + T(41)) + 2*dx*h3*T_inf3 ) / ( 4*k + 2*dx*h3 );

    % Internal nodes
    T(31) = (T(28) + T(34) + T(30) + T(32)) / 4;
    T(34) = (T(31) + T(37) + T(33) + T(35)) / 4;
    T(37) = (T(34) + T(40) + T(36) + T(38)) / 4;

    % Top-edge node (Conv3)
    T(40) = ( k*(2*T(37) + T(39) + T(41)) + 2*dx*h3*T_inf3 ) / ( 4*k + 2*dx*h3 );

    % Top-corner nodes (Conv3 on two faces)
    T(39) = ( k*(T(36) + T(40)) + 2*dx*h3*T_inf3 ) / ( 2*k + 2*dx*h3 );   % left-top corner
    T(41) = ( k*(T(38) + T(40)) + 2*dx*h3*T_inf3 ) / ( 2*k + 2*dx*h3 );   % right-top corner

    % Re-enforce Dirichlet BC
    T(1:7) = T_bottom;

    % Convergence check 
    err  = max(abs(T - T_old));
    iter = iter + 1;

end

T_C = T - 273.15;   % Convertit to Celsius

if err <= tol
    fprintf('Converged in %d iterations  (residual = %.2e K)\n', iter, err);
else
    fprintf('ERROR: did NOT converge after %d iterations (residual = %.2e K)\n', iter, err);
end

fprintf('\n%-8s %-12s\n', 'Node', 'Temp (°C)');
fprintf('%-8s %-12s\n', '----', '---------');
for i = 1:N
    fprintf('Node %2d:  %8.3f °C\n', i, T_C(i));
end
