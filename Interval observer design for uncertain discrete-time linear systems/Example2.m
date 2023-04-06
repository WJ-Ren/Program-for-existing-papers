% This script simulates the interval observer design method proposed in
% "Interval observer design for uncertain discrete-time linear systems"
% Zhenhua Wang, 2018, Systems & Control Letters
%
% /* Example 2 */
%
% This script uses semidefinite programming (SDP) technique. Please add
% Yalmip toolbox and appropriate SDP solvers to your MATLAB path or put
% them in the same direction with this script.
%
% Additional Toolbox Needed:  Yalmip
% Additional Solver Needed:   SeDuMi
% Download Yalmip: https://yalmip.github.io/
% Download SeDuMi: https://yalmip.github.io/solver/sedumi/
%
%
% Comparisons will be added in the near future...

% Version:              1.5
% Author:               Weijie Ren
% Contact:              weijie.ren@outlook.com
% Initial modified:     Jul. 26, 2020
% Last modified:        Apr. 06, 2023

clc, clear
% close all

%% Simulation settings
Nk = 61;   % total simulation steps

%% System model
% System matrix
A=[-0.54  0.45  0.36            0             0;
    0.63  0.45  0.18         0.36             0;
    0.09  0.45  0.27         0.09          0.18;
       0     0  0.25  0.25*sqrt(2) -0.25*sqrt(2);
       0     0     0  0.25*sqrt(2) -0.25*sqrt(2)];
% Input matrix
B=[-1; 0; 0; 0; 1];
% Output/Measurement matrix
C=[1 0 0 0 0;
   0 0 0 1 0];

% Dimension
n = size(A,1);  % dimension of the state vector
p = size(B,2);  % ...       of the unknown disturbance
m = size(C,1);  % ...       of the measured output vector

%% Compute interval observer gain
[T,N,L,gamma] = gain(A,C);
% Note: In our function "gain", we minimized the performance level \gamma
% instead of specifying a fixed value as in the paper.

% Or one can use data from the paper as follows:
% T = [ 0.2690  0  0  -0.1062  0
%      -0.3737  1  0  -0.0683  0
%      -0.5762  0  1   0.2108  0
%       0.1032  0  0  -0.0313  0
%       0.7532  0  0  -1.0255  1];
%  N = [ 0.7310   0.1062
%        0.3737   0.0683
%        0.5762  -0.2108
%       -0.1032   1.0313
%       -0.7532   1.0255];
%  L = [-0.2291  -0.0974
%        0.7722   0.3010
%        0.2740   0.0900
%       -0.1230  -0.0790
%       -0.4762  -0.0859];

%% Interval information
% Express a matrix as the difference of two nonnegative matrices
% One can type "help ApAm" into the Command Window for detailed information
% Or see the function ApAm in the same direction of this script
[TB_p,TB_m] = ApAm(T*B);    % abbreviations _p, _m: _plus, _minus
[L_p,L_m] = ApAm(L);
[N_p,N_m] = ApAm(N);

% Interval bounds of disturbances and noises
w_upperBound = 1;               % upper bound of the unknown disturbance w
w_lowerBound = -w_upperBound;   % lower ...
v_upperBound = [0.1; 0.1];      % upper bound of the measurement noise v
v_lowerBound = -v_upperBound;   % lower ...

% Delta's upper & lower bounds, equations (9) & (10), on page 42
Delta_u = TB_p*w_upperBound - TB_m*w_lowerBound - (L_p*v_lowerBound - L_m*v_upperBound)...
          + N_p*v_upperBound - N_m*v_lowerBound;
      
Delta_l = TB_p*w_lowerBound - TB_m*w_upperBound - (L_p*v_upperBound - L_m*v_lowerBound)...
          + N_p*v_lowerBound - N_m*v_upperBound;
% Please note that (9) and (10) in the paper may be INCORRECT (as follows)!
% Delta_u = TB_p*w_upperBound - TB_m*w_lowerBound + L_p*v_upperBound - L_m*v_lowerBound...
%           + N_p*v_upperBound - N_m*v_lowerBound;
%       
% Delta_l = TB_p*w_lowerBound - TB_m*w_upperBound + L_p*v_lowerBound - L_m*v_upperBound...
%           + N_p*v_lowerBound - N_m*v_upperBound;

%% Declare variables
% System variables
x = zeros(n,Nk);        % system state
y = zeros(m,Nk);        % measurable output
w = zeros(p,Nk);        % unknown disturbance
v = zeros(m,Nk);        % measurement noise
% Observer variables
x_u = zeros(n,Nk);      % upper estimation bound 
x_l = zeros(n,Nk);      % lower ...
zeta_u = zeros(n,Nk);   % upper bound of intermediate variable
zeta_l = zeros(n,Nk);   % lower ...

%% Initialization
x(:,1) = [-0.3 -0.5 0.6 0.9 -0.2]'; % initial system state
x_u(:,1) = ones(n,1);               % initial upper estimation bound
x_l(:,1) = -x_u(:,1);               % ...     lower ...

%% main loop
for k=1:Nk
    % Disturbance signal
    w(:,k) = w_upperBound .* (2 * rand(p,1) - 1);
    % Noise signal
    v(:,k) = v_upperBound .* (2 * rand(m,1) - 1);
    
    % System <-- eq. (1)
    x(:,k+1) = A * x(:,k) + B * w(:,k);     % state equation
    y(:,k) = C * x(:,k) + v(:,k);           % output equation
    
    % Interval observer <-- eq. (8)
    if k >= 2
        x_u(:,k) = zeta_u(:,k) + N * y(:,k);
        x_l(:,k) = zeta_l(:,k) + N * y(:,k);
    end
    zeta_u(:,k+1) = T * A * x_u(:,k) + L * (y(:,k) - C * x_u(:,k)) + Delta_u;
    zeta_l(:,k+1) = T * A * x_l(:,k) + L * (y(:,k) - C * x_l(:,k)) + Delta_l;  
end

%% Plot results
Nks = 0:Nk-1;

for i = 1:n
    figure
    plot(Nks,x(i,1:Nk),'k')
    hold on
    plot(Nks,x_u(i,1:Nk),'r-.')
    plot(Nks,x_l(i,1:Nk),'b-.')
    xlabel('Time (samples)')
    ylabel(['x_' num2str(i)])
    title(['Responses of x_' num2str(i) ' and its interval estimation in Example 2'])
end
