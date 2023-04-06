function [T,N,L,gamma] = gain(A,C)
% GAIN computes interval observer gain L, matrix T, and matrix N, by
% H_{\infty} technique (i.e., bounded real lemma).

% dimensions
n=size(A,1);
m=size(C,1);
zt = 1e-9;

alpha1 = [eye(n); zeros(m,n)];
alpha2 = [zeros(n,m); eye(m)];
Theta = [eye(n);C];
Psi = eye(n+m) - Theta*pinv(Theta);
phi1 = pinv(Theta)*alpha1;
phi2 = pinv(Theta)*alpha2;
eta1 = Psi*alpha1;
eta2 = Psi*alpha2;

P = sdpvar(n,n,'diagonal');
Y = sdpvar(n,n+m,'full');
W = sdpvar(n,m,'full');
gamma2 = sdpvar(1);

psi11 = [eye(n)-P,     zeros(n,n)',     zeros(m,n)',     zeros(m,n)';
         zeros(n,n),   -gamma2*eye(n),  zeros(m,n)',     zeros(m,n)';
         zeros(m,n),   zeros(m,n),      -gamma2*eye(m),  zeros(m,m)';
         zeros(m,n),   zeros(m,n),      zeros(m,m),      -gamma2*eye(m)];

psi21 = [P*phi1*A+Y*eta1*A-W*C,  P,  W,   P*phi2+Y*eta2];

psi22 = -P;

F1 = [psi11, psi21';
      psi21, psi22];

F2 = P*phi1*A+Y*eta1*A-W*C;

constraints = [F1<=0,F2>=0,gamma2>=0,P>=zt,norm(Y,'fro')<=100];

objective = gamma2;
ops = sdpsettings('solver','sedumi','sedumi.eps',1e-12,'verbose',0);
sol = optimize(constraints,objective,ops);

P = double(P);
Y = double(Y);
W = double(W);
gamma = sqrt(double(gamma2));

L = P \ W;
S = P \ Y;
T = phi1 + S*eta1;
N = phi2 + S*eta2;
end

