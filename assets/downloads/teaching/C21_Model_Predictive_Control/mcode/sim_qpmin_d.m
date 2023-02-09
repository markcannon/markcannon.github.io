function [t,z,u,y,J,Jrun,info] = ...
    sim_qpmin_d(x0,Bd,d,dbnd,N,s,p,w,c,opt_flag,options)
%sim_qpmin Simulate closed-loop response for QP-based control law.
%   [t,z,u,y,J,info] = sim_qpmin(x0,s,p,w,c,options)
%   Input arguments:
%    x0 -- initial plant state
%    s -- plant state space model
%    p -- structure containing prediction model parameters:
%         p.nx, p.nu, p.nc, p.Phi, p.K, p.umax (see predmodel)
%    w -- structure containing cost matrices
%         w.Qcost, w.Rcost, w.P, w.W (see predmodel)
%         w.Pf (finite horizon cost, see fh_cost)
%    c -- structure containing linear constraint matrices
%         c.A,c.Bx,c.b (see linconstr)
%    opt_flag -- 0 => finite horizon cost (u=0 assumed in mode2 predictions)
%                1 => infinite horizon cost
%    options -- set options.Display = 'off' to supress messages
%   Returns:
%    t,z,u,y -- sample times and responses of state, inputs and outputs
%    J -- infinite horizon closed-loop cost
%    info -- optimization status

if nargin < 11
  options = optimset('quadprog');
  %options.Display = 'off';
  options.LargeScale = 'off';
end 
if (nargin < 7 || isempty(opt_flag)), opt_flag = 1; end;

% objective for QP objective
if opt_flag
  H = [];
  for i = 1:p.nc 
    H = blkdiag(H,w.W);
  end
  H = 0.5*(H+H');
  G = zeros(p.nu*p.nc,p.nx);
  F = w.P;
else
  H = w.Pf(p.nx+1:end,p.nx+1:end);
  G = w.Pf(p.nx+1:end,1:p.nx);
  F = w.Pf(1:p.nx,1:p.nx);
end

if size(dbnd,2) > 0
  AA = []; BBx = []; bb = [];
  for i = 1:size(dbnd,2)
    bb = [bb;c.b+c.Bd*dbnd(:,i)];
    AA = [AA;c.A]; BBx = [BBx;c.Bx];
  end
else
  AA = c.A; BBx = c.Bx; bb = c.b;
end

x_k = x0; Jrun = 0; flag = 1; k = 1; z = []; u = zeros(p.nu,0); y = []; J = [];
while flag && k <= N+1
  [c_k,obj,eflag1,out,lam] = quadprog(H,G*x_k,...
                                      AA,BBx*x_k+bb,...
                                      [],[],[],[],[],options);
  if max(lam.ineqlin) == 0, % unconstrained optimal is feasible
    [N,x_lq,u_lq,y_lq,J_lq] = sim_lq(x_k,Bd,d,Jrun,...
                                     s,p,w,F,-p.umax,p.umax,N-k+1);
    z = [z,[x_lq;zeros(p.nc,N)]]; u = [u,u_lq]; y = [y,y_lq];
    J = [J,J_lq];
    Jrun = Jrun + x_k'*w.P*x_k;
    info(:,k) = eflag1;
    flag = 0;
  elseif eflag1 == -1
    Jrun = -1;
    flag = 0;
  else

%    u_k = p.K*[x_k;c_k];
    u_k = satu(p.K*[x_k;c_k],-p.umax,p.umax);
    Jrun = Jrun + x_k'*w.Qcost*x_k + u_k'*w.Rcost*u_k;
    Jpred = 2*obj + x_k'*F*x_k;
  
    z(:,k) = [x_k;c_k];
    u(:,k) = u_k;
    y(:,k) = s.C*x_k;
    J(:,k) = [Jpred;Jrun];
    info(:,k) = eflag1;

    x_k = s.A*x_k + s.B*u_k + Bd*d;
    k = k+1;
    
    
  end
end
t = 0:(size(z,2)-1);

%------------------------------------------------------------------------------
function [N,x,u,y,J] = sim_lq(x0,Bd,d,Jrun,s,p,w,F,umin,umax,NN)

% Closed-loop response under LQ feedback
x_k = x0; r0 = 1e-3*norm(x_k,2); k = 1;
while (norm(x_k,2) >= r0 || k <= NN) && k < 100
  u_k = satu(p.K(:,1:p.nx)*x_k,umin,umax);
%  u_k = p.K(:,1:p.nx)*x_k;
  Jrun = Jrun + x_k'*w.Qcost*x_k + u_k'*w.Rcost*u_k;
  Jpred = x_k'*F*x_k;
  x(:,k) = x_k;
  u(:,k) = u_k;
  y(:,k) = s.C*x_k;
  J(:,k) = [Jpred;Jrun];
  x_k = s.A*x_k + s.B*u_k + Bd*d;
  k = k+1;
end
N = k-1;

%------------------------------------------------------------------------------
function [x,u,y] = sim_pred(x0,c0,s,p,N,opt_flag)

if nargin < 6 || isempty(opt_flag), opt_flag = 1; end;
if nargin < 5 || isempty(N) || N < p.nc, N = p.nc; end;

% Predicted response
x_k = x0; z_k = [x0;c0]; k = 1;
while k <= N
  if (k <= p.nc || opt_flag)
    u_k = p.K*z_k;
  else % FH mode 2
    u_k = 0;
  end
  x(:,k) = x_k;
  u(:,k) = u_k;
  y(:,k) = s.C*x_k;
  if (k <= p.nc || opt_flag)
    z_k = p.Phi*z_k;
    x_k = z_k(1:p.nx);
  else % FH mode 2
    x_k = s.A*x_k;
  end
  k = k+1;
end
x(:,k) = x_k;
y(:,k) = s.C*x_k;

%------------------------------------------------------------------------------
function u = satu(u,umin,umax)

if u > umax
  u = umax;
elseif u < umin
  u = umin;
end