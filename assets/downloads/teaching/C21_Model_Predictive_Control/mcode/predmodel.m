function [s,p,w,robsv,info] = predmodel(A,B,C,umax,nc,Qcost,Rcost)
%predmodel Constructs augmented prediction model.
%   [s,param,rstate,Rcost] = init_param(attr,rstate,Qcost,Rcost) returns:
%   State space plant model s
%   Augmented prediction model parameters are returned in structure p:
%    p.nx    plant state dimension = length(x)         }
%    p.nu    plant input dimension = length(u)         } z[k] = [x[k];c[k]]
%    p.nc    prediction horizon: p.nu*p.nc = length(c) }
%    p.Phi   prediction system state transition matrix z[k+1] = Phi*z[k]
%    p.K     prediction system feedback gain u[k] = K*z[k]
%    p.umax  input constraint (symmetric) |u[k]| < umax, k = 0,1...
%   State and control penalty matrices in LQ cost for determining p.K
%   are returned in structure w:
%    w.Qcost = Qcost     -- numeric Qcost
%    w.Rcost = Rcost     -- numeric Rcost
%    w.Qcost = s.C'*s.C  -- Qcost = 'out' or ''
%    w.Qcost = eye(p.nx) -- Qcost = 'eye'
%    w.Rcost = eye(p.nu) -- Rcost = ''
%   Infinite Horizon Prediction cost matrices are returned in w:
%    w.P,w.W  -- x'*w.P*x + c'*w.W*c = sum_k x[k]'*Qcost*x[k]+u[k]'*Rcost*u[k]

% Plant model
s = ss(A,B,C,0,1);

[p.nx,p.nu] = size(s.B);
p.nc = nc;
p.umax = umax;

% Cost weights
if nargin < 7 || isempty(Rcost)
  w.Rcost = eye(p.nu);
else
  w.Rcost = Rcost;
end
if nargin < 6 || isempty(Qcost)
  w.Qcost = s.C'*s.C;
elseif ~isnumeric(Qcost)
  if isempty(Qcost) || ~strcmp(Qcost(1),'e')
    w.Qcost = s.C'*s.C;
  else
    w.Qcost = eye(p.nx);
  end
else
  w.Qcost = Qcost;
end

% LQ feedback gain and weighting
[P,lam,K,info] = dare(s.A,s.B,w.Qcost,w.Rcost,'report');
w.P = P;
w.W = s.B'*P*s.B + w.Rcost;

% Prediction system: p.Phi, p.K
if p.nc > 0
  p.K = [-K,eye(p.nu),zeros(p.nu,p.nu*(p.nc-1))];
  ZZ = zeros(p.nx,p.nu*p.nc);
  T = diag(ones(1,p.nu*(p.nc-1)),p.nu);
  p.Phi = [([s.A,ZZ]+s.B*p.K);ZZ',T];
else
  p.K = -K;
  p.Phi = s.A + s.B*p.K;
end
robsv = svd(obsv(p.Phi,p.K)); robsv = robsv(end)/robsv(1);
