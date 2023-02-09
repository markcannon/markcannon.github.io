function c = linconstr_d(p,Bd,N)
%linconstr Computes maximal admissible set for the augmented prediction
%   system state z corresponding to input predictions:
%    u[k]= p.K*z[k]u[k] = , z[k+1] = p.Phi*z[k]
%   subject to limits: p.umin <= u[k] <= p.umax
%   [c] = linconstr(p,N) returns constraint matrices in structure c:
%    c.A,c.Bx,c.b -- c.A*c <= c.Bx*x + c.b + c.Bd*d
%   Input arguments
%   p -- structure containing Augmented prediction model parameters
%   N -- number of samples between termination check (default = 5)

if(nargin < 3),
  N = 5;
end;

nz = p.nx + p.nu*p.nc;   % nz = dimension of z
umax = p.umax;
umin = -p.umax;          % symmetric constraints assumed

% Linear constraints
A = []; Bx = []; bmax = []; bmin = []; K_loc = p.K;
flag = 0; k = 0;
options = optimset('display','off'); zbnd = 1e3*ones(nz,1);
while ~flag
  for i=1:N
    A = [A;K_loc(:,p.nx+1:end)];
    Bx = [Bx;-K_loc(:,1:p.nx)];  % u[k] = K_loc*z[k]
    bmax = [bmax;umax]; bmin = [bmin;umin];
    K_loc = K_loc*p.Phi; % u[k+1] = K_loc*z[k]
  end
  j = 1; flag = 1;
  while (flag && j <= 2*p.nu)
    jj = mod(j-1,p.nu)+1; sgn = sign(j-(p.nu+0.5));
    [z,ujmax] = ...
      linprog(sgn*K_loc(jj,:),[[-Bx,A];[Bx,-A]],[bmax;-bmin], ...
		[],[],-zbnd,zbnd,[],options);
    ujmax = sgn*ujmax;
    fprintf(1,'%d %d %f\n',k,j,ujmax);
    if ( (sgn < 0 && ujmax < umax(jj)) || (sgn > 0 && ujmax > umin(jj)) )
      j = j+1; % check max/min of next element of u
    else
      flag = 0; % need to add hyperplane K_loc*z = umax or umin
                % to admissible set boundary
    end
  end
  k = k+1;
end
c.A = [A;-A];
c.Bx = [Bx;-Bx];
c.b = [bmax;-bmin];

nd = size(Bd,2);
BBd = Bx*Bd;
BBd = [zeros(p.nu,nd);BBd(1:end-p.nu,:)];
BBd = cumsum(BBd);
c.Bd = [BBd;-BBd];