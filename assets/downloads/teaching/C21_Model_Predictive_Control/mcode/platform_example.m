%% Script: runs platform simulation example

g = 9.81; M = 10; Kv = 7; m = 0.5; umax = 1;

T = 0.1; lam = 1e-4;
A = [1,T;0,1]; B = (Kv/(2*M))*[T^2/2;T]; Bd = -(g/(2*M))*[T^2/2;T];
C = [1,0];

s = ss(A,B,C,0,T);
K = -dlqr(A,B,C'*C,lam);

%% Compute closed loop response with unconstrained optimal controller
%  No integral action 
if ~exist('plotflag','var'), plotflag = 1; end;
if plotflag
  [y,x] = dstep(A+B*K,Bd*m,C,0);
  x = x'; u = K*x; t = 0:(length(y)-1);
  figure;
  subplot(2,2,1); plot(t*T,y,'-o'); ax1 = axis;
  axis([0,1,ax1(3:4)]);
  ylabel('e (metres)');
  title('No integral action');
  subplot(2,2,3); plot(t*T,u,'o'); hold on;
  ylabel('u (Volts)');
  xlabel('t (seconds)');
  [tt,uu] = plot_u(t,u); plot(tt*T,uu,'-');
  ax2 = axis; axis([0,1,ax2(3:4)]);
  hold off;
end

%% Compute closed loop response with unconstrained optimal controller
%  Integral action included
AA = [A,zeros(2,1); C, 1]; BB = [B;0]; BBd = [Bd;0]; CC = [C,0];
sa = ss(AA,BB,CC,0,T);
KK = -dlqr(AA,BB,blkdiag(C'*C,1),lam);

if plotflag
  [y,x] = dstep(AA+BB*KK,BBd*m,CC,0);
  x = x'; u = KK*x; t = 0:(length(y)-1);
  %figure;
  subplot(2,2,2); plot(t*T,y,'-o'); ax1 = axis;
  axis([0,1,ax1(3:4)]);
  ylabel('e (metres)');
  title('With integral action');
  subplot(2,2,4); plot(t*T,u,'o'); hold on;
  [tt,uu] = plot_u(t,u); plot(tt*T,uu,'-');
  ax2 = axis; axis([0,1,ax2(3:4)]);
  ylabel('u (Volts)');
  xlabel('t (seconds)');
  hold off;
end

%% Compute closed loop response with MPC
%  Integral action included
[ssa,pa,wa] = predmodel(AA,BB,CC,umax,5,diag([1,0,1]),1e-4);
c = linconstr_d(pa,BBd,5);
N = length(c.b)/2;

% Plot tightened input constraints (just for illustration)
if plotflag
  figure;
  [tt,uu] = plot_u(0:N-1,(c.Bd(1:N,:)*m+c.b(1:N))');
  plot(tt,uu,'r--'); hold on;
  plot(0:N-1,ones(N,1),'r-.')
  [tt,uu]=plot_u(0:N-1,(-c.Bd(N+1:end,:)*m-c.b(N+1:end))');
  plot(tt,uu,'b--');
  plot(0:N-1,-ones(N,1),'b-.')
  [tt,uu] = plot_u(0:N-1,(min(c.Bd(1:N,:)*m+c.b(1:N),1))');
  plot(tt,uu,'r');
  [tt,uu]=plot_u(0:N-1,(max(-c.Bd(N+1:end,:)*m-c.b(N+1:end),-1))');
  plot(tt,uu,'b');
  xlabel('prediction time'); ylabel('u')
end

% Simulate closed loop system with MPC
[ta,za,ua,ya,Ja,Jruna,info] = sim_qpmin_d([0;0;0],BBd,m,[0,m],20,ssa,pa, ...
                                          wa,c,1);
if plotflag
  figure; N = min(15,length(ta));
  subplot(2,1,1);
  plot(T*ta(1:N),ya(1:N),'b-o'); hold on;
  ylabel('e (metres)');
  subplot(2,1,2);
  [tt,uu] = plot_u(ta(1:N),ua(1:N));
  plot(T*tt,uu,'b-'); hold on;
  plot(T*ta(1:N),ua(1:N),'bo'); hold on;
  plot([0,N*T],[1,1],'--');
  plot([0,N*T],-[1,1],'--');
  ylabel('u (Volts)');
  xlabel('t (seconds)');
end

% % Simulate closed loop system with MPC
% [ta0,za0,ua0,ya0,Ja0,Jruna0,info] = sim_qpmin_d([0;0;0],BBd,0,[0,m],20,ssa,pa, ...
%                                           wa,c,1);
% if plotflag
%   figure; N = min(20,length(ta0));
%   subplot(2,1,1);
%   plot(T*ta0(1:N),ya0(1:N),'b-o'); hold on;
%   ylabel('e (metres)');
%   subplot(2,1,2);
%   [tt,uu] = plot_u(ta0(1:N),ua0(1:N));
%   plot(T*tt,uu,'b-'); hold on;
%   plot(T*ta0(1:N),ua0(1:N),'bo'); hold on;
%   plot([0,N*T],[1,1],'--');
%   plot([0,N*T],-[1,1],'--');
%   ylabel('u (Volts)');
%   xlabel('t (seconds)');
% end