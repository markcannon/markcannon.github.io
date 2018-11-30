function xdot = ex4b_ode(t,x)

x1dot = x(1)^2+x(1)*x(2);
x2dot = 0.5*x(2)^2+x(1)*x(2);
xdot = [x1dot;x2dot];