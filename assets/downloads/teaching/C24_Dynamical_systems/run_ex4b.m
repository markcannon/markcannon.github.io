clear; close all;

[x1,x2] = meshgrid([-2:0.2:2,-2:0.2:2]);
x1dot = x1.^2+x1.*x2;
x2dot = 0.5*x2.^2+x1.*x2;
quiver(x1,x2,x1dot,x2dot);

[~,x] = ode23(@(t,x) ex4bode(t,x),[0 50],[-0.01,0.2]);
hold on
plot(x(:,1),x(:,2))
[~,x] = ode23(@(t,x) ex4bode(t,x),[0 60],[-0.01,0.3]);
plot(x(:,1),x(:,2))
[~,x] = ode23(@(t,x) ex4bode(t,x),[0 60],[-0.01,0.1]);
plot(x(:,1),x(:,2))
[~,x] = ode23(@(t,x) ex4bode(t,x),[0 60],[0.25,-0.1]);
plot(x(:,1),x(:,2))
[~,x] = ode23(@(t,x) ex4bode(t,x),[0 60],[0.1,-0.1]);
plot(x(:,1),x(:,2))
[~,x] = ode23(@(t,x) ex4bode(t,x),[0 60],[0.2,-0.1]);
plot(x(:,1),x(:,2))
[~,x] = ode23(@(t,x) ex4bode(t,x),[0 60],[0.01,0.01]);
plot(x(:,1),x(:,2))
[~,x] = ode23(@(t,x) ex4bode(t,x),[0 60],[-2,-2]);
plot(x(:,1),x(:,2))
axis([-2,2,-2,2])