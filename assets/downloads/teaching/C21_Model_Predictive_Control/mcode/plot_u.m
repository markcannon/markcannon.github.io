function [tt,uu] = plot_u(t,u,cline)

if nargin == 1
  u = t; t = 0:(size(u,2)-1); 
end
tt = reshape(ones(2,1)*t,1,2*length(t));
tt = [tt(2:end),tt(end)+1];

uu = zeros(size(u,1),2*size(u,2));
for i = 1:size(u,1)
  uu(i,:) = reshape(ones(2,1)*u(i,:),1,2*size(u,2));
end

if nargout == 0
  if(nargin < 3), cline = ''; end
  plot(tt,uu,cline);
end
