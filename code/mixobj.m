function [f, g] = mixobj (L, x, e)
  p = numel(x);
  y = L*x + e;
  if any(y <= 0)
    f = Inf;
    g = zeros(p,1);  
  else 
    f = -sum(w * log(y));
    % TO DO: Compute g.
  end
