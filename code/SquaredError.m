function [f, g] = SquaredError(w,X,y)
  r = X*w-y;
  f = sum(r.^2);
  g = 2*(X.'*r);
