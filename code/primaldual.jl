# Compute maximum-likelihood estimates of the mixture weights by
# solving the (primal) simplex-constrained optimization problem with
# MOSEK.
function simplexIP(L)

  # Get the number of rows (n) and columns (m) of the likelihood matrix.
  n, m = size(L);

  # Define the constrained optimization problem in JuMP, and solve it
  # using MOSEK.
  mod = Model(solver = MosekSolver(QUIET = true));
  @variable(mod,x[1:m] >= 0,start = 1/m);
  @NLobjective(mod,Min,-sum(log(sum(L[i,j]*x[j] for j = 1:m)) for i = 1:n)/n);
  @constraint(mod,sum(x) == 1)
  solve(mod);

  # Output the solution.
  return getvalue(x)
end

# Compute maximum-likelihood estimates of the mixture weights by
# solving the reformulated, non-negatively-constrained optimization
# problem with MOSEK.
function nonnegIP(L)

  # Get the number of rows (n) and columns (m) of the likelihood matrix.
  n, m = size(L);

  # Define the constrained optimization problem in JuMP, and solve it
  # using MOSEK.
  mod = Model(solver = MosekSolver(QUIET = true));
  @variable(mod,x[1:m] >= 0,start = 1/m);
  @NLobjective(mod,Min,-sum(log(sum(L[i,j]*x[j] for j = 1:m)) for i = 1:n)/n +
                         sum(x[i] for i = 1:m));
  solve(mod);

  # Output the solution.
  return getvalue(x)
end

# Compute maximum-likelihood estimates of the mixture weights by
# solving the dual problem.
function dualIP(L)

  # Get the number of rows (n) and columns (m) of the likelihood matrix.
  n, m = size(L);

  # Define the constrained optimization problem in JuMP, and solve it
  # using MOSEK.
  mod = Model(solver = MosekSolver(QUIET = true));
  @variable(mod,y[1:n] >= 0,start = 1/n);
  @NLobjective(mod,Min,-sum(log(y[i]) for i = 1:n)/n);
  @constraint(mod,ic,L'*y .<= 1)
  solve(mod);

  # Output the solution to the primal problem (contained in the dual
  # variables).
  return -getdual(ic)
end

# Compute maximum-likelihood estimates of the mixture weights by
# solving the (primal) simplex-constrained optimization problem using
# an SQP method.
function simplexSQP(L; maxiter = 1000, convtol = 1e-8, eps = 1e-8)

  # Get the number of rows (n) and columns (m) of the likelihood matrix.
  n, m = size(L);

  # This is the initial estimate of the solution.
  x = ones(m)/m;

  # Repeat until we reach the maximum number of iterations, or until
  # convergence is reached.
  for i = 1:maxiter

    # Compute ghe gradient and Hessian.
    D = 1./(L*x + eps);
    g = -L'*D/n;
    H = L'*Diagonal(D.^2)*L/n + eps*eye(m);
    
    # Check convergence.
    if minimum(g + 1) >= -convtol
      break;
    end
    
    # Define the constrained quadratic program in JuMP, and solve it
    # using MOSEK.
    mod = Model(solver = MosekSolver(QUIET = true));
    a, b = findn(H);
    @variable(mod,y[1:m] >= 0);
    @objective(mod,Min,QuadExpr(y[a],y[b],H[:]/2,AffExpr(y,2*g,0)))
    @constraint(mod,sum(y) == 1);
    solve(mod);
    x = getvalue(y);
    x[x .< 0] = 0;
  end
    
  # Output the solution.
  return x
end

# Compute maximum-likelihood estimates of the mixture weights by
# solving the non-negatively constrained optimization problem.
function nonnegSQP(L; maxiter = 1000, convtol = 1e-8, eps = 1e-8)

  # Get the number of rows (n) and columns (m) of the likelihood matrix.
  n, m = size(L);

  # This is the initial estimate of the solution.
  x = ones(m)/m;

  # Repeat until we reach the maximum number of iterations, or until
  # convergence is reached.
  for i = 1:maxiter

    # Compute ghe gradient and Hessian.
    D = 1./(L*x + eps);
    g = -L'*D/n;
    H = L'*Diagonal(D.^2)*L/n + eps*eye(m);

    # Check convergence.
    if minimum(g + 1) >= -convtol
      break;
    end

    # Define the constrained quadratic program in JuMP, and solve it
    # using MOSEK.
    mod = Model(solver = MosekSolver(QUIET = true));
    a, b = findn(H);
    @variable(mod,y[1:m] >= 0);
    @objective(mod,Min,QuadExpr(y[a],y[b],H[:]/2,AffExpr(y,2*g + 1,0)));
    solve(mod);
    x = getvalue(y);
    x[x .< 0] = 0;
  end

  # Output the solution.
  return x
end

function sqp_dual(L; maxiter = 1000, verbose = false)
  n, m = size(L);
  x = ones(n)/n;
  z = zeros(m);
  if verbose
    @printf "%4d %0.2e\n"
  end
  for i = 1:maxiter
    D = 1./x;
    mod = Model(solver = MosekSolver(QUIET = true));
    @variable(mod,y[1:n] >= 0);
    @objective(mod,Min,QuadExpr(y,y,D.^2/2/n,AffExpr(y,-2*D/n,0)))
    @constraint(mod,ic,L'*y .<= 1);
    solve(mod);
    x = getvalue(y);
    z = -getdual(ic);
    x[x .< 0] = 0;
    if verbose
      @printf "\n"
    end
    if minimum(L*z - x) > 0
      break
    end
  end
  return z
end
