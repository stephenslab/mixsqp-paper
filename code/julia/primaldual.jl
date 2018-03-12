# vanilla interior point method for the simplex
function ip_simp(L)
    n,m = size(L);
    mod = Model(solver=MosekSolver(QUIET = true));
    @variable(mod, x[1:m] >= 0, start = 1/m);
    @NLobjective(mod, Min, -sum(log(sum(L[i,j]*x[j] for j = 1:m)) for i = 1:n)/n);
    @constraint(mod, sum(x) == 1)
    solve(mod);
    x = getvalue(x);
    x[x.< 1e-3] = 0;
    return sparse(x/sum(x));
end

# vanilla interior point method for the box constrained
function ip_box(L)
    n,m = size(L);
    mod = Model(solver=MosekSolver(QUIET = true));
    @variable(mod, x[1:m] >= 0, start = 1/m);
    @NLobjective(mod, Min, -sum(log(sum(L[i,j]*x[j] for j = 1:m)) for i = 1:n)/n + sum(x[i] for i = 1:m));
    solve(mod);
    x = getvalue(x);
    x[x.< 1e-3] = 0;
    return sparse(x/sum(x));
end

# vanilla interior point method for the dual
function ip_dual(L)
    n,m = size(L);
    mod = Model(solver=MosekSolver(QUIET = true));
    @variable(mod, y[1:n] >= 0, start = 1/n);
    @NLobjective(mod, Min, -sum(log(y[i]) for i = 1:n)/n);
    @constraint(mod, ic, L'*y .<= 1)
    solve(mod);
    x = -getdual(ic);
    x[x.< 1e-3] = 0;
    return sparse(x/sum(x));
end

function sqp_box(L)
    n,m = size(L);
    x = ones(m)/m;
    for i = 1:100
        D = 1./(L*x + 1e-8);
        g = -L'*D/n;
        H = L'*Diagonal(D.^2)*L/n + 1e-8 * speye(m);
        
        if minimum(g+1) >= -1e-6
          break;
        end
        
        mod = Model(solver=MosekSolver(QUIET = true));
        a,b = findn(H);
        @variable(mod, y[1:m] >= 0);
        @objective(mod, Min, QuadExpr(y[a],y[b],H[:]/2,AffExpr(y, 2*g+1, 0)) )
        solve(mod);
        x = getvalue(y);
        x[x .< 0] = 0;
        
    end
    x[x .< 1e-3] = 0;
    return x
end

function sqp_simp(L)
    n,m = size(L);
    x = ones(m)/m;
    for i = 1:100
        D = 1./(L*x + 1e-8);
        g = -L'*D/n;
        H = L'*Diagonal(D.^2)*L/n + 1e-8 * speye(m);
        
        if minimum(g+1) >= -1e-8
          break;
        end
        
        mod = Model(solver=MosekSolver(QUIET = true));
        a,b = findn(H);
        @variable(mod, y[1:m] >= 0);
        @objective(mod, Min, QuadExpr(y[a],y[b],H[:]/2,AffExpr(y, 2*g, 0)) )
        @constraint(mod, sum(y) == 1);
        solve(mod);
        x = getvalue(y);
        x[x .< 0] = 0;

    end
    x[x .< 1e-3] = 0;
    return x
end

function sqp_dual(L)
    n,m = size(L);
    x = ones(n)/n; lambda = zeros(m);
    for i = 1:100
        D = 1./x;
        
        mod = Model(solver=MosekSolver(QUIET = true));
        @variable(mod, y[1:n] >= 0);
        @objective(mod, Min, QuadExpr(y,y,D.^2/2/n,AffExpr(y, -2*D/n, 0)) )
        @constraint(mod, ic, L'*y .<= 1);
        solve(mod);
        x = getvalue(y);
        lambda = -getdual(ic);
        x[x .< 0] = 0;
        if minimum(L *lambda - x) > 0
            break;
        end
    end
    lambda[lambda .< 1e-3] = 0;
    lambda = lambda/sum(lambda);
    
    return lambda
end