function REBayes(L)
    @rput L;
    R"require(REBayes);
    t_rebayes = system.time(res <- KWDual(L, rep(1,dim(L)[2]), rep(1,dim(L)[1])/dim(L)[1]))[3];
    res$f[res$f < 1e-3] = 0
    x_rebayes = res$f / sum(res$f)"
    @rget x_rebayes;
    @rget t_rebayes;
    return x_rebayes, t_rebayes
end
