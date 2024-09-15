function EstMdl = estimate_var(y, p)
    Mdl = varm(size(y,2), p);
    EstMdl = estimate(Mdl, y);
end