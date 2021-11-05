%% Dirichelt random variable
function out = sy_dir(a)
    out = random('gam', a, 1);
    out = out / sum(out);
end
