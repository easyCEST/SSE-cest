function y_fit = lorentz_N(par, x)
% Lorentzian fitting model
denum = 1 + 4 * ((x - par(1)) ./ par(3)).^2;
y_fit = par(4) - par(2) ./ denum;
end
