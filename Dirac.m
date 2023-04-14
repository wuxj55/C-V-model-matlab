function delta_h = Dirac(phi, epsilon)
%  compute the smooth Dirac function
delta_h=(epsilon/pi)./(epsilon^2+ phi.^2);