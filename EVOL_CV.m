function phi = EVOL_CV(phi0,Img,lambda1,lambda2,mu,nu,pu,timestep,epsilon,numIter)
%   input:
%       Img: input image
%       phi0: level set function to be updated
%       mu: weight for length term
%       nu: weight for area term, default value 0
%       pu：the weight parameter of the penalizing term
%       lambda_1:  weight for c1 fitting term
%       lambda_2:  weight for c2 fitting term
%       timestep: time step
%       epsilon: parameter for computing smooth Heaviside and dirac function
%       numIter: number of iterations
%   output:
%       phi: updated level set function
phi = phi0;
for k=1:numIter
    phi = NeumannBoundCond(phi);
    delta_h = Dirac(phi,epsilon);
    H = Heaviside(phi,epsilon);
    Curv = curvature(phi,1);
    [C1,C2]=binaryfit_CV(Img,H);
    F_CV = (-lambda1*(Img-C1).^2+lambda2*(Img-C2).^2);
    penalizingTerm=pu*(4*del2(phi)-Curv);
    lengthTerm = delta_h.*Curv ;
    phi=phi+timestep*(delta_h.* (F_CV-nu)+ mu.*lengthTerm+penalizingTerm);
end
