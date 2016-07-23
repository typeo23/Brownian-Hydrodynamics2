function [ force ] = CalcForce(parts,eint,rc )
%Calc lennard johns repulsion between parts
N = length(parts);
NN = 2*N;

kb = 1.3806488e-23;
T=298;
kT = kb*T;
Eint = eint*kT;
%rc = .1e-7;

fac = 12*Eint*rc^11;
forceTmp = zeros(NN,N);
for i = 1:N
    for j =i+1:N
        r1 = parts(i,:);
        r2 = parts(j,:);
        rij2 =(r1-r2)*(r1-r2)';% norm(r1-r2);
        rijv = (r1-r2);
        tmp = fac*(rijv./rij2^6);
        forceTmp(2*i-1,j) = tmp(1);
        forceTmp(2*i,j) = tmp(2);
        forceTmp(2*j-1,i) = -tmp(1);
        forceTmp(2*j,i) = -tmp(2);
    end
end
force = sum(forceTmp')';


end

