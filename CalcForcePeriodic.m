function [ force ] = CalcForcePeriodic(parts,eint,rc,boxSize)
%Calc lennard johns repulsion between parts; Periodic version
N = round(length(parts)/2);
NN = 2*N;

kb = 1.3806488e-23;
T=298;
kT = kb*T;
Eint = eint*kT;
%rc = 1e-7;

fac = 12*Eint*rc^11;
forceTmp = zeros(NN,N);
for i = 1:N
    for j =i+1:N
        r1 = parts((2*i-1):2*i);
        r20 = parts((2*j-1):2*j);
        for x = -1:1
            for y=-1:1
                
                r2 = r20 + [x,y]'*boxSize;
                rij2 =(r1-r2)'*(r1-r2);% norm(r1-r2);
                rijv = (r1-r2);
                tmp = fac*(rijv./rij2^6);
                forceTmp(2*i-1,j) = forceTmp(2*i-1,j)+tmp(1);
                forceTmp(2*i,j) =forceTmp(2*i,j)+ tmp(2);
                forceTmp(2*j-1,i) =forceTmp(2*j-1,i) -tmp(1);
                forceTmp(2*j,i) =forceTmp(2*j,i)  -tmp(2);
            end
        end
    end
end
force = sum(forceTmp')';


end

