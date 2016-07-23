function D = BuildDiffusionMatrix2(parts,Eta_m,Kappa,kT,boxSize)
% coder.extrinsic('bessely');
% coder.extrinsic('besselj');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Build the full diffusion matrix for all of the paticles
%@in parts : Nx2 - particle position array
%@in Eta_m membrene viscosity
%@Kappa SD length
%@in kT = tempertaure (times kb)
%@out D= 2Nx2N diffusion matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(parts,1)/2;
D = zeros(2*N,2*N);

for i=1:2:2*N
    for j=(i+2):2:2*N
        Dtmp = BuildDiffusionMatrixInteraction(parts(i:i+1),parts(j:j+1),Eta_m,Kappa,kT,boxSize);
        D(i:(i+1),j:(j+1)) = Dtmp;
        D(j:(j+1),i:(i+1)) = Dtmp;
    end
end
for i = 1:2:2*N
    D(i:i+1,i:i+1) = BuildDiffusionMatrixSelf(Eta_m,Kappa,kT);
end
end

function [D] = BuildDiffusionMatrixInteraction(r1,r2,Eta_m,Kappa,kT,boxSize)
%coder.extrinsic('bessely');
% Build the 3x3 diffusion interaction matrix for the two particles
D = zeros(2,2);


r = r1-r2; 
rn = norm(r);
if (rn < 0.5e-7)
    'bukita'
end
kr = Kappa*rn;
if (kr > 1)
    return
end
%H0 = struveH0AproxSum(kr,10);%struve_mex(0,kr,cast(20,'int32'));%interp1q(k_r,H0_int',kr);
%H0 = StruveH0(kr);%struve_mex(0,kr,cast(20,'int32'));
%H0=H0look(kr);
if coder.target('MATLAB')
    H0 = StruveH0(kr);%struveH0Aprox4_mex(kr);
    H1 = StruveH1(kr);%real(struveH1Aprx(kr));%struve_mex(1,kr,cast(20,'int32'));%interp1q(k_r,H1_int',kr);
    Y0 = real(bessely(0,kr));
    Y2 = real(bessely(2,kr));
    %display(H0);
else
    H0=0.0;
    H1=0.0;
    Y0=0.0;
    Y2=0.0;
    coder.ceval('STVH0',kr,coder.wref(H0));
    coder.ceval('STVH1',kr,coder.wref(H1));
    coder.ceval('BESSY',cast(0,'int64'),kr,coder.wref(Y0));
    coder.ceval('BESSY',cast(2,'int64'),kr,coder.wref(Y2));
    %fprintf('%f\n',H0);
end
for i = 1:2
    for j = 1:2
        D(i,j) = kT/(4*Eta_m)*(...
            (H0-H1/kr - 0.5*(Y0-Y2) + 2/(pi()*kr^2))*(i==j)...
            -(H0-2*H1/kr + Y2+4/(pi()*kr^2))*(r(i)*r(j)/rn^2)...
            ); %Haim
    end
end
end

% function [ out ] = struveH0Aprox( z )
% %UNTITLED2 Summary of this function goes here
% %   Detailed explanation goes here
% %coder.extrinsic('besselj');
% besfun =@(x)(-1)^x*besselj(2*x+1,z).*sin(2*x+1)*(pi()/2)./(2*x+1);
% tmp = 1:10;
% % for i = 0:11
% %     tmpi=0;
% %     tmpi = real(besselj(2*i+1,z,1));
% %     tmp = tmp+((-1)^i * tmpi) .*(sin((2*i+1)*(pi()/2))./(2*i+1));
% % end
% out = (4/pi())*sum(arrayfun(besfun,tmp));
% end


function [D] = BuildDiffusionMatrixSelf(Eta_m,Kappa,kT)
% Build the 2x2 self diffusion matrix for the two particles
a=0.5e-9;
D = eye(2,2)*kT/(4*pi()*Eta_m)*(-log((Kappa*a)/2)-0.58);%*2.3e-12; %Yaels diffususion coef;
end

