function  SimPartsNparticlsNonPeriodic(filename,N,tend,dt,SaveInterval)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulate brownian dynamics of N parts in a membrene
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Eq = false;
kb = 1.3806488e-23;
T=298;
kT = kb*T;
boxSize=6e-6;
Eta_m=1e-9;%.5;%.5e-9;
Eta_f=1e-3;
Kappa = (2*Eta_f)/Eta_m;
% dt = 1e-6;
% tend = 1e-0;
steps = round(tend/dt);
%SaveInterval = 1e2;
countSave=1;
saveFileName = sprintf('SimPartsN_%s_N=%d_tend=%3.3f_dt=%1.1e_%s.mat',filename,N,tend,dt,date);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set up initial position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N=5; %number of particles in the box
NN= 2*N;
sq2dt = sqrt(2*dt);
%parts = (rand(N,2)-0.5)*1e-6;%[0,0];
parts = arrangePartsGrid2(sqrt(N),6e-6);
saveRnd = zeros(round(steps/SaveInterval),2*N);
saveParts = zeros(round(steps/SaveInterval),2*N);
saveVel = saveParts;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run a few timesteps without random force to make sure the particles do
% not overlap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (Eq)
    for i = 1:round(steps/1e5)
        F = CalcForce_mex(parts);
        parts = parts + reshape(F,2,N)'*(dt/1e2)*1e14;
        parts(parts(:,1) > boxSize/2,1) = -mod(parts(parts(:,1) > boxSize/2,1),boxSize/2) ;
        parts(parts(:,2) > boxSize/2,2) = -mod(parts(parts(:,2) > boxSize/2,2),boxSize/2) ;
        parts(parts(:,1) < -boxSize/2,1) = mod(parts(parts(:,1) < -boxSize/2,1),boxSize/2) ;
        parts(parts(:,2) < -boxSize/2,2) = mod(parts(parts(:,2) < -boxSize/2,2),boxSize/2) ;
%         scatter(parts(:,1),parts(:,2));
%         axis([-boxSize/2 boxSize/2 -boxSize/2 boxSize/2]);
%         drawnow();
         fprintf('%03.3f %% done\r',((i/round(steps/1e5))*100));
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% running the simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=0;
 D = rand(2*N,2*N)*3;
    D = D*D';
for i = 1:steps
   % F = CalcForcePeriodic_mex(parts,boxSize); %calculate lenard johns repulsion
    
    %if the time step is small enough we dont need to calculate the
    %diffusion matrix ebery timestep
   
    if (mod(i-1,10)==0)
        D = BuildDiffusionMatrix2Periodic_mex(parts,Eta_m,Kappa,kT,boxSize);
    end
        %Checks if the matrix is potive definite, Oseen tensoe may not be for
        %distances less then the particle raduis, thou this should not happen
        %because of the repulsive interaction
        try
            %L= chol(D,'lower');
            ranF=mvnrnd(zeros(2*N,1),D)*sq2dt;
        catch
            fprintf( 'zibi\n');
          
            
        end
    
    
    %random force on the particles
    %ranF=(F'*D/kT)'*dt+(L*randn(NN,1).*sq2dt);
    %ranF= F'*D/kT + ranF;
    %taking care of periodic boudary conditions
    parts = parts + ranF';
    parts(parts > boxSize/2) = parts(parts > boxSize/2) - boxSize;
    parts(parts < -boxSize/2) = parts(parts < -boxSize/2) + boxSize;
    
    
%     parts(parts(:,1) > boxSize/2,1) = parts(parts(:,1) > boxSize/2,1) - boxSize;
%     parts(parts(:,2) > boxSize/2,2) = parts(parts(:,2) > boxSize/2,2) - boxSize ;
%     parts(parts(:,1) < -boxSize/2,1)= parts(parts(:,1) < -boxSize/2,1)+ boxSize;
%     parts(parts(:,2) < -boxSize/2,2) = parts(parts(:,2) < -boxSize/2,2) + boxSize ;
    
    
    
    if (mod(i,SaveInterval) == 0)
        
        %plotting stuff every few intervals
        scatter(parts(1:2:2*N),parts(2:2:2*N));
        axis([-boxSize/2 boxSize/2 -boxSize/2 boxSize/2]);
        drawnow();
       
        saveParts(countSave,:) =parts;% parts(:,1);
        saveVel(countSave,:)= ranF/dt;
        %saveParts(countSave,:,2) = ranF(2:2:end)';%parts(:,2);
        saveRnd(countSave,:)=ranF;
        
            
       if (mod(countSave,100) ==0)
        save(char(saveFileName));
        fprintf('%03.3f %% done\r',((i/steps)*100));
       end
        countSave = countSave+1;
        
    end
end
D;
end


