function [ velocity ] = ShouldItMove( t, positions, charge, center, ProbScat, ProbMoveWhenHappy, velocity, mass)
%ShouldItMove Summary of this function goes here
%   Detailed explanation goes here

k = physconst('Boltzmann'); %Use of constants in matlab
T = 300; % temperature in Kalvin
mass = 9.109E-31; %in kg
v_th = sqrt(k*T/mass);


for n=1:length(positions)
    prob = rand();
    if ProbScat > prob
        %rethermalize the particle's velocity by assigning new Vx
        %and Vy from the MB distribution
        RandAngle = rand(1).*2.*pi;
        RandVelX = rand(1).*v_th.* cos(RandAngle);
        RandVelY = rand(1).*v_th.* sin(RandAngle);
        velocity(n, :) = [RandVelX, RandVelY];
        scatterTime(n,1) = 0;
        if charge < 0 && positions(n, 1)< center && ProbScat < prob
            RandAngle = rand(1).*2.*pi;
            RandVelX = rand(1).*v_th.* cos(RandAngle);
            RandVelY = rand(1).*v_th.* sin(RandAngle);
            velocity(n, :) = [RandVelX, RandVelY];
            scatterTime(n,1) = 0;
            
        elseif charge<0 && positions(n, 1)>= center && ProbMoveWhenHappy > prob
            RandAngle = rand(1).*2.*pi;
            RandVelX = rand(1).*v_th.* cos(RandAngle);
            RandVelY = rand(1).*v_th.* sin(RandAngle);
            velocity(n, :) = [RandVelX, RandVelY];
            scatterTime(n,1) = 0;
            TempPos(n, 1) = positions(n, 1) + velocity(n, 1)*t;
            EnergyParticle = k*T+0.5*mass*(velocity(1)^2 +velocity(2)^2);
            
            
            if TempPos(n, 1) < center && EnergyParticle < 2%required energy of partciel is achieved
                velocity(n, :) = [0, 0];
            end
            
        elseif charge > 0 && positions(n, 1)> center && ProbScat < prob
            
            RandAngle = rand(1).*2.*pi;
            RandVelX = rand(1).*v_th.* cos(RandAngle);
            RandVelY = rand(1).*v_th.* sin(RandAngle);
            velocity(n, :) = [RandVelX, RandVelY];
            scatterTime(n,1) = 0;
            
        elseif charge > 0 && positions(n, 1)<= center && ProbMoveWhenHappy > prob
            RandAngle = rand(1).*2.*pi;
            RandVelX = rand(1).*v_th.* cos(RandAngle);
            RandVelY = rand(1).*v_th.* sin(RandAngle);
            velocity(n, :) = [RandVelX, RandVelY];
            scatterTime(n,1) = 0;
            
            TempPos(n, 1) = positions(n, 1) + velocity(n, 1)*t;
            EnergyParticle = k*T+0.5*mass*(velocity(1)^2 +velocity(2)^2);
            
            
            if TempPos(n, 1) > center && EnergyParticle < 2%required energy of partciel is achieved
                velocity(n, :) = [0, 0];
            end
        end
    end
    
end

