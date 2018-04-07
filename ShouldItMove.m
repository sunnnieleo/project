function [ velocity ] = ShouldItMove( t, positions, charge, center, ProbScat, ProbMoveWhenHappy, velocity, mass, positionsOther, chargeOther)
%ShouldItMove Summary of this function goes here
%   Detailed explanation goes here

k = physconst('Boltzmann'); %Use of constants in matlab
T = 300; % temperature in Kalvin
v_th = sqrt(k*T/mass);
minSep = 1e-9; %0.1 nm seperation for Foreces to effect particle
kcon = 9e9;

for n =1:length(positions)-1
    rx = abs(positions(n,1)-positions(:,1));
    ry = abs(positions(n,2)-positions(:,2));
    rxo = abs(positions(n,1)-positionsOther(:,1));
    ryo = abs(positions(n,2)-positionsOther(:,2));
end

for n=1:length(positions)-1
    if rx(n) <= minSep && ry(n) <=minSep
        
        accel = (kcon.*charge.^2)/(mass.*sqrt(rx(n)^2+ry(n)^2));
        velocity(n, :) = velocity(n, :) + accel.*t;
        positions(n, :) = positions(n, :) + velocity(n, :)*t;
        
    elseif rxo(n) <= minSep && ryo(n) <=minSep
        accel = (kcon.*charge.*chargeOther)./(mass.*sqrt(rx(n)^2+ry(n)^2));
        velocity(n, :) = velocity(n, :) + accel*t;
        positions(n, :) = positions(n, :) + velocity(n, :)*t;
        
    else
        
        prob = rand();
        
        if charge < 0 && positions(n, 1)< center && ProbScat < prob
            RandAngle = rand(1).*2.*pi;
            RandVelX = rand(1).*v_th.* cos(RandAngle);
            RandVelY = rand(1).*v_th.* sin(RandAngle);
            velocity(n, :) = [RandVelX, RandVelY];
            scatterTime(n,1) = 0;
            
        elseif charge<0 && positions(n, 1)>= center && ProbScat < prob && ProbMoveWhenHappy > prob
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
            
        elseif charge > 0 && positions(n, 1)<= center && ProbScat < prob && ProbMoveWhenHappy > prob
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

end

