function [ positions ] = updatePosition( numP, positions, velocity, t, xmin, xmax, ymin, ymax)
%updatePosition Summary of this function goes here
%   Detailed explanation goes here

%Boundary conditions
for j=1:2
    for k=1:numP
        revx(k) = 0;
        %move particles a each time step based on their velocity
        positions(k, j) = positions(k, j) + velocity(k, j)*t;
        
        %restrications of x-cordinate of each particle BC
        if j == 1 %(for all x cordinates)
            if positions(k, 1) <= xmin || positions(k, 1)>= xmax
                
                velocity(k, 1) = -1*velocity(k, 1);%just negate to reflect
                
                revx(k) = 1;
            end
        end
        
        %y parmaters of region BC reflective
        if j == 2
            if (positions(k, 2) <= ymin || positions(k, 2) >= ymax) %---&& %engery less than required to penetrate
                velocity(k, 2) = -1*velocity(k, 2);%just negate y component fr reflection
                
            end
        end
        
    end
end
end

