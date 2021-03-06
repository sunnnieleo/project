function [ positions ] = updatePosition( vth, numP, positions, velocity, t, xmin, xmax, ymin, ymax, charge)
%updatePosition Summary of this function goes here
%   Detailed explanation goes here

%Boundary conditions
for j=1:2
    for k=1:numP
        Prevposition = positions(k, :);
        
        %move particles by randome velocit generated in ShouldItMove
        positions(k, j) = positions(k, j) + velocity(k, j)*t;
        
        %             %restrications of x-cordinate of each particle BC
        %             if j == 1 %(for all x cordinates)
        %                 if positions(k, 1) <= xmin
        %
        %                     positions(k, 1) = Prevposition(1, 1) + vth*t;
        %                 end
        %                 if positions(k, 1)>= xmax
        %                     positions(k, 1) = Prevposition(1, 1) - vth*t;
        %                 end
        %             end
        %
        %             if j == 2
        %                 if positions(k, 2) <= ymin
        %                     positions(k, 2) = Prevposition(1, 2)  + vth*t;;
        %                 end
        %                 if  positions(k, 2) >= ymax %---&& %engery less than required to penetrate
        %                     positions(k, 2) = Prevposition(1, 2)  - vth*t;;
        %                 end
        %             end
        
        
        if positions(k, 1) <= xmin || positions(k, 1)>= xmax || positions(k, 2) <= ymin || positions(k, 2) >= ymax
            %generate new position within bounds
            if charge > 0
            [positions(k, :), velocity(k, :) ] = randInit( 1, xmax/2, xmax, ymax, vth );
            else
                [positions(k, :), velocity(k, :) ] = randInit( 1, xmin, xmax/2, ymax, vth );
            end
        end
    end
end

