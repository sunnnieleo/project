function [ position, velocity ] = randInit( numP, xmin, xmax, ymax, v_th )
%RANDINIT Randome initalization of particles in box
%   Proveided a number of particles numP, and box maxiimum dimensions
%   outputs the required number of particles inside the box centered at
%   origin
%   INPUTS: numP = number of particles 
%           xmax = max x value (starting at origin)
%           ymax = max y value (starting at origin)
%           v_th = thermal velocity of particles
%   OUTPUTS: position = matrix of x and y cordiates for numP particle
%            positions; x in collumn 1, y in collumn 2
%            velocity = matrix of x and y cordiates for numP particle
%            velocities; x in collumn 1, y in collumn 2


%initalize randome positions for particles
%randmone number between 0 and xmax or ymax to be within box
xlocations = rand(numP, 1).*(xmax-xmin)+xmin;
ylocations = rand(numP, 1).*ymax;
position = [xlocations, ylocations];

% randome velocity angle with magnitude v_th
% generate a randome inital angle for each particle in radians
angle = rand(numP, 1).*2.*pi;
velocityX = v_th.* cos(angle);
velocityY = v_th.* sin(angle);
velocity = [velocityX, velocityY];


end

