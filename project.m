%% Project Sonya Stuhec- Leonard 100963181

%define electron parameters
mass = 9.109E-31; %in kg
charge = -1602E-19; %in C
TauMN = 0.2E-12;

k = physconst('Boltzmann'); %Use of constants in matlab
T = 300; % temperature in Kalvin

%define thermal velocity (source:
%https://en.wikipedia.org/wiki/Thermal_velocity)
v_th = sqrt(k*T/mass);

numP = 100; %number of particles

%box definitions
xmax = 75;
xmin = 0;
ymax = 50;
ymin = 0;
center = xmax/2;


%use 100 steps to get across the region xmax long
t = (200e-9/v_th)/100;


%particle initalization
xmaxNano = 75e-9;
ymaxNano = 50e-9;

% %% Creating G-matrix
% %Create G matrix
% G = sparse(xmax*ymax);
% %B matrix defines boundary conditions
% B = zeros(1, xmax*ymax);
% %Applied oltage
% V0 = 0;
% 
% %populate G-matrix
% for i = 1:xmax
%     for j = 1:ymax
%         %numbering scheme for G
%         n = j+(i-1)*ymax;
%         
%         if i==1
%             %at x=0 Bc = V0, and djagonal of BC = 1
%             G(n, n) = V0;
%             B(n) = V0;
%             
%         elseif i==xmax
%             %at x=xmax BC=0, and djagonal of BC = 1
%             G(n, n) = V0;
%             B(n) = V0;
%             
%         elseif j==1
%             
%             G(n, n) = 1;
%             B(n) = 0;
%             
%         elseif j==ymax
%             G(n, n) = 1;
%             B(n) = 0;
%             
%         else
%             %solve FD equations or put 1 an d4 in the row coresponding to
%             %the j, i posjtjon
%             xmax_Minus1 = j+(i-2)*ymax;
%             xmax_Plus1 = j+ (i)*ymax;
%             ymax_Minus1 = j-1+(i-1)*ymax;
%             ymax_Plus1 = j+1+(i-1)*ymax;
%             
%             G(n, n)= -4;
%             G(n, xmax_Minus1) = 1;
%             G(n, xmax_Plus1) = 1;
%             G(n, ymax_Minus1) = 1;
%             G(n, ymax_Plus1) = 1;
%             
%         end
%     end
% end
% 
% %useful for quick verification of G
% figure (1)
% spy(G)
% title('Spying on G matrix')
% 
% V = (G)\B'; %inverse(G)*B
% 

%% Creating conductivity map and applying to G matrix
%populate conductivity map

%conductiity equation
condEq = charge^2*TauMN*numP/mass;
background = 1;

%initalize regions
[conductivtyN, posN, velN] = regionInit( 'electron', numP);
[conductivtyP, posP, velP] = regionInit( 'hole', numP);

% % positions = vertcat(posN, posP);
% % velocity = vertcat(velN, velP);
% 
% for i=1:xmax
%     for j=1:ymax
%         Map(i, j) = background;
%         %set conductivity of regions (NP from left to right)
%         if i <= center
%             Map(i, j) = conductivtyN;
%         else %(i>=45 &&<=75)
%             Map(i, j) = conductivtyP;
%         end
%     end
% end
% 
% %plot conductivity
% figure(2)
% surf(Map)
% title('Conductivity map')
% colorbar
% 
% % put conductivity map into n space
% for i=1:xmax
%     for j=1:ymax
%         n=j+(i-1)*ymax;
%         map_n(n) = Map(i, j);
%     end
% end
% 
% %conductivity acounted for G matrjx
% G2 = G.*map_n;
% 
% 
% 
% % Define Eigen vector of nx by ny plane
% V_vec = zeros(xmax, ymax);
% 
% for i=1:xmax
%     for j=1:ymax
%         n=j+(i-1)*ymax;
%         V_vec(i, j) = V(n);
%     end
% end
% 
% figure (3)
% surf(V_vec)
% title('Plot of E-feild G-matrix method')

%% Main loops for producing the "movie" of particles


%Probability of scattering
ProbScat = 1- exp(-t/TauMN);
scatterTime = zeros(numP, 1);
ProbMoveWhenHappy = 0.1;
iterations = 20;

v = VideoWriter('try1.avi');
figure (4)
axis([xmin, xmaxNano, ymin, ymaxNano])
plot(posN(:, 1), posN(:, 2), '.b')
hold on
plot(posP(:, 1), posN(:, 2), '.r')
hold off
pause(0.2)
title ('Simulation of Electron Trajectories')

for iter =1:iterations
    scatterTime= scatterTime+t*iter;
    
    %Probability of electrons scattering - should it move and how much?
    
    velN = ShouldItMove( t, posN, charge, center, ProbScat, ProbMoveWhenHappy, velN, mass);
    velP = ShouldItMove( t, posP, charge, center, ProbScat, ProbMoveWhenHappy, velP, mass);
    
    %Boundary conditions
    
    posN = updatePosition( numP, posN, velN, t, xmin, xmaxNano, ymin, ymaxNano);
    posP = updatePosition( numP, posP, velP, t, xmin, xmaxNano, ymin, ymaxNano);
    

    
    figure (4)
    axis([xmin, xmaxNano, ymin, ymaxNano])
    plot(posN(:, 1), posN(:, 2), '.b')
    hold on
    plot(posP(:, 1), posN(:, 2), '.r')
    hold off
    pause(0.2)
    title ('Simulation of Electron Trajectories')
    
    %Plot updated conductivity fo rmoving electrons in each resion
    
    CondMapUpdated = ConductivityCal(posN, 'electron', numP);
    
    
    % %         open(v);
    %
    %         figure(5)
    %         subplot(1, 2, 1)
    %         surf(CondMapUpdated)
    %         title ('Conductivity map changes with time')
    %         colorbar
    %         subplot(1, 2, 2)
    %         plot(posN(1:length(posN), 1), posN(1:length(posN), 2), '.b')
    %         hold on
    %         plot(posP(1:length(posP), 1), posN(1:length(posP), 2), '.r')
    %         axis([xmin, xmaxNano, ymin, ymaxNano])
    %         pause(0.2)
    %         title ('Simulation of Electron Trajectories')
    %
    % %         frame = getframe(gcf);
    % %         writeVideo(v,frame);
    % %         close(v)
    
    close(v);
    
    %temperature map
    %convert velocieites into temperatures, then use hist3 to bin and plot
    Temperature = (velN.^2).*(mass/k);
    
    %     figure(4)
    %     tempDist = hist3(Temperature, [50, 50]);
    %     pcolor(tempDist')
    %     title ('Electron temerature map')
    %     colorbar
    
    
    %calculate the mean free path of the electrons. The time between
    %collisions  is incimetned each iteratin at the top of the iter loop.
    MFP = mean(scatterTime(:, 1));
    
end
