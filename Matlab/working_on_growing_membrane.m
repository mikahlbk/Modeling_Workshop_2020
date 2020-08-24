%% This code uses SCE model to simulate one plant cell
%%It models the cell wall nodes interacting with cytoplasm nodes
%%Author: Mikahl Banwarth-Kuhn
clear all
clc
close all

%initial wall node placement
%number of nodes on membrane at t=0
InitmembrNodeCount = 100;
%height of cell
a_h = -1.25;
b_h = 1.25;
%width of cell
a_w = -1.25;
b_w = 1.25;
%to equally distribute nodes along side walls
dh = (b_h-a_h)/(InitmembrNodeCount/4);
%to equally distribute nodes along end walls
dw = (b_w-a_w)/(InitmembrNodeCount/4);
x(1:InitmembrNodeCount/4) = a_w;
x = [x a_w:dw:b_w-dw];
wall_x = [x -x];
y = [a_h:dh:b_h-dh];
y(length(y)+1:length(y) +InitmembrNodeCount/4) = b_h;
wall_y = [y -y];
%keeps track of total number of wall nodes
wN = length(wall_x);
%keeps track of nodes in the four corners
corner_vec = [1 26 51 76];

%initial cytoplasm node placement
InitCellRadius = .9375;
InitCellNodeCount = 20;
random_angles = 2*pi*rand(InitCellNodeCount,1);
random_radius = InitCellRadius*rand(1,InitCellNodeCount);
cyt_x = random_radius.*cos(random_angles)';
cyt_y = random_radius.*sin(random_angles)';
cN =length(cyt_x);


%spring parameters
%neutral angles for rotational springs
%angles for corner versus flank nodes
theta_90 = pi/2;
theta_180 = pi;
%spring constant for rotational spring
k_theta = 6;
k_theta_flank = 6;
%spring constant for linear springs
%different on flanks vs. ends
k_spring_flank = 6;
k_spring_ends = 30;
%linear spring neutral length
spring_l_not = .0625;

%parameters for morse potential functions
%internal-internal nodes
U_II = .488;
W_II = .146484;
xsi_II = .3125;
gamma_II = 1.25;
%internal-membrane nodes
U_MI = 0.78125;
W_MI = 0.00;
xsi_MI = .125;
gamma_MI = .625;

%time vector
dt = .0003;
num_steps = 10;
t = [0:dt:num_steps];

%halfway through time vector
halfway = ceil(length(t)/2);
time_points = ceil(length(t)/100);
num_plots = 100;
mult = 4;
plotter = 0;
growth = 0;

%array to hold changing node positions of wall
wall_nodes_x = zeros(mult*wN, 1);
wall_nodes_y = zeros(mult*wN, 1);
%initial wall node placement
wall_nodes_x(1:wN,1) = wall_x';
wall_nodes_y(1:wN,1) = wall_y';
wall_nodes_x_plot = zeros(mult*wN,time_points);
wall_nodes_y_plot = zeros(mult*wN,time_points);

%array to hold changing node positions of cytoplasm
cyt_nodes_x = zeros(mult*cN, 1);
cyt_nodes_y = zeros(mult*cN, 1);
%initial cytoplasm node placement
cyt_nodes_x(1:cN,1) =cyt_x';
cyt_nodes_y(1:cN,1) =cyt_y'; 
cyt_nodes_x_plot = zeros(mult*cN, time_points);
cyt_nodes_y_plot = zeros(mult*cN, time_points);

%arrays to hold changing angles at each time step
thetas = zeros(mult*wN,1);
%holds angles to the left(right) of each node
thetas_left = zeros(mult*wN,1);
thetas_right = zeros(mult*wN,1);
%holds cos of angle at each node
cos_thetas = zeros(6*wN,1);
%holds cos of angle to the left(right) of each node
cos_thetas_left = zeros(mult*wN,1);
cos_thetas_right = zeros(mult*wN,1);
%holds cross product of each angle
crossZ = zeros(1,mult*wN);

%arrays to hold changing energy 
%holds the bending energy at each node
Energy_angle_individual = zeros(mult*wN, 1);
%holds the bending energy of the whole system
Energy_angle_sum = zeros(1,1);
%holds linear spring energy at each node
Energy_spring_individual = zeros(mult*wN, 1);
%holds linear spring energy of whole system
Energy_spring_sum = zeros(1,1);
%holds morse potential energy between each set of II nodes
Energy_morse_pairs_cyt_to_cyt = zeros(mult^2*cN*cN, 1);
%holds morse potential energy between each set of MI nodes
Energy_morse_pairs_cyt_to_wall = zeros(mult^2*cN*wN, 1);
%holds morse potential energy between each set of MI nodes
Energy_morse_pairs_wall_to_cyt = zeros(mult^2*wN*cN, 1);
%holds morse potential energy at each wall node
Energy_morse_individual_wall_nodes = zeros(mult*wN, 1);
%holds morse potential energy at each cytoplasm node
Energy_morse_individual_cyt_nodes = zeros(mult*cN, 1);
%holds morse potential energy of whole system
Energy_morse_sum = zeros(1,1);
%holds energy of the whole system
Energy_total = zeros(1,1);

%arrays to hold vectors needed for energy and force functions
%arrays used for linear and rotational spring functions
%holds x component of vector to left of each node
left_diff_x = zeros(mult*wN,1);
left_left_diff_x = zeros(mult*wN,1);
%holds x component of vector to right of each node
right_diff_x = zeros(mult*wN,1);
right_right_diff_x = zeros(mult*wN,1);
%holds y component of vector to left of each node
left_diff_y = zeros(mult*wN,1);
left_left_diff_y = zeros(mult*wN,1);
%holds y component of vector to right of each node
right_diff_y = zeros(mult*wN,1);
right_right_diff_y = zeros(mult*wN,1);
%holds length of vector to left of each node
left_lengths = zeros(mult*wN, 1);
left_left_lengths = zeros(mult*wN,1);
%holds length of vector to right of each node
right_lengths = zeros(mult*wN,1);
right_right_lengths = zeros(mult*wN,1);

%arrays used for morse functions
%aray(1:cN) are the vectors pointing from all cyt nodes to wall node 1 etc
wall_to_cyt_vecs_x = zeros(mult^2*wN*cN, 1);
wall_to_cyt_vecs_y = zeros(mult^2*wN*cN, 1);
wall_to_cyt_lengths = zeros(mult^2*wN*cN, 1);

%aray(1:wN) are the vectors pointing from all wall nodes to cyt node 1 etc
cyt_to_wall_vecs_x = zeros(mult^2*wN*cN, 1);
cyt_to_wall_vecs_y = zeros(mult^2*wN*cN, 1);
cyt_to_wall_lengths = zeros(mult^2*wN*cN, 1);

%aray(1:cN) are the vectors pointing from all cyt nodes to cyt node 1 etc
cyt_to_cyt_vecs_x = zeros(mult^2*cN*cN, 1);
cyt_to_cyt_vecs_y = zeros(mult^2*cN*cN, 1);
cyt_to_cyt_lengths = zeros(mult^2*cN*cN, 1);


%arrays to hold linear spring forces
%holds force in x direction of spring to right of each node
F_spring_x_right = zeros(mult*wN, 1);
%holds force in y direction of spring to right of each node
F_spring_y_right = zeros(mult*wN, 1);
%holds force in x direction of spring to left of each node
F_spring_x_left = zeros(mult*wN, 1);
%holds force in y direction of spring to left of each node
F_spring_y_left = zeros(mult*wN, 1);



%arays to hold rotational spring forces
%holds four terms derived by taking gradient of potential function
%x components of each term
term1x = zeros(mult*wN, 1);
term2x = zeros(mult*wN, 1);
term3x = zeros(mult*wN, 1);
term4x = zeros(mult*wN, 1);
term1lx = zeros(mult*wN,1);
term2lx = zeros(mult*wN,1);
term1rx = zeros(mult*wN,1);
term2rx = zeros(mult*wN,1);
%y components of each term
term1y = zeros(mult*wN,1);
term2y = zeros(mult*wN,1);
term3y = zeros(mult*wN,1);
term4y = zeros(mult*wN,1);
term1ly = zeros(mult*wN,1);
term2ly = zeros(mult*wN,1);
term1ry = zeros(mult*wN,1);
term2ry = zeros(mult*wN,1);
%constant out front of sum of terms for force on each node
constant = zeros(mult*wN,1);
constant_left = zeros(mult*wN,1);
constant_right = zeros(mult*wN,1);
%holds constant*sum of terms 1-4 x component
F_angle_x = zeros(mult*wN,1);
%holds constant*sum of terms 1-4 y component
F_angle_y = zeros(mult*wN,1);

%is there a better solution than this
morse_force_pairs_cyt_to_cyt_x = zeros(mult^2*cN*cN,1);
morse_force_pairs_cyt_to_cyt_y = zeros(mult^2*cN*cN,1);
morse_force_pairs_cyt_to_wall_x= zeros(mult^2*cN*wN,1);
morse_force_pairs_cyt_to_wall_y= zeros(mult^2*cN*wN,1);
morse_force_pairs_wall_to_cyt_x= zeros(mult^2*cN*wN,1);
morse_force_pairs_wall_to_cyt_y= zeros(mult^2*cN*wN,1);


%arrays to hold morse potential forces
F_morse_wall_nodes_x = zeros(mult*wN,1);
F_morse_wall_nodes_y = zeros(mult*wN,1);
F_morse_cyt_nodes_x = zeros(mult*cN,1);
F_morse_cyt_nodes_y = zeros(mult*cN,1);

%holds sum of linear,rotational and morse forces
F_sum_wall_x = zeros(mult*wN,1);
F_sum_wall_y = zeros(mult*wN,1);
F_sum_cyt_x = zeros(mult*cN,1);
F_sum_cyt_y = zeros(mult*cN,1);


%run it
for m = 1:length(t)-1
     if (m < 16666)
         growth = 1;
     else
         growth = 0;
     end
     if (growth)
    if(mod(m,200) == 0)
         wN = wN+1;
         [M,I] = max(right_lengths);
         if (I < corner_vec(2))
             corner_vec(2:4) = corner_vec(2:4) +1;
         elseif (I >= corner_vec(2))&&(I< corner_vec(3))
             corner_vec(3:4) = corner_vec(3:4) +1;
         elseif (I >= corner_vec(3)) && (I<corner_vec(4))
             corner_vec(4) = corner_vec(4) +1;
         end
         if (I == wN-1)
             new_wall_node_x = (wall_nodes_x(1,1)+wall_nodes_x(wN-1,1))/2;
             new_wall_node_y = (wall_nodes_y(1,1)+wall_nodes_y(wN-1,1))/2;
             wall_nodes_x(wN,1) = new_wall_node_x;
             wall_nodes_y(wN,1) = new_wall_node_y;
         else
            new_wall_node_x = (wall_nodes_x(I+1,1)+wall_nodes_x(I,1))/2;
            new_wall_node_y = (wall_nodes_y(I+1,1)+wall_nodes_y(I,1))/2;
            wall_nodes_x(I+2:wN) = wall_nodes_x(I+1:wN-1);
            wall_nodes_y(I+2:wN) = wall_nodes_y(I+1:wN-1);
            wall_nodes_x(I+1,1) = new_wall_node_x;
            wall_nodes_y(I+1,1) = new_wall_node_y;
         end
    end
  
    
    if (mod(m,1666) == 0)
        cN = cN+1;
        CurrentCellRadius = min(wall_nodes_x(corner_vec(4),1), wall_nodes_y(corner_vec(2),1));
        random_angles = 2*pi*rand(1,1);
        random_radius = CurrentCellRadius*rand(1,1);
        cyt_nodes_x(cN,1) = random_radius.*cos(random_angles);
        cyt_nodes_y(cN,1) = random_radius.*sin(random_angles);
    end
    end
%%%%compute updated vectors
    %vectors for linear and rotational spring functions
    %these are the x values of the vectors to the left of each wall node
    left_diff_x(2:wN,1) = wall_nodes_x(1:wN-1,1) - wall_nodes_x(2:wN,1);
    left_diff_x(1,1) = wall_nodes_x(wN,1) - wall_nodes_x(1,1);
    %these are the x values of the vectors to the right of each wall node
    right_diff_x(1:wN-1,1) = wall_nodes_x(2:wN,1) - wall_nodes_x(1:wN-1,1);
    right_diff_x(wN,1) = wall_nodes_x(1,1) - wall_nodes_x(wN,1);
    %these are the y values of the vectors to the left of each wall node
    left_diff_y(2:wN,1) = wall_nodes_y(1:wN-1,1) - wall_nodes_y(2:wN,1);
    left_diff_y(1,1) =  wall_nodes_y(wN,1) - wall_nodes_y(1,1);
    %these are the y values of the vectors to the right of each wall node
    right_diff_y(1:wN-1,1) = wall_nodes_y(2:wN,1) - wall_nodes_y(1:wN-1,1);
    right_diff_y(wN,1) = wall_nodes_y(1,1) - wall_nodes_y(wN,1);
    
%%%%get lengths of updated vectors
    left_lengths(1:wN,1) = sqrt(left_diff_x(1:wN,1).^2 + left_diff_y(1:wN,1).^2);
    right_lengths(1:wN,1) = sqrt(right_diff_x(1:wN,1).^2 + right_diff_y(1:wN,1).^2);
    

    %reshape to get left_left and right_right
    %these are the x values of the vectors two places to the left of each wall node
    left_left_diff_x(2:wN,1) = left_diff_x(1:wN-1,1);
    left_left_diff_x(1,1) = left_diff_x(wN,1);
    %these are the y values of the vectors two places to the left of each wall node
    left_left_diff_y(2:wN,1) = left_diff_y(1:wN-1,1);
    left_left_diff_y(1,1) = left_diff_y(wN,1);
    %these are the x values of the vectors two places to the right of each wall node
    right_right_diff_x(1:wN-1,1) = right_diff_x(2:wN,1);
    right_right_diff_x(wN,1) = right_diff_x(1,1);
    %these are the y values of the vectors two places to the right of each wall node
    right_right_diff_y(1:wN-1,1) = right_diff_y(2:wN,1);
    right_right_diff_y(wN,1) = right_diff_y(1,1);
    %reshape to get left_left_lengths and right_right_lengths
    left_left_lengths(2:wN,1) = left_lengths(1:wN-1,1);
    left_left_lengths(1,1) = left_lengths(wN,1);
    right_right_lengths(1:wN-1,1) = right_lengths(2:wN,1);
    right_right_lengths(wN,1) = right_lengths(1,1);
    
   
    %vectors for morse potential functions
    for i = 1:cN
        cyt_to_cyt_vecs_x(cN*(i-1)+1:cN*i,1) = cyt_nodes_x(i,1) - cyt_nodes_x(1:cN,1);
        cyt_to_cyt_vecs_y(cN*(i-1)+1:cN*i,1) = cyt_nodes_y(i,1) - cyt_nodes_y(1:cN,1);
        cyt_to_wall_vecs_x(wN*(i-1)+1:wN*i,1) = cyt_nodes_x(i,1) - wall_nodes_x(1:wN,1);
        cyt_to_wall_vecs_y(wN*(i-1)+1:wN*i,1) = cyt_nodes_y(i,1) - wall_nodes_y(1:wN,1);
        cyt_to_cyt_lengths(cN*(i-1)+1:cN*i,1) = sqrt(cyt_to_cyt_vecs_x(cN*(i-1)+1:cN*i,1).^2 + cyt_to_cyt_vecs_y(cN*(i-1)+1:cN*i,1).^2);
        cyt_to_wall_lengths(wN*(i-1)+1:wN*i,1) = sqrt(cyt_to_wall_vecs_x(wN*(i-1)+1:wN*i,1).^2 + cyt_to_wall_vecs_y(wN*(i-1)+1:wN*i,1).^2);
    end
    for i = 1:wN
        wall_to_cyt_vecs_x(cN*(i-1)+1:cN*i,1) = wall_nodes_x(i,1) - cyt_nodes_x(1:cN,1);
        wall_to_cyt_vecs_y(cN*(i-1)+1:cN*i,1) = wall_nodes_y(i,1) - cyt_nodes_y(1:cN,1);
        wall_to_cyt_lengths(cN*(i-1)+1:cN*i,1) = sqrt(wall_to_cyt_vecs_x(cN*(i-1)+1:cN*i,1).^2 + wall_to_cyt_vecs_y(cN*(i-1)+1:cN*i,1).^2);
    end
    
    
%%%get updated angle between nodes
    for i = 1:wN
        [crossZ(i),thetas(i,1), cos_thetas(i,1)] = get_angle([left_diff_x(i,1) left_diff_y(i,1)], [right_diff_x(i,1) right_diff_y(i,1)]);
    end
    
    %reshape to left and right angles
    thetas_left(2:wN,1) = thetas(1:wN-1,1);
    thetas_left(1,1) = thetas(wN,1);
    thetas_right(1:wN-1,1) = thetas(2:wN,1);
    thetas_right(wN,1) = thetas(1,1);
    cos_thetas_left(2:wN,1) = cos_thetas(1:wN-1,1);
    cos_thetas_left(1,1) = cos_thetas(wN,1);
    cos_thetas_right(1:wN-1,1) = cos_thetas(2:wN,1);
    cos_thetas_right(wN,1) = cos_thetas(1,1);


%%%%compute energy based on updated angle and vector lenghts
    
% %Angle potential function for node j is V = .5*k_ijk(theta_ijk-theta_eq)^2
%     %%Flank nodes have equilibrium angle of 180 deg
%     Energy_angle_individual(:,1) = .5*k_theta_flank.*(thetas(:,1)-theta_180).^2;
%     %%Corner nodes have equilibrium angle of 90 deg
%     Energy_angle_individual(corner_vec(:),1) = .5*k_theta.*(thetas(corner_vec(:),1)-theta_90).^2;
%     %Energy of all rotational springs
%     Energy_angle_sum(1) = sum(Energy_angle_individual(:,1),1);
   
%     %%Spring potential for r_ij is V = .5*k_ij(||r_ij|| -length_eq)^2
%     Energy_spring_individual(:,1) = .5.*k_spring_flank.*(right_lengths(:,1)-spring_l_not).^2;
%     Energy_spring_individual(corner_vec(2):corner_vec(3)-1,1) = .5.*k_spring_ends.*(right_lengths(corner_vec(2):corner_vec(3)-1,1)-spring_l_not).^2;
%     Energy_spring_individual(corner_vec(3):wN,1) = .5.*k_spring_ends.*(right_lengths(corner_vec(3):wN,1)-spring_l_not).^2;
%     Energy_spring_sum(1) = sum(Energy_spring_individual(:,1),1);
%     
%     %%morse potential between nodes i and j is V = U*exp(-r/xsi) - W*exp(-r/gamma)
%     %Energy between II nodes
%     Energy_morse_pairs_cyt_to_cyt(:,1) = U_II*exp(-cyt_to_cyt_lengths(:,1)/xsi_II) - W_II*exp(-cyt_to_cyt_lengths(:,1)/gamma_II);
%     %Energy between MI, cyt to wall nodes
%     Energy_morse_pairs_cyt_to_wall(:,1) = U_MI*exp(-cyt_to_wall_lengths(:,1)/xsi_MI) - W_MI*exp(-cyt_to_wall_lengths(:,1)/gamma_MI);
%     %Energy between MI, wall to cyt nodes 
%     Energy_morse_pairs_wall_to_cyt(:,1) = U_MI*exp(-wall_to_cyt_lengths(:,1)/xsi_MI) - W_MI*exp(-wall_to_cyt_lengths(:,1)/gamma_MI);
%     
%     %Energy on individual node is the sum
%     %for cytoplasm nodes
%     for i = 1:cN
%         Energy_morse_individual_cyt_nodes(i,1) = sum(Energy_morse_pairs_cyt_to_cyt(cN*(i-1)+1:i*cN,1),1) + sum(Energy_morse_pairs_cyt_to_wall(wN*(i-1)+1:i*wN,1),1); 
%     end
%     
%     %for wall nodes
%     for i = 1:wN
%         Energy_morse_individual_wall_nodes(i,1) = sum(Energy_morse_pairs_wall_to_cyt(cN*(i-1)+1:i*cN,1),1);
%     end
%     
%     %energy of the whole cytoplasm system
%     Energy_morse_sum(1) = sum(Energy_morse_individual_wall_nodes(:,1),1) + sum(Energy_morse_individual_cyt_nodes(:,1),1);
%     
%     %Energy of the whole system is sum of linear, cytoplasm, and rotational
%     Energy_total(1) = Energy_spring_sum(1) + Energy_angle_sum(1) +  Energy_morse_sum(1);

%%%%compute the forces
    %BENDING
    %constant out front in terms of theta and equilibrium angle
    checker = 1-cos_thetas(1:wN,1).^2;
    zeroes = find(~checker);
    constant(1:wN,1) = k_theta_flank*(thetas(1:wN,1)-theta_180)./sqrt(1-(cos_thetas(1:wN,1)).^2);
    constant(corner_vec(:),1) = k_theta*(thetas(corner_vec(:),1)-theta_90)./sqrt(1-(cos_thetas(corner_vec(:),1)).^2);
    

     for i = 1:length(zeroes)
            constant(zeroes(i),1) = 0;%k_theta_flank*(thetas(i,1)-theta_180)./(sqrt(1-(cos_thetas(i,1)).^2) + .0001);
     end
    
    %reshape to get constant_left and constant_right
    constant_left(2:wN,1) = constant(1:wN-1,1);
    constant_left(1,1) = constant(wN,1);
    constant_right(1:wN-1,1) = constant(2:wN,1);
    constant_right(wN,1) = constant(1,1);
        
    %force for each node is computed with four terms
    %do the x component first
    term1x(1:wN,1) = (-left_diff_x(1:wN,1))./(left_lengths(1:wN,1).*right_lengths(1:wN,1));
    term2x(1:wN,1) = cos_thetas(1:wN,1).*left_diff_x(1:wN,1)./(left_lengths(1:wN,1).^2);
    term3x(1:wN,1) = (-right_diff_x(1:wN,1))./(left_lengths(1:wN,1).*right_lengths(1:wN,1));
    term4x(1:wN,1) = cos_thetas(1:wN,1).*right_diff_x(1:wN,1)./(right_lengths(1:wN,1).^2);
    term1lx(1:wN,1) = (left_left_diff_x(1:wN,1))./(left_left_lengths(1:wN,1).*left_lengths(1:wN,1));
    term2lx(1:wN,1) = -cos_thetas_left(1:wN,1).*-left_diff_x(1:wN,1)./(left_lengths(1:wN,1).^2);
    term1rx(1:wN,1) = (right_right_diff_x(1:wN,1))./(right_right_lengths(1:wN,1).*right_lengths(1:wN,1));
    term2rx(1:wN,1) = -cos_thetas_right(1:wN,1).*-right_diff_x(1:wN,1)./(right_lengths(1:wN,1).^2);
    %do the y component second
    term1y(1:wN,1) = (-left_diff_y(1:wN,1))./(left_lengths(1:wN,1).*right_lengths(1:wN,1));
    term2y(1:wN,1) = cos_thetas(1:wN,1).*left_diff_y(1:wN,1)./(left_lengths(1:wN,1).^2);
    term3y(1:wN,1) = (-right_diff_y(1:wN,1))./(left_lengths(1:wN,1).*right_lengths(1:wN,1));
    term4y(1:wN,1) = cos_thetas(1:wN,1).*right_diff_y(1:wN,1)./(right_lengths(1:wN,1).^2);
    term1ly(1:wN,1) = (left_left_diff_y(1:wN,1))./(left_left_lengths(1:wN,1).*left_lengths(1:wN,1));
    term2ly(1:wN,1) = -cos_thetas_left(1:wN,1).*-left_diff_y(1:wN,1)./(left_lengths(1:wN,1).^2);
    term1ry(1:wN,1) = (right_right_diff_y(1:wN,1))./(right_right_lengths(1:wN,1).*right_lengths(1:wN,1));
    term2ry(1:wN,1) = -cos_thetas_right(1:wN,1).*-right_diff_y(1:wN,1)./(right_lengths(1:wN,1).^2);
    
    %sum the terms
    F_angle_x(1:wN,1) = constant(1:wN,1).*(term1x(1:wN,1) + term2x(1:wN,1) + term3x(1:wN,1) + term4x(1:wN,1))+ constant_left(1:wN,1).*(term1lx(1:wN,1) + term2lx(1:wN,1)) + constant_right(1:wN,1).*(term1rx(1:wN,1) + term2rx(1:wN,1));
    F_angle_y(1:wN,1) = constant(1:wN,1).*(term1y(1:wN,1) + term2y(1:wN,1) + term3y(1:wN,1) + term4y(1:wN,1))+ constant_left(1:wN,1).*(term1ly(1:wN,1) + term2ly(1:wN,1)) + constant_right(1:wN,1).*(term1ry(1:wN,1) + term2ry(1:wN,1));
    
   
    %to get the right forces for concave angles
    for i=1:wN
        if crossZ(i) == 1
            F_angle_x(i,1) = -F_angle_x(i,1);
            F_angle_y(i,1) = -F_angle_y(i,1);
        end
    end 
    
    %LINEAR FORCES
    %spring right
    F_spring_x_right(1:wN,1) = (k_spring_flank.*(right_lengths(1:wN,1)-spring_l_not)).*(right_diff_x(1:wN,1))./right_lengths(1:wN,1);
    F_spring_y_right(1:wN,1) = (k_spring_flank.*(right_lengths(1:wN,1)-spring_l_not)).*(right_diff_y(1:wN,1))./right_lengths(1:wN,1);
    %change spring constant on ends
    F_spring_x_right(corner_vec(2):corner_vec(3)-1,1) = (k_spring_ends.*(right_lengths(corner_vec(2):corner_vec(3)-1,1)-spring_l_not)).*(right_diff_x(corner_vec(2):corner_vec(3)-1,1))./right_lengths(corner_vec(2):corner_vec(3)-1,1);
    F_spring_y_right(corner_vec(2):corner_vec(3)-1,1) = (k_spring_ends.*(right_lengths(corner_vec(2):corner_vec(3)-1,1)-spring_l_not)).*(right_diff_y(corner_vec(2):corner_vec(3)-1,1))./right_lengths(corner_vec(2):corner_vec(3)-1,1);
    F_spring_x_right(corner_vec(4):wN,1) = (k_spring_ends.*(right_lengths(corner_vec(4):wN,1)-spring_l_not)).*(right_diff_x(corner_vec(4):wN,1))./right_lengths(corner_vec(4):wN,1);
    F_spring_y_right(corner_vec(4):wN,1) = (k_spring_ends.*(right_lengths(corner_vec(4):wN,1)-spring_l_not)).*(right_diff_y(corner_vec(4):wN,1))./right_lengths(corner_vec(4):wN,1);


    %spring left
    F_spring_x_left(1:wN,1) = (k_spring_flank.*(left_lengths(1:wN,1)-spring_l_not)).*(left_diff_x(1:wN,1))./left_lengths(1:wN,1);
    F_spring_y_left(1:wN,1) = (k_spring_flank.*(left_lengths(1:wN,1)-spring_l_not)).*(left_diff_y(1:wN,1))./left_lengths(1:wN,1);
    %change spring constant on ends
    F_spring_x_left(corner_vec(2)+1:corner_vec(3),1) = (k_spring_ends.*(left_lengths(corner_vec(2)+1:corner_vec(3),1)-spring_l_not)).*(left_diff_x(corner_vec(2)+1:corner_vec(3),1))./left_lengths(corner_vec(2)+1:corner_vec(3),1);
    F_spring_y_left(corner_vec(2)+1:corner_vec(3),1) = (k_spring_ends.*(left_lengths(corner_vec(2)+1:corner_vec(3),1)-spring_l_not)).*(left_diff_y(corner_vec(2)+1:corner_vec(3),1))./left_lengths(corner_vec(2)+1:corner_vec(3),1);
    F_spring_x_left(corner_vec(4)+1:wN,1) = (k_spring_ends.*(left_lengths(corner_vec(4)+1:wN,1)-spring_l_not)).*(left_diff_x(corner_vec(4)+1:wN,1))./left_lengths(corner_vec(4)+1:wN,1);
    F_spring_y_left(corner_vec(4)+1:wN,1) = (k_spring_ends.*(left_lengths(corner_vec(4)+1:wN,1)-spring_l_not)).*(left_diff_y(corner_vec(4)+1:wN,1))./left_lengths(corner_vec(4)+1:wN,1);
    F_spring_x_left(corner_vec(1),1) = (k_spring_ends.*(left_lengths(corner_vec(1),1)-spring_l_not)).*(left_diff_x(corner_vec(1),1))./left_lengths(corner_vec(1),1);
    F_spring_y_left(corner_vec(1),1) = (k_spring_ends.*(left_lengths(corner_vec(1),1)-spring_l_not)).*(left_diff_y(corner_vec(1),1))./left_lengths(corner_vec(1),1);

   
    %MORSE FORCES
    %Morse sum on wall node
    for i = 1:wN
        for j= 1:cN
            if (i == j)
               morse_force_pairs_wall_to_cyt_x(cN*(i-1) +j, 1) = 0;
               morse_force_pairs_wall_to_cyt_y(cN*(i-1) +j, 1) = 0;
            else
                morse_force_pairs_wall_to_cyt_x(cN*(i-1) +j, 1) = get_morse_force(U_MI, W_MI,xsi_MI, gamma_MI, wall_to_cyt_lengths(cN*(i-1) + j,1), -wall_to_cyt_vecs_x(cN*(i-1) +j,1));
                morse_force_pairs_wall_to_cyt_y(cN*(i-1) +j, 1) = get_morse_force(U_MI, W_MI,xsi_MI, gamma_MI, wall_to_cyt_lengths(cN*(i-1) + j,1), -wall_to_cyt_vecs_y(cN*(i-1) +j,1));
            end 
        end
        F_morse_wall_nodes_x(i,1) = sum(morse_force_pairs_wall_to_cyt_x(cN*(i-1)+1:cN*i, 1),1);
        F_morse_wall_nodes_y(i,1) = sum(morse_force_pairs_wall_to_cyt_y(cN*(i-1)+1:cN*i, 1),1);
    end
    for i = 1:cN
        for j= 1:cN
            if (i == j)
               morse_force_pairs_cyt_to_cyt_x(cN*(i-1) +j, 1) = 0;
               morse_force_pairs_cyt_to_cyt_y(cN*(i-1) +j, 1) = 0;
            else
            morse_force_pairs_cyt_to_cyt_x(cN*(i-1) +j, 1) = get_morse_force(U_II, W_II,xsi_II, gamma_II, cyt_to_cyt_lengths(cN*(i-1) + j,1), -cyt_to_cyt_vecs_x(cN*(i-1) +j,1));
            morse_force_pairs_cyt_to_cyt_y(cN*(i-1) +j, 1) = get_morse_force(U_II, W_II,xsi_II, gamma_II, cyt_to_cyt_lengths(cN*(i-1) + j,1), -cyt_to_cyt_vecs_y(cN*(i-1) +j,1));
        
            end
        end
        
        for j = 1:wN
            if (i == j)
               morse_force_pairs_cyt_to_wall_x(cN*(i-1) +j, 1) = 0;
               morse_force_pairs_cyt_to_wall_y(cN*(i-1) +j, 1) = 0;
            else
            morse_force_pairs_cyt_to_wall_x(wN*(i-1) +j, 1) = get_morse_force(U_MI, W_MI,xsi_MI, gamma_MI, cyt_to_wall_lengths(wN*(i-1) + j,1), -cyt_to_wall_vecs_x(wN*(i-1) +j,1));
            morse_force_pairs_cyt_to_wall_y(wN*(i-1) +j, 1) = get_morse_force(U_MI, W_MI,xsi_MI, gamma_MI, cyt_to_wall_lengths(wN*(i-1) + j,1), -cyt_to_wall_vecs_y(wN*(i-1) +j,1));
            end
        end
        F_morse_cyt_nodes_x(i,1) = sum(morse_force_pairs_cyt_to_cyt_x(cN*(i-1)+1:cN*i, 1),1)+sum(morse_force_pairs_cyt_to_wall_x(wN*(i-1)+1:wN*i, 1),1);
        F_morse_cyt_nodes_y(i,1) = sum(morse_force_pairs_cyt_to_cyt_y(cN*(i-1)+1:cN*i, 1),1)+sum(morse_force_pairs_cyt_to_wall_y(wN*(i-1)+1:wN*i, 1),1);
    end
        
   
        
    %sum
    F_sum_wall_x(1:wN,1)= F_angle_x(1:wN,1) + F_spring_x_right(1:wN,1) + F_spring_x_left(1:wN,1)+ F_morse_wall_nodes_x(1:wN,1);
    F_sum_wall_y(1:wN,1)= F_angle_y(1:wN,1) + F_spring_y_right(1:wN,1) + F_spring_y_left(1:wN,1)+ F_morse_wall_nodes_y(1:wN,1);
    F_sum_cyt_x(1:cN,1) = F_morse_cyt_nodes_x(1:cN,1);
    F_sum_cyt_y(1:cN,1) = F_morse_cyt_nodes_y(1:cN,1);

    %update position
    wall_nodes_x(1:wN,1) = wall_nodes_x(1:wN,1) + F_sum_wall_x(1:wN,1)*dt;
    wall_nodes_y(1:wN,1) = wall_nodes_y(1:wN,1) + F_sum_wall_y(1:wN,1)*dt;
    cyt_nodes_x(1:cN,1)= cyt_nodes_x(1:cN,1) + F_sum_cyt_x(1:cN,1)*dt;
    cyt_nodes_y(1:cN,1)= cyt_nodes_y(1:cN,1) + F_sum_cyt_y(1:cN,1)*dt;
    
   
    
      if (mod(m,num_plots)==1)
        plotter = plotter+1;
        wall_nodes_x_plot(1:wN,plotter) = wall_nodes_x(1:wN,1);
        wall_nodes_y_plot(1:wN,plotter) = wall_nodes_y(1:wN,1);
        cyt_nodes_x_plot(1:cN,plotter) = cyt_nodes_x(1:cN,1);
        cyt_nodes_y_plot(1:cN,plotter) = cyt_nodes_y(1:cN,1);
        
      end


end

%PLOTTING
%  wall_nodes_x_plot = [wall_nodes_x_plot(1:wN,1) ; wall_nodes_x_plot(1,1)];
%  wall_nodes_y_plot = [wall_nodes_y_plot(1:wN,1) ; wall_nodes_y_plot(1,1)];

%to plot with static axis size
figure
h = plot(NaN, NaN, '-go'); %// empty plot
hold on
g = plot(NaN, NaN, 'r*');
%axis([-2 2 -2 2])
%axis manual %// this line freezes the axes
for i = 1:1000
x = wall_nodes_x_plot(:,i);
y = wall_nodes_y_plot(:,i);
p = cyt_nodes_x_plot(1:cN, i);
m= cyt_nodes_y_plot(1:cN,i);
set(h, 'XData', x, 'YData', y)

set(g, 'XData', p, 'YData',m)
pause(.1);
end

%%to plot changing axis size
% figure
% for i=1:length(t)
%   plot(wall_nodes_x_plot(:,i), wall_nodes_y_plot(:,i), '-go')
%   hold on
%   plot(cyt_nodes_x(:,i), cyt_nodes_y(:,i), 'r*')
%   pause(1);
%   hold off
%      
% end

%morse Potential Plots between select cytoplasm nodes
% for i = 1:cN
% figure 
% plot(cyt_to_cyt_lengths(2*cN+i,:), Energy_morse_pairs_cyt_to_cyt(2*cN+i,:), 'go');
% end
% for i = 1:cN
% figure 
% plot(wall_to_cyt_lengths(2*cN+i,:), Energy_morse_pairs_wall_to_cyt(2*cN+i,:), 'r*');
% end

%plot energy change over time
% figure
%         plot([0:1:length(t)-1],Energy_spring_sum(:), 'b')
%         hold on
%         plot([0:1:length(t)-1],Energy_angle_sum(:),'k')
%         hold on 
%         plot([0:1:length(t)-1],Energy_morse_sum(:),'g')
%         hold on 
%         plot([0:1:length(t)-1],Energy_total(:),'r')
%         hold off

% % % % % % Make a simulation to save
% myVideo = VideoWriter('myfile.avi');
% myVideo.FrameRate = 50;
% open(myVideo);
% for i=1:plotter
% plot(wall_nodes_x_plot(1:wN,i), wall_nodes_y_plot(1:wN,i), '-go')
% hold on
% plot(cyt_nodes_x_plot(1:cN,i), cyt_nodes_y_plot(1:cN,i), 'r*')
% axis([-10 10 -10 10])
% F = getframe(gcf);
% writeVideo(myVideo,F);
% hold off
% end
% close(myVideo)