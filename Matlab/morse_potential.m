clear all
close all
clc

%%This code simulates an initial distribution of nodes 
%%moving based on morse potential

%%initial cytoplasm node placement
cyt_x = rand(1,7);
cyt_y = rand(1,7);
cN =length(cyt_x);
U_II = .488;
W_II = .146;
xsi_II = .3125;
gamma_II = 1.25;

%%time vector
dt = .01;
num_steps = 10000;
t = [0:dt:num_steps];
%halfway through time vector
halfway = ceil(length(t)/2);
time_points = ceil(length(t)/1000);
plot_points = 1000;
plotter = 0;

%holds morse potential energy between each pair of nodes
Energy_morse_pairs_cyt_to_cyt = zeros(16*cN*cN, time_points);
%holds morse potential energy at each cytoplasm node
Energy_morse_individual_cyt_nodes = zeros(4*cN, time_points);
%holds morse potential energy of whole system
Energy_morse_sum = zeros(1,time_points);
%holds energy of the whole system
Energy_total = zeros(1,time_points);


%aray(1:cN) are the vectors pointing from all cyt nodes to cyt node 1 etc
cyt_to_cyt_vecs_x = zeros(16*cN*cN, time_points);
cyt_to_cyt_vecs_y = zeros(16*cN*cN, time_points);
cyt_to_cyt_lengths = zeros(16*cN*cN, time_points);

%hold morse force between each pair of nodes
morse_pairs_x = zeros(16*cN*cN, time_points);
morse_pairs_y = zeros(16*cN*cN, time_points);
%holds total morse force on each cyt node at each time step
F_morse_cyt_nodes_x = zeros(4*cN*cN, time_points);
F_morse_cyt_nodes_y = zeros(4*cN*cN, time_points);
 
%holds sum of linear,rotational and morse forces
F_sum_x = zeros(4*cN,time_points);
F_sum_y = zeros(4*cN, time_points);
 
%array to hold changing node positions of cytoplasm
cyt_nodes_x = zeros(4*cN, time_points);
cyt_nodes_y = zeros(4*cN, time_points);
%initial cytoplasm node placement
cyt_nodes_x(1:cN,1) =cyt_x';
cyt_nodes_y(1:cN,1) =cyt_y'; 

%run it
for m = 1:length(t)-1
    
    %get the vector at the current time step between each cytoplasm node
    for i = 1:cN
        cyt_to_cyt_vecs_x(cN*(i-1)+1:cN*i,1) = cyt_nodes_x(i,1) -cyt_nodes_x(1:cN,1);
        cyt_to_cyt_vecs_y(cN*(i-1)+1:cN*i,1) = cyt_nodes_y(i,1) -cyt_nodes_y(1:cN,1);
    end
    
    %get the length of each vector
    cyt_to_cyt_lengths(:,1) = sqrt(cyt_to_cyt_vecs_x(:,1).^2 + cyt_to_cyt_vecs_y(:,1).^2);

    %Morse potential between nodes i and j is V = U*exp(-r/xsi) - W*exp(-r/gamma)
    
    %Energy between each pair of nodes
    Energy_morse_pairs_cyt_to_cyt(1:cN,1) = U_II*exp(-cyt_to_cyt_lengths(1:cN,1)/xsi_II) - W_II*exp(-cyt_to_cyt_lengths(1:cN,1)/gamma_II);
    
    %Energy on individual node is the sum
    for i = 1:cN
        Energy_morse_individual_cyt_nodes(i,1) = sum(Energy_morse_pairs_cyt_to_cyt(cN*(i-1)+1:i*cN,1),1);
    end
    
    %Energy of the whole cytoplasm system
    Energy_morse_sum(1) = sum(Energy_morse_individual_cyt_nodes(1:cN,1),1);
    %energy of the whole system
    Energy_total(1) = Energy_morse_sum(1);

    
    %morse force on each cytoplasm node 
    for i = 1:cN
        for j = 1:cN
            if (i == j)
               morse_pairs_x(cN*(i-1) +j, 1) = 0;
               morse_pairs_y(cN*(i-1) +j, 1) = 0;
            else
                morse_pairs_x(cN*(i-1) +j, 1) = get_morse_force(U_II, W_II,xsi_II, gamma_II, cyt_to_cyt_lengths(cN*(i-1) +j,1), -cyt_to_cyt_vecs_x(cN*(i-1) +j,1));
                morse_pairs_y(cN*(i-1) +j, 1) = get_morse_force(U_II, W_II,xsi_II, gamma_II, cyt_to_cyt_lengths(cN*(i-1) +j,1), -cyt_to_cyt_vecs_y(cN*(i-1) +j,1));
            end
        end
    
        
        F_morse_cyt_nodes_x(i,1) = sum(morse_pairs_x(cN*(i-1)+1:cN*i, 1),1);
        F_morse_cyt_nodes_y(i,1) = sum(morse_pairs_y(cN*(i-1)+1:cN*i, 1),1);
    end
    
    %the force doesnt change here because only modeling cytoplasm
    F_sum_x(1:cN,1) = F_morse_cyt_nodes_x(1:cN,1);
    F_sum_y(1:cN,1) = F_morse_cyt_nodes_y(1:cN,1);
    
    %update position
    cyt_nodes_x(1:cN,1)= cyt_nodes_x(1:cN,1) + F_sum_x(1:cN,1)*dt;
    cyt_nodes_y(1:cN,1)= cyt_nodes_y(1:cN,1) + F_sum_y(1:cN,1)*dt;
    
    %keep for plots
    if (mod(m,plot_points) ==1)
        plotter = plotter +1;
        cyt_nodes_x_plot(1:cN,plotter) = cyt_nodes_x(1:cN,1);
        cyt_nodes_y_plot(1:cN,plotter) = cyt_nodes_y(1:cN,1);
    end

end
