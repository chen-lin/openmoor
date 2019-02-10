% Example to use libMoorApi in Matlab on MacOS
freq = 0.2;
time = 0:0.02:10/freq;
time_step = time(2:end) - time(1:end-1);
n_time = length(time);
displacement = zeros(n_time,6);
velocity = displacement;
mooring_load = zeros(n_time,6);
mooring_load_ptr = libpointer('doublePtr', mooring_load(1,:));
direction_id = 1;
amplitude = 2*pi*freq*1;

for j=1:n_time
    % Let the cable start from zero displacements
    if (time(j)<1/freq/2) 
        velocity(j, direction_id)=-0.5*amplitude*sin(2.0*pi*freq*time(j));
    else
        velocity(j, direction_id)= -amplitude*sin(2.0*pi*freq*time(j));
    end
    if j >= 2
        displacement(j, direction_id) = displacement(j-1, direction_id)+...
            (time(j)-time(j-1))*velocity(j, direction_id);
    end
end

plot(time, displacement)
%% Load library
loadlibrary('libMoorApi','moorapi.h');

%% Initialization.
input_file = 'CaseOC3.xml';
calllib('libMoorApi','initialize',input_file)
calllib('libMoorApi','get_cable_fairlead_force', 2); 

%% Time stepping.
tic
for j = 2:n_time
    % Cable state is saved every 0.1 s.
    calllib('libMoorApi','update', mooring_load_ptr, displacement(j,:)', ...
        velocity(j,:)', time(j), time_step(j-1), 0.1); 
    mooring_load(j,:) = mooring_load_ptr.value;
    display(time(j));
    display(mooring_load(j,:))
    % Check fairlead force if desired.
    % calllib('libMoorApi','get_cable_fairlead_force',0); 
end
toc
%% Finish simulation and clear variables.
% Clear memory
calllib('libMoorApi','finish');
unloadlibrary('libMoorApi');