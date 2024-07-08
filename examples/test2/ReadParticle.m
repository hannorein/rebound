%% read binary file from rebound simulation
function [N_particle, time, N_colDidy, N_colDimor, N_escape, r_dust, particle_all, data_c, data_p] = ReadParticle(file_particle, file_collide)

% read particle file
fileID = fopen(file_particle,'r');

i = 1;
while ~feof(fileID)
    N_particle_new = fread(fileID,1,'int'); % Number of particles
    if isempty(N_particle_new)
        break;
    end
    time_new = fread(fileID,1,'double'); % time, s
    r_dust = fread(fileID,1,'double'); % dust particle size, m
    data_new = fread(fileID,[7,N_particle_new],'double'); 
    data_new = data_new';

    N_particle(i) = N_particle_new;
    time(i) = time_new;
    data_p{i} = data_new;
    i = i+1;
end

fclose(fileID);

% read collision file

fileID = fopen(file_collide,'r');

i = 1;
while ~feof(fileID)
    flag_remove = fread(fileID,1,'int');
    if isempty(flag_remove)
        break;
    end
    p_id = fread(fileID,1,'int');
    time_new = fread(fileID,1,'double');
    data_c(i,:) = [flag_remove, p_id, time_new]; 
    i = i+1;
end

% process data into individual particle information
particle = struct('ID', 0, 'id_collide', 0, 'time', 0, 'pos', [0,0,0], 'vel', [0,0,0]);
N_output = size(N_particle,2);
particle_all = repmat(particle,N_particle(1),1);
for n = 1:N_particle(1)  % initial condition
    particle_all(n).ID = data_p{1}(n,1);
    if abs( n - particle_all(n).ID ) > 0.1
        error('The hash of the particles are not in the expected order: 1, 2, 3, ...');
    end
    particle_all(n).id_collide = 0;
    particle_all(n).time = time(end);
    particle_all(n).pos = data_p{1}(n,2:4);
    particle_all(n).vel = data_p{1}(n,5:7);
end


for i = 2:N_output
    clc;
    fprintf( 'No. %d of output frame has been analyzed (%.1f%%)...\n', i, i/N_output*100.0 );
    for n = 1:N_particle(i)
        p_id = data_p{i}(n,1);
        length_0 =  size(particle_all(p_id).pos,1) + 1;
        particle_all(p_id).pos(length_0,:) = data_p{i}(n,2:4);
        particle_all(p_id).vel(length_0,:) = data_p{i}(n,5:7);
    end
end


i_colDidy = 1;
i_colDimor = 1;
i_escape = 1;
N_colDidy = [];
N_colDimor = [];
N_escape = [];
for i = 1:size(data_c,1)
    particle_all(data_c(i,2)).id_collide = data_c(i,1);
    particle_all(data_c(i,2)).time = data_c(i,3);
    switch data_c(i,1)
        case 1
            N_colDidy(i_colDidy, 1) = data_c(i,3);
            N_colDidy(i_colDidy, 2) = i_colDidy;
            i_colDidy = i_colDidy + 1;
        case 2
            N_colDimor(i_colDimor, 1) = data_c(i,3);
            N_colDimor(i_colDimor, 2) = i_colDimor;
            i_colDimor = i_colDimor + 1;
        case 3
            N_escape(i_escape, 1) = data_c(i,3);
            N_escape(i_escape, 2) = i_escape;
            i_escape = i_escape + 1;
    end
end



