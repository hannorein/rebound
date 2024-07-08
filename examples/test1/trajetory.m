% post-processing of rebound sim results
% Yun Zhang (2023 May)

clear all;

mass_system = 5.5e11;
mu_system = 6.6743e-11*5.5e11;
vol_didy = 0.2295409644951028e9;
vol_dimor = 0.001830603200702610e9;
mass_didy = mass_system*vol_didy/(vol_didy+vol_dimor);
mass_dimor = mass_system*vol_dimor/(vol_didy+vol_dimor);

[N_particle, time, N_colDidy, N_colDimor, N_escape, r_dust, particle, data_c, data_p] = ReadParticle('particles70000.txt','collide70000.txt');


% particles exclude escaper
N_particle_cor = N_particle;
for i = 1:length(N_escape)
    index_t = find(time>=N_escape(i,1));
    N_particle_cor(index_t(1):end) = N_particle_cor(index_t(1):end) + 1;
end

R_hill = 70000;  % Hill radius, m
time_hill = [];  % time reach the boundary of the hill radius, s
N_hill = [];  % particle id
ii = 1;
for i = 4:size(particle,1)
    for j = 1:size(particle(i).pos,1)
        dis_Didy = norm( particle(i).pos(j) - particle(1).pos(j) );
        if dis_Didy > R_hill
            time_hill(ii) = time(j);
            N_hill(ii) = particle(i).ID;
            ii = ii + 1;
            if particle(i).id_collide < 3
                length_col = length(data_c);
                data_c(length_col,1) = 3;
                data_c(length_col,1) = particle(i).ID;
                data_c(length_col,1) = time(j);
            end
            break;
        end
    end
end

time_sort = sort(time_hill);
N_sort = 1:1:length(time_sort);

% particles include new-defined escaper
N_particle_cor_esc = N_particle_cor;
for i = 1:length(N_sort)
    index_t = find(time>=time_sort(i));
    N_particle_cor_esc(index_t(1):end) = N_particle_cor_esc(index_t(1):end) - 1;
end

figure;
plot([0;N_colDidy(:,1);time(end)]/24/3600, [0;N_colDidy(:,2);N_colDidy(end,2)]/N_particle(1)*100, 'LineWidth',1.5, 'DisplayName','Didymos collider');
hold on;
plot([0;N_colDimor(:,1);time(end)]/24/3600, [0;N_colDimor(:,2);N_colDimor(end,2)]/N_particle(1)*100, 'LineWidth',1.5, 'DisplayName','Dimorphos collider');
plot([0;time_sort';time(end)]/24/3600, [0;N_sort';N_sort(end)]/N_particle(1)*100, 'LineWidth',1.5, 'DisplayName','Escaped ejecta');
plot(time/24/3600, N_particle_cor_esc/N_particle(1)*100, 'LineWidth',1.5, 'DisplayName','Remaining ejecta');
set(gca, 'FontName','Times','FontSize', 15);
set(gcf, 'Position',[100,500,400,300]);
xlabel('Time  [days]');
ylabel('Ejecta type percentage  [%]');
title('Dust particle radius  {\itr} = 1 mm');
legend;


figure; % Dimorphos orbit
plot3(particle(2).pos(:,1)-particle(1).pos(:,1), particle(2).pos(:,2)-particle(1).pos(:,2), particle(2).pos(:,3)-particle(1).pos(:,3), 'LineWidth',2.0,'Color','k');
hold on;
N_p = 200;
plot3(particle(N_p).pos(:,1)-particle(1).pos(1:size(particle(N_p).pos,1),1), ...
    particle(N_p).pos(:,2)-particle(1).pos(1:size(particle(N_p).pos,1),2),...
    particle(N_p).pos(:,3)-particle(1).pos(1:size(particle(N_p).pos,1),3));


figure;
plot([0;N_colDidy(:,1);time(end)]/24/3600, [0;N_colDidy(:,2);N_colDidy(end,2)]/N_particle(1)*100, 'LineWidth',1.5, 'DisplayName','Didymos collider');
hold on;
plot([0;N_colDimor(:,1);time(end)]/24/3600, [0;N_colDimor(:,2);N_colDimor(end,2)]/N_particle(1)*100, 'LineWidth',1.5, 'DisplayName','Dimorphos collider');
plot([0;N_escape(:,1);time(end)]/24/3600, [0;N_escape(:,2);N_escape(end,2)]/N_particle(1)*100, 'LineWidth',1.5, 'DisplayName','Escaped ejecta');
plot(time/24/3600, N_particle/N_particle(1)*100, 'LineWidth',1.5, 'DisplayName','Remaining ejecta');

set(gca, 'FontName','Times','FontSize', 15);
set(gcf, 'Position',[100,500,400,300]);
xlabel('Time  [days]');
ylabel('Ejecta type percentage  [%]');
title('Dust particle radius  {\itr} = 1 mm');
legend;


figure;
Didy_bary = [0.5*(vol_didy*data_p{1}(1,2:7) + vol_dimor*data_p{1}(2,2:7))/(vol_didy+vol_dimor)];
for i = 4:N_particle
    [coe(i-3,:),~]=r2e(data_p{1}(i,2:7)-Didy_bary,mu_system);
end

for i = 1:300
    for j = 1:20
        a_p(j+(i-1)*20) = 500.0 + (i-1)*20.0;
        e_p(j+(i-1)*20) = (j-1)/20;
    end
end

figure;
histogram(a_p(data_c(data_c(:,1)==1,2)-3),'BinWidth',200, 'DisplayName','Didymos collider');
hold on;histogram(a_p(data_c(data_c(:,1)==2,2)-3),'BinWidth',200, 'DisplayName','Dimorphos collider');
hold on;histogram(a_p(data_c(data_c(:,1)==3,2)-3),'BinWidth',200, 'DisplayName','Escaped ejecta');
hold on;histogram(a_p(data_p{end}(4:end,1)-3),'BinWidth',200, 'DisplayName','Remaining ejecta')
set(gca, 'FontName','Times','FontSize', 15);
set(gcf, 'Position',[100,500,400,300]);
xlabel('Semimajor axis  [m]');
ylabel('Number');
title('Dust particle radius  {\itr} = 1 mm');


figure;
histogram(e_p(data_c(data_c(:,1)==1,2)-3),'BinWidth',0.05, 'DisplayName','Didymos collider');
hold on;histogram(e_p(data_c(data_c(:,1)==2,2)-3),'BinWidth',0.05, 'DisplayName','Dimorphos collider');
hold on;histogram(e_p(data_c(data_c(:,1)==3,2)-3),'BinWidth',0.05, 'DisplayName','Escaped ejecta');
hold on;histogram(e_p(data_p{end}(4:end,1)-3),'BinWidth',0.05, 'DisplayName','Remaining ejecta')
set(gca, 'FontName','Times','FontSize', 15);
set(gcf, 'Position',[100,500,400,300]);
xlabel('Eccentricity');
ylabel('Number');
title('Dust particle radius  {\itr} = 1 mm');

figure;
color_code  = [particle(4:end).id_collide];
plot(e_p(color_code==0), a_p(color_code==0),'LineStyle','none','Marker','o','LineWidth',2.0);
hold on;
plot(e_p(color_code==1), a_p(color_code==1),'LineStyle','none','Marker','v','LineWidth',2.0);
plot(e_p(color_code==2), a_p(color_code==2),'LineStyle','none','Marker','p','LineWidth',2.0);
plot(e_p(color_code==3), a_p(color_code==3),'LineStyle','none','Marker','diamond','LineWidth',2.0);
title('Dust particle radius  {\itr} = 1 mm');

figure;
N_t = 12;
scatter3(data_p{N_t}(3:end,3), data_p{N_t}(3:end,4), data_p{N_t}(3:end,5), 10, data_p{N_t}(3:end,2),'filled');
