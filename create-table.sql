create table if not exists simulations
       (lambda real not null
       ,number_of_sites integer not null
       ,energy_gap real not null
       ,multisweep_convergence_criterion real not null
       ,bandwidth_increase_convergence_criterion real not null
       ,timestamp text not null
       );
create unique index if not exists simulation_parameters on simulations (lambda,number_of_sites);
