function [reduced_fabric,reduced_sample] = reduceFB(fabric,sample,new_n_machines)

N_machines = fabric.numFabricPorts/2;

reduced_sample = sample;

reduced_fabric = network_elements.FabricSkeleton(new_n_machines);
for m = 1:new_n_machines
    reduced_fabric.setIngress(m,1);
    reduced_fabric.setEgress(m,1);
end

for c = reduced_sample
    for f = c.flows
        f.source.id = mod(f.source.id - 1,new_n_machines) + 1;
        f.destination.id = mod(f.destination.id - N_machines - 1,new_n_machines) + new_n_machines + 1;
        f.links = [f.source.id,f.destination.id];
    end
    c.update(new_n_machines);
end

end

