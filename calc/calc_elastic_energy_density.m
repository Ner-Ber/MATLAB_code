function E=calc_elastic_energy_density(struct)


%E=(-struct.Sxx.*struct.Uxx-struct.Syy.*struct.Uyy+struct.Sxy.*struct.Uxy)/2; %compresion in stress is negative-> '-' is needed
E=(-struct.Sxx.*struct.Uxx+struct.Sxy.*struct.Uxy)/2; %compresion in stress is negative-> '-' is needed