function zPrev = so_standard_map_inverse_cylinder(z, cfg)
%SO_STANDARD_MAP_INVERSE_CYLINDER Inverse with x wrapped only at output.
zPrev = so_to_cylinder(so_standard_map_inverse_lifted(z, cfg));
end

