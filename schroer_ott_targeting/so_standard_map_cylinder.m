function zNext = so_standard_map_cylinder(z, cfg)
%SO_STANDARD_MAP_CYLINDER Cylinder representation: wrap x only.
zNext = so_to_cylinder(so_standard_map_lifted(z, cfg));
end

