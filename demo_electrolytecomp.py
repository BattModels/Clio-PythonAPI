from ElectrolyteComposition import ElectrolyteComposition

#Mixing volumes of solutions
rho_feeders={
	"DMC|100|LiTFSI|2":1.321,
	"DMC_EC|50_50|LiTFSI|2":1.436,
	"DMC|100":1.071,
	"DMC_EC|50_50":1.192,
}
sample_recipe = {
	"DMC|100|LiTFSI|2":500,
	"DMC_EC|50_50|LiTFSI|2":500,
	"DMC|100":500,
	"DMC_EC|50_50":500,
}
el = ElectrolyteComposition.by_solution_volume(volumes = sample_recipe, densities = rho_feeders)
print(el)
el = ElectrolyteComposition.by_solution_volume(volumes = sample_recipe, densities = rho_feeders,solvent_precision=1000)
print(el)
el = ElectrolyteComposition.by_solution_volume(volumes = sample_recipe, densities = rho_feeders,solvent_precision=1000,salt_decimals=3)
print(el)

#Specifying via CompositionID
el = ElectrolyteComposition.by_CompositionID("DMC_EC|74_26|LiTFSI|0.87")
print(el)

#Specifying via mass fraction and molality
solvents_by_mf = {
	"DMC":74,
	"EC":26,
}
salts_by_mol = {
	"LiTFSI":0.87
}
el = ElectrolyteComposition.by_mass_fraction_and_molality(solvents = solvents_by_mf, salts = salts_by_mol)
print(el)

