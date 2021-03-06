# ElectrolyteComposition

ElectrolyteComposition converts all electrolytes to a fixed string representation in solvent mass fraction and salt molality.
This is called a CompositionID'and looks like the following, with delimiter ''_'' and separator ''|''

$DELIMITED_SOLVENTS | $DELIMITED_SOLVENT_MASS_FRACS | $DELIMITED_SALTS | $DELIMITED_MOLALITIES

Thus, a DMC:EC 50:50 by wt solution of 2 molality LiPF6 looks like:
DMC_EC|50_50|LiPF6|2

Conversion happens with two specified precisions: 
1) solvent mass fraction precision (set in 100,1000, etc for whole percentages and single decimal percentages respectively)
2) salt molality precision (set in decimal places, i.e. 2 yields 2.00 for a 2 molality solution)

All constituent species must be added to the /data/ folder in either solventDB.csv or saltDB.csv - molar mass is required for conversions.

Three constructors are provided:
1) by_CompositionID
2) by_solution_volume
3) by_mole_fraction
4) by_mass_fraction_and_molality

See worked through example in demo_electrolytecomp.py

# Experiment

Class ''Clio'' with Experiment.py handles all I/O to the Clio web-server located on the lab CPU - 
LabView VIs translate these HTTP requests into orchestration commands.

An object is instantiated with specification of a few kwargs:
1) run_type: production, testing, etc. - this is passed along in the postdata and is used to segregate data downstream.
2) mix_volume: overall sample testing volume
3) person: a pre-fix for the RunID to specify who ran the script, used to segregate data downtream.

# vector_make_lookup

A script that creates grids of volumes across a set offeeder solutions,
then converts the grid to composition axes in a fast, vectorized manner.

Key parameters:
1) ''all_solvents'' : list of solvents present in ''densities_all'' for separation/calculation
2) ''salts'' : dict of salt species as key and molar mass as value
3) ''densities_all'' : dictionary of feeder solutions, keys as CompositionIDs, values as densities of the solutions
4) ''volume_step'', ''total'' : discretization of volume grid and total mix volume respectively.
