% modification date:  6/16/19, cylindrical puff target
run_name = mfilename;
 
%  data(particle#,property)
% INITIAL CONDITION CALCULATIONS
% NOTES

%use ray tracing to include ohmic heating on side wall of thruster

%% =========================common changes to the input file moved up to the top==============================================
eV = 11605; 
ma = 1.6605e-27;
kg = 1;
grams = .001*kg;

% target parameters
target_outer_radius = 0.02;
target_inner_radius = 0;
target_length = .02;
hsize = (target_outer_radius - target_inner_radius)/5;

%initial conditions
Target_temperature = 100*eV;
Target_mass = 10*grams;

%final conditions
nozzle_radius = 1;  %making this up
Final_time = 3*nozzle_radius/(sqrt(1.67 * 8814.5/235 * Target_temperature));

% %n = 3.17e22;  %initial ion density
% P0 = 1e7;  %%stagnation pressure in plenum
% MW = 1; %molecular weight of gas
% T0 = 2500; %K chamber temperature
% %placeholder for mass flow rate, or 1D assumption about mass flow rate into the engine
% 
% Npts = 75;  %resolution
% hsize = 2/Npts;
% 
% final_time = 1e-3; 
% output_time_step = final_time/100;

%% derived inputs
% R = 8314.5/MW;
Target_volume = 4/3*pi*target_outer_radius^3;
Target_density = Target_mass/Target_volume;


%turn on desired physics models
% physics.two_temperature = 1;  %split electron and ion temperatures, uses energy_e and energy_i instead of energy.m
% physics.conductivity = 1;  %electrical conductivity  
% physics.ionization = 1;  %ionization state, number of free electrons per ion
% physics.radiation = 1;
% physics.radiation_model = 'optically_thin';  %or thin
% physics.thermal_conduction=1;
% physics.ion_slip = 0;
%physics.radiation_model = 'thin_or_thick';  %or thin
%physics.radiation_model = 'optically_thick';  %or thin
%solver.timer_radiation_start = 21e-6;  %for turning on radiation later

%% -----------Creation of electrode and target geometry

%create a structured array (i use 'input') with fields length, hsize, xc, yc, zc, outer_radius, inner_radius
inputs(1).function_handle = @spherical_shell;  %target

inputs(1).outer_radius = target_outer_radius;
inputs(1).inner_radius = 0;  %one or both of the radii can be a constant value
inputs(1).hsize = hsize;
inputs(1).xc = 0;
inputs(1).yc = 0;
inputs(1).zc = 0;
inputs(1).fig = 1; %figure number that the plot is generated in
inputs(1).color = [0 119 200]/256;  %charger blue
inputs(1).hold = 'on';  %allows you to add other blocks to the same plot
inputs(1).name = 'target';

%make sure the cylindrical shell function is in your path
addpath([repmat(['..' filesep],1,sum(strfind(pwd,filesep)>strfind(pwd,'SPFMax'))) 'geometry']); 

%call the relevant function using this script below
for ii = 1:numel(inputs)
    %    [p, t, hsize, volume] = cylindrical_shell(inputs(ii));
    [p, t, hsize, volume] = feval(inputs(ii).function_handle,inputs(ii));
    h = hsize;
    N = size(p,1);
    x = p(:,1); y = p(:,2); z = p(:,3);
    inputs(ii).name = [inputs(ii).name num2str(N) '.mat'];
    save(inputs(ii).name,'h','N','p','x','y','z','volume');
    inputs(ii).name = [pwd '/' inputs(ii).name];
    
    %set the geometry file name for each of the items (i.e. the computational blocks)
    item(ii).inputs = inputs(ii);
    item(ii).type = 'matfile';
    item(ii).filename = item(ii).inputs.name; %sphere
end

%% parameters describing different materials, gases, etc (i.e. the thermodynamic and transport properties and models of each material)

% 1 - lithium/lithium deuteride
matter(1).name = 'Li-LiD';  %homogeneous approximation of lithium with LiD drops suspended in it
matter(1).relative_permittivity = 1; 
matter(1).relative_permeability = 1;
matter(1).conductivity = 1.078e7;  %1/ohm-m
matter(1).thermal_conductivity = 84.8; % [W/mK]
matter(1).MW = (6*4+4)/5;  %kg/kmol
matter(1).g = 1.67;
matter(1).opacity_Planck_abs = 'Bremsstrahlung';
matter(1).Z = 2;
matter(1).hydrostatic = 'stationary';  %keep it still for now
matter(1).species.XD = 1/(6+2);  %mole fractions of each species.  they all default to zero unless specified in the input file.  The format is 'X' with the species appended.  Must be consistent with physics.species_fields
matter(1).species.XLi6 = 1-matter(1).species.XD;

% 2 - Lead
matter(2).name = 'Pb';
matter(2).phase = 'solid';
matter(2).eos = 'solid';
matter(2).hydrostatic = 'stationary';
matter(2).thermal_conductivity = 50.2;
matter(2).specific_heat = 490; %J/kgK
matter(2).density = 7600;
matter(2).T = 293.15;   % Initial temperature of lead
matter(2).species.XPb = 1;  %mole fraction of lead

% 3 Uranium
matter(3).name = 'U238';  %uranium
%matter(3).phase = 'ablation';
%matter(3).T_melting = 1356.15+720.235;
%matter(3).T_boiling = 2840.15+720.235;
%matter(3).relative_permittivity = 1;
%matter(3).relative_permeability = 1;
matter(3).eos = 'solid';    % EOS Solid change
matter(3).conductivity = 3.57e7 ;  %1/ohm-m
matter(3).hydrostatic = 'stationary';
matter(3).thermal_conductivity = 27.5;   %W/m-K
matter(3).density = 19.1e3 * 3; %shocked to this density, but not by charger because it would take more than 2 MA for a 5 mm outer radius
matter(3).specific_heat = 116;
matter(3).MW = 238;
matter(3).opacity_Planck_abs = 'nist';  %hey mitchell, this sets the absorption coefficient, but matter(ii).name has to match the Xe, W, Water, etc material you input, as it does above for Cu
matter(3).species.XU238 = 1;  %mole fraction of uranium
matter(3).Z = 1;

% 4 Tungsten
matter(4).name = 'W';  %Tungsten
%matter(4).phase = 'ablation';
%matter(4).T_melting = 1356.15+720.235;
%matter(4).T_boiling = 2840.15+720.235;
%matter(4).relative_permittivity = 1;
%matter(4).relative_permeability = 1;
matter(4).eos = 'solid';    % EOS Solid change
matter(4).conductivity = 1/52.8e-9 ;  %1/ohm-m
matter(4).hydrostatic = 'stationary';
matter(4).thermal_conductivity = 173;   %W/m-K
matter(4).density = 19.3e3 ; 
matter(4).specific_heat = 134;
matter(4).MW = 183.84;
matter(4).opacity_Planck_abs = 'nist';  %hey mitchell, this sets the absorption coefficient, but matter(ii).name has to match the Xe, W, Water, etc material you input, as it does above for Cu

% 1 - lithium deuteride
matter(5).name = 'LiD';  %homogeneous approximation of lithium with LiD drops suspended in it
matter(5).relative_permittivity = 1; 
matter(5).relative_permeability = 1;
matter(5).conductivity = 1.078e7;  %1/ohm-m
matter(5).thermal_conductivity = 84.8; % [W/mK]
matter(5).MW = (2 + 6)/2;  %kg/kmol
matter(5).g = 1.67;
matter(5).opacity_Planck_abs = 'Bremsstrahlung';
matter(5).Z = 2;
matter(5).hydrostatic = 'stationary';  %keep it still for now
matter(5).species.XD = .5;  %mole fractions of each species.  they all default to zero unless specified in the input file.  The format is 'X' with the species appended.  Must be consistent with physics.species_fields
matter(5).species.XLi6 = .5;

% Plutonium
matter(6).name = 'Pu';  %good stuff
matter(6).relative_permittivity = 1; 
matter(6).relative_permeability = 1;
matter(6).conductivity = 1/1.46e-6;  %1/ohm-m
matter(6).thermal_conductivity = 6.74; % [W/mK]
matter(6).MW = 239;  %kg/kmol
matter(6).g = 1.67;
matter(6).opacity_Planck_abs = 'nist';
matter(6).Z = 1;
matter(6).hydrostatic = 'stationary';  %keep it still for now
matter(6).species.XPu = 1;  %mole fractions of each species.  they all default to zero unless specified in the input file.  The format is 'X' with the species appended.  Must be consistent with physics.species_fields
matter(6).T_boiling = 3505;
matter(6).bulk_modulus = 100e9;
matter(6).solid_density = 19816;
matter(6).eos = 'idealgas';

% 7 Uranium 235
matter(7).name = 'U235';  %uranium
%matter(7).phase = 'ablation';
%matter(7).T_melting = 1356.15+720.235;
%matter(7).T_boiling = 2840.15+720.235;
%matter(7).relative_permittivity = 1;
%matter(7).relative_permeability = 1;
% matter(7).eos = 'solid';    % EOS Solid change
% matter(7).conductivity = 3.57e7 ;  %1/ohm-m
% %matter(7).hydrostatic = 'stationary';
% matter(7).thermal_conductivity = 27.5;   %W/m-K
% matter(7).density = 19.1e3 * 3; %shocked to this density, but not by charger because it would take more than 2 MA for a 5 mm outer radius
% matter(7).specific_heat = 116;
% matter(7).MW = 235;
% matter(7).opacity_Planck_abs = 'nist';  %hey mitchell, this sets the absorption coefficient, but matter(ii).name has to match the Xe, W, Water, etc material you input, as it does above for Cu
% matter(7).species.XU235 = 1;  %mole fraction of uranium
% matter(7).Z = 1;
matter(7).relative_permittivity = 1; 
matter(7).relative_permeability = 1;
matter(7).conductivity = 3.57e7;  %1/ohm-m
matter(7).thermal_conductivity = 27.5; % [W/mK]
matter(7).MW = 235;  %kg/kmol
matter(7).g = 1.67;
matter(7).opacity_Planck_abs = 'nist';
matter(7).Z = 1;
%matter(7).hydrostatic = 'stationary';  %keep it still for now
matter(7).species.XU235 = 1;  %mole fractions of each species.  they all default to zero unless specified in the input file.  The format is 'X' with the species appended.  Must be consistent with physics.species_fields
matter(7).T_boiling = 3505;
matter(7).bulk_modulus = 100e9;
matter(7).solid_density = 18900;
matter(7).eos = 'idealgas';

matter(8).name = 'Cu';  %copper
matter(8).phase = 'ablation';
matter(8).T_melting = 1356.15+720.235;
matter(8).T_boiling = 2840.15+720.235;
matter(8).relative_permittivity = 1;
matter(8).relative_permeability = 1;
matter(8).eos = 'solid';    % EOS Solid change
matter(8).conductivity = 5.95e7 ;  %1/ohm-m
matter(8).hydrostatic = 'stationary';
matter(8).thermal_conductivity = 401;   %IF THIS SHOWS 1, REMEMBER TO CHANGE THIS LATER
matter(8).density = 8960;
matter(8).specific_heat = 386;


%% item definition, where each item is a block with particle geometry, material, and initial conditions specified 
% define the item, initial conditions, state, and particle resolution
item(1).matter = 'U235';
item(1).T = Target_temperature;  
item(1).rho = Target_density;   
% 
% item(2).matter = 'Cu';
% item(2).T = 300;  
% item(2).rho = 8960;   
% item(2).isstatic = 1;
% 
% item(3).matter = 'Cu';
% item(3).T = 300;  
% item(3).rho = 8960;   
% item(3).isstatic = 1;

%% physics models (except for nuclear because it is more convenient to group them together elsewhere)
%physics.thermal_conduction = 1;
%physics.phase = 1;
physics.conductivity = 0;  %needed for circuit model
physics.circuit = 0;  %needed for circuit model
%physics.current_density = 1;
%physics.update_current_manually = 1;
physics.rays = 1;  %ray tracing
physics.radiation = 1;  %turn on radiation
physics.radiation_model = 'Bremsstrahlung'; %specify the type of radiation model
physics.compressibility_correction = 0;  %turn on compressibility to account for compression of solids, needs a lot of inputs (solid_density, bulk_modulus, boiling temperature, all set in matter, and can be material dependent)


%setup a block of rays that are distributed according to some convenient method.  see initialize_rays for all the options
rayblock.type = 'snap2particles';  %snaps a number of rays to particles in a way that is hopefully evenly spaced. 
rayblock.cycle_ray_positions = 1;  %for snap2particles, this cycles the ray positions every hydrostep, by snapping to different particles each step or every few steps.
rayblock.fraction_of_particles = .001;  %1/(Njets + target)
rayblock.raytype = '4pi';
%rayblock.Rc = .08;  %multiplier times 'h' to model how far away from particle source the radiation is modeled
rayblock.nh = 20;  %multiplier times 'h' to model how far away from particle source the radiation is modeled for when ray changes size
rayblock.nR = 20; %
rayblock.nPol = 8;
rayblock.nAz = 16;  %2x as many as nPol gives even angular resolution in polar and azimuthal directions
rayblock.use_ion = 1; %for attenuation of fast ions or other particles not photons
rayblock.use_neutron = 0; %probably not going to be used, 
rayblock.use_em_radiation = 0;  %for attenuation of radiation


%% Fusion and fission parameters
% For simulations of fusion reactions, fission reactions, and scattering with homogeneous mixtures of species, ormore than 1 material which may react, need mole fractions, requiring user input for matter(ii).species.Xabc and physics.species_fields{mm} = 'abc';  Not literally abc, but a fusion/fission species.
%
%cell array of all species present.  We only need the unique list, no duplicates.
% Rules:
% 1. The string needs to be consistent with a name left of X in each of matter(xx).species.XPb
% 2. set matter(nn).species.XLi or other (replacing 'Li6' with D, Pb, etc) to the appropriate mole fraction.  Any species listed below in physics.species_field{mm} not in included in this matter(nn) will be set to 0. 
%
% physics.species_fields and matter.species.Xabc are used in define_particles to set the particle mole fractions for each species
%                           'string here'
%physics.species_fields{1} = 'Li6'; %Lithium
%physics.species_fields{2} = 'D';  %Deuterium
%physics.species_fields{3} = 'Pb'; %Lead
%physics.species_fields{3} = 'U235';  %Uranium
physics.species_fields{1} = item(1).matter;  %Uranium

physics.fusion = 0; %turn on fusion reactions
physics.fission = 1; %turn on fission reactions with neutrons
physics.neutrons = 1;  %neutron diffusion model set to 'on'.  will be set automatically to 'on' in 'define_fission_reactions, but go ahead, turn it on here too if you like.  whatevs.
physics.neutron_source = 1;  %this sets a controlled neutron source rate
physics.neutron_source_rate = 1e10;  %cannned neutron source rate.

%thermonuclear fusion reaction
physics.fusion_thermal = 'DD'; %specifies the thermonuclear reaction mechanism, built in 'define_fusion_reactions'
%stopping power for ion deposition
physics.stopping_power =  'HarrisMiley'; %set a value, or string, one of Bethe, HarrisMiley, LiPetrasso

%% Output and solver conditions 
solver.final_time = Final_time;
output.time_step = solver.final_time/100; % 1e-7;

% solver.gpu = 1;  %setting this to 0 uses the cpu instead.
% solver.gpuidx = 0;  %this can be chosen automatically, but you can force it to a specific gpu device if there is more than one.

%% ==============  Electromagnetic field and circuit options
% constant flux, or flux is controlled over time.  Basically it calculates a scalar magnetic field for affecting pressure 
% and transport quantities without the fuss of solving maxwell's equations. 
% you do lose a lot of physics but it is useful for exploring possibilities in 3D if B is given.  Kevin Schillo's dissertation is based on it.

%------------------Setting constant flux, which does not use circuit or electromagnetic field solver
% solver.magnetic_flux = 1;  %must set this
% item(1).Bx = 1;  % basically, you need to set up the magnetic field in the items you want, and the orientation does not matter
% item(1).By = 0;  % you could set it up as By, Bz, Bx, or any combination.  It's the magnitude that matters
% item(1).Bz = 0;  % With scalar B initiated as sqrt(Bx.^2 + By^2 + Bz^2), the intial flux is computed and is preserved, or programmed to decay

%-----------Use solver, phyiscs, and circuit inputs below to use the circuit and electromagnetic equation solvers
solver.Maxwell_equations = 0;

% %rays are now used to collect charge density and current density sources for propagating field
% physics.rays = 1;  %ray tracing  

physics.circuit = 0;
circuit(1).type = 'rlc'; %'rlc';
circuit(1).positive_lead_block = 2;  %anode, the integer refers to the corresponding item number above
circuit(1).positive_lead_side  = 'zmin';  %could be faces of brick (xmin, xmax, ymin, etc) or some custom thing
circuit(1).negative_lead_block = 3;  %cathode, the integer refers to the corresponding item number above
circuit(1).negative_lead_side  = 'zmax';
circuit(1).current = 0; %amps %constant if constant model, otherwise integrated 
circuit(1).voltage = 1e7;  %volts
circuit(1).R = 1e-3;
circuit(1).L = 50e-9;
circuit(1).C = 50e-6;

%unit normal for current going into/out of electrode
circuit(1).pos_nx = 0;  
circuit(1).pos_ny = 0;
circuit(1).pos_nz = 1; %+1 means current flows in positive z direction
circuit(1).neg_nx = 0;
circuit(1).neg_ny = 0;
circuit(1).neg_nz = 1; %+1 current flows in positive z direction

%bounding box for selecting coil particles connected to circuit
% circuit(1).pos_x0 = -4;  circuit(1).pos_dx = 8;
% circuit(1).pos_y0 = 0.1; circuit(1).pos_dy = 10
% circuit(1).pos_z0 = -50; circuit(1).pos_dz = 6;
% 
% circuit(1).neg_x0 = -4;   circuit(1).neg_dx = 8;
% circuit(1).neg_y0 =-10.1; circuit(1).neg_dy = 10;
% circuit(1).neg_z0 = -50;  circuit(1).neg_dz = 6;


