steps=[1 2 3];

%Set-up cluster parameter
cluster=generic('np',2);
execution_dir=pwd; %add execution path where simulations are run
execution_dir='/thayerfs/home/f006fk7/ENGS199/execution'; %add execution pah where simulations are run
cluster.executionpath=execution_dir;
%Run Steps
org=organizer('repository','./Models','prefix','Thermal_','steps',steps);

if perform(org,'Mesh_generation'),% {{{

	md=model();
	%Use bamg to create a mesh
	md=bamg(model,'domain','Domain.exp','hmax',1000);  % domain path is working directory
	%add a name for your model
	md.miscellaneous.name='Lab1_Model_Valentine';

	savemodel(org,md);
end % }}}
if perform(org,'Parameterization'),% {{{

	md=loadmodel(org,'Mesh_generation');

	%Set the mask for gorunded ice only
	md=setmask(md,'','');    %% I am a little confused on this, I think I need

	disp('      creating model parameters');
	%Fill in information for all the missing parameters
	%Remember that most fields are defined on the vertices of the mesh
	md.geometry.surface= 300*ones(md.mesh.numberofvertices, 1); 
	md.geometry.thickness= 300*ones(md.mesh.numberofvertices,1);
	md.geometry.bed= md.geometry.surface - md.geometry.thickness; %% okay
	md.geometry.base= md.geometry.bed;%%%% okay
	md.initialization.vx= 0.*ones(md.mesh.numberofvertices, 1); 
	md.initialization.vy= 0*ones(md.mesh.numberofvertices, 1); 
	md.initialization.vz= 0*ones(md.mesh.numberofvertices, 1); 
	md.initialization.vel= 0*ones(md.mesh.numberofvertices, 1); 
	md.initialization.temperature= 230*ones(md.mesh.numberofvertices, 1); 
	md.materials.rheology_law='None';
	md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);
	md.materials.rheology_B=1.e-10*ones(md.mesh.numberofvertices,1);
	md.friction.coefficient=0.256.*ones(md.mesh.numberofvertices,1); % I should look up what a reasonable value is--> 
	md.friction.p=1*ones(md.mesh.numberofelements,1);
	md.friction.q=1*ones(md.mesh.numberofelements,1);
	md=SetIceSheetBC(md);
 
	md.mask.ice_levelset= -1.*ones(md.mesh.numberofvertices,1);   %%%I think this number needs to be negative (ice present?) 
	md.mask.ocean_levelset= 1.*ones(md.mesh.numberofvertices,1); %%% I think this needs to be positive (ice is grounded?) 
	md.basalforcings.geothermalflux= 0.0086.*ones(md.mesh.numberofvertices,1); % GHF in W/m^2 , noV
	md.basalforcings.groundedice_melting_rate= 0.*ones(md.mesh.numberofvertices,1) ; % noV
	md.basalforcings.floatingice_melting_rate= 0.*ones(md.mesh.numberofvertices,1) ; %noV
    
    md.smb.mass_balance = 0.1;

	md=setflowequation(md,'SSA','all');
	%Extrude the model to create a 3d model
    md = extrude(md, 50, 1);

	%Use Dirichlet conditions at the surface
	pos=find(md.mesh.vertexonsurface)
    %md.thermal.spctemperature = nan(md.mesh.numberofvertices,1);  
	md.thermal.spctemperature(pos)= 273-30; % starting temp
    
	%md.thermal.spctemperature = [md.initialization.temperature; 1];

	savemodel(org,md);
end% }}}
if perform(org,'Steadystate'),% {{{

	md=loadmodel(org,'Parameterization');

	%Change the time step
	md.timestepping.time_step=0; % set equal to zero, this is a steady-state solution

	md.cluster=cluster;
	md.verbose=verbose('convergence',true,'solution',true,'module',true);
	
    %Solve the thermal steady-state
	md=solve(md,'th');

	savemodel(org,md);
end% }}}
