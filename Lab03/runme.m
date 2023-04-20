steps=[4];
% step 4 is checkconfined
interactive=0;
execution_dir='/jumbo/ice/ISSM/issm_anna/Lab03/LabFiles/execution'; %add execution path where simulations are run
cluster=generic('np',5,'executionpath',execution_dir);
cluster.interactive=interactive;

%Run Steps
org=organizer('repository','./Models','prefix','IceShelf.','steps',steps);

%Confined ice shelf
if perform(org,'Mesh_generation'),% {{{1

	md=triangle(model,'./Domain.exp',4000);

	savemodel(org,md);
end % }}}
if perform(org,'ParameterizationConfined'),% {{{1

	md=loadmodel(org,'Mesh_generation');
	md.miscellaneous.name='Confined';
	md=setmask(md,'','');

	disp('      creating thickness');
	md.geometry.bed=-100-abs(md.mesh.x)/1000;
	md.geometry.base=md.geometry.bed;
	md.geometry.surface=max(md.geometry.bed+500,10);
	di=md.materials.rho_ice/md.materials.rho_water;
	pos=find(-di/(1-di)*md.geometry.surface > md.geometry.bed);
	md.geometry.base(pos)=di/(di-1)*md.geometry.surface(pos);
	md.geometry.thickness=md.geometry.surface-md.geometry.base;  % changing thickness here! 
	md.mask.ocean_levelset=md.geometry.thickness+md.geometry.bed/di; 

	disp('      creating drag');
	md.friction.coefficient=sqrt(10^7)*ones(md.mesh.numberofvertices,1); %q=1.
	md.friction.p=3*ones(md.mesh.numberofelements,1);
	md.friction.q=zeros(md.mesh.numberofelements,1);

	disp('      creating flow law paramter');
	md.materials.rheology_B=0.5*1/((10^-25)^(1/3))*ones(md.mesh.numberofvertices,1);
	md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);
	md.materials.rheology_law='None';

	disp('      boundary conditions for stressbalance model');
	md=SetIceShelfBC(md,'./Front.exp');
	pos=find(md.mesh.x>750000);    %% this is what you change for surface!! 
	md.mask.ice_levelset(pos)=1;
	md.stressbalance.spcvx(:)=NaN;
	md.stressbalance.spcvy(:)=NaN;
	pos=find(md.mesh.y==50000 | md.mesh.y==0);
	md.stressbalance.spcvy(pos)=0;
	md.stressbalance.spcvx(pos)=0;
	pos2=find(md.mesh.x<0.1);
	md.stressbalance.spcvx(pos2)=0;

	disp('      forcing conditions');
	md.smb.mass_balance=0.5*ones(md.mesh.numberofvertices,1);
	md.basalforcings.geothermalflux=0.5*ones(md.mesh.numberofvertices,1);

	md.thermal.spctemperature=NaN*ones(md.mesh.numberofvertices,1);
	md.groundingline.migration='SubelementMigration';
	md.groundingline.friction_interpolation='SubelementFriction1'; 
	md.groundingline.melt_interpolation='NoMeltOnPartiallyFloating'; 

	%Parameters
	md.transient.isthermal=0;
	md.transient.isgroundingline=1;

	%Initialization
	md.initialization.vx=ones(md.mesh.numberofvertices,1);
	md.initialization.vy=ones(md.mesh.numberofvertices,1);
	md.initialization.vz=ones(md.mesh.numberofvertices,1);
	md.initialization.vel=sqrt(2)*ones(md.mesh.numberofvertices,1);
	md.initialization.temperature=273*ones(md.mesh.numberofvertices,1);

	md=setflowequation(md,'ssa','all');

	savemodel(org,md);
end% }}}
if perform(org,'SteadyStateConfined'),% {{{1

	md=loadmodel(org,'ParameterizationConfined');
	md.verbose=verbose('solution',true,'module',true,'convergence',false);
	md.timestepping.time_step=1;
	md.timestepping.final_time=50000;
	md.settings.output_frequency=100;
	md.cluster=cluster;
	md=solve(md,'tr');

	savemodel(org,md);
end% }}}
if perform(org,'CheckSteadyStateConfined'),% {{{1

	md=loadmodel(org,'SteadyStateConfined');
    

	%Put the results of the steady-state as initial conditions (check all the fields to be updated in the results)
	md.geometry.surface= md.results.TransientSolution(end).Surface;
	md.geometry.thickness= md.results.TransientSolution(end).Thickness;
    %change thickness just on the floating ice
    %pos=find(md.mask.ocean_levelset < 0);
	%md.geometry.thickness(pos)=md.geometry.thickness(pos)+1000;
	md.geometry.base= md.geometry.surface - md.geometry.thickness;
    %%%% Keep going on initialization woohoo! ( i need to look at what is
    md.initialization.vx= md.results.TransientSolution(end).Vx;
	md.initialization.vy= md.results.TransientSolution(end).Vy; 
	md.initialization.vel= md.results.TransientSolution(end).Vel;
	md.initialization.pressure = md.results.TransientSolution(end).Pressure;

    %ocean levelset
    md.mask.ocean_levelset = md.results.TransientSolution(end).MaskOceanLevelset;

    %smb
    md.smb.mass_balance = md.results.TransientSolution(end).SmbMassBalance;
    
    % change the front
    pos=find(md.mesh.x>725000);  % originally 750000
	md.mask.ice_levelset(pos)=1; % positive for more no-ice


	%Update cluster and solve
	md.verbose=verbose('solution',true,'module',true,'convergence',false);
	md.cluster=cluster;
	md=solve(md,'sb');

	savemodel(org,md);
end% }}}

%Unconfined ice shelf
if perform(org,'ParameterizationUnconfined'),% {{{1

	md=loadmodel(org,'Mesh_generation');
	md.miscellaneous.name='Unconfined';
	md=setmask(md,'','');

	disp('      creating thickness');
	md.geometry.bed=-100-abs(md.mesh.x)/1000;
	md.geometry.base=md.geometry.bed;
	md.geometry.surface=max(md.geometry.bed+500,10);
	di=md.materials.rho_ice/md.materials.rho_water;
	pos=find(-di/(1-di)*md.geometry.surface > md.geometry.bed);
	md.geometry.base(pos)=di/(di-1)*md.geometry.surface(pos);
	md.geometry.thickness=md.geometry.surface-md.geometry.base;  %% change thickness
	md.mask.ocean_levelset=md.geometry.thickness+md.geometry.bed/di;

	disp('      creating drag');
	md.friction.coefficient=sqrt(10^7)*ones(md.mesh.numberofvertices,1); %q=1.
	md.friction.p=3*ones(md.mesh.numberofelements,1);
	md.friction.q=zeros(md.mesh.numberofelements,1);

	disp('      creating flow law paramter');
	md.materials.rheology_B=0.5*1/((10^-25)^(1/3))*ones(md.mesh.numberofvertices,1);
	md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);
	md.materials.rheology_law='None';

	disp('      boundary conditions for stressbalance model');
	md=SetIceShelfBC(md,'./Front.exp');
	pos=find(md.mesh.x>750000);
	md.mask.ice_levelset(pos)=1;
	md.stressbalance.spcvx(:)=NaN;
	md.stressbalance.spcvy(:)=NaN;
	pos=find(md.mesh.y==50000 | md.mesh.y==0);
	md.stressbalance.spcvy(pos)=0;
	pos2=find(md.mesh.x<0.1);
	md.stressbalance.spcvx(pos2)=0;

	disp('      forcing conditions');
	md.smb.mass_balance=0.5*ones(md.mesh.numberofvertices,1);
	md.basalforcings.geothermalflux=0.5*ones(md.mesh.numberofvertices,1);

	md.thermal.spctemperature=NaN*ones(md.mesh.numberofvertices,1);
	md.groundingline.migration='SubelementMigration';
	md.groundingline.friction_interpolation='SubelementFriction1'; 
	md.groundingline.melt_interpolation='NoMeltOnPartiallyFloating'; 

	%Parameters
	md.transient.isthermal=0;
	md.transient.isgroundingline=1;

	%Initialization
	md.initialization.vx=ones(md.mesh.numberofvertices,1);
	md.initialization.vy=ones(md.mesh.numberofvertices,1);
	md.initialization.vz=ones(md.mesh.numberofvertices,1);
	md.initialization.vel=sqrt(2)*ones(md.mesh.numberofvertices,1);
	md.initialization.temperature=273*ones(md.mesh.numberofvertices,1);

	md=setflowequation(md,'ssa','all');

	savemodel(org,md);
end% }}}
if perform(org,'SteadyStateUnconfined'),% {{{1

	md=loadmodel(org,'ParameterizationUnconfined');
	md.verbose=verbose('solution',true,'module',true,'convergence',false);
	md.timestepping.time_step=1;
	md.timestepping.final_time=50000;
	md.settings.output_frequency=100;
	md.cluster=cluster;
	error
	md=solve(md,'tr');

	savemodel(org,md);
end% }}}
if perform(org,'CheckSteadyStateUnconfined'),% {{{1

	md=loadmodel(org,'SteadyStateUnconfined');

	%Put the results of the steady-state as initial conditions (check all the fields to be updated in the results)
	md.geometry.surface= md.results.TransientSolution(end).Surface;
	md.geometry.thickness= md.results.TransientSolution(end).Thickness;
    %change thickness just on the floating ice
    %pos=find(md.mask.ocean_levelset < 0);
	%md.geometry.thickness(pos)=md.geometry.thickness(pos)+500;
	md.geometry.base= md.geometry.surface - md.geometry.thickness;
    md.initialization.vx= md.results.TransientSolution(end).Vx;
	md.initialization.vy= md.results.TransientSolution(end).Vy; 
	md.initialization.vel= md.results.TransientSolution(end).Vel;
	md.initialization.pressure = md.results.TransientSolution(end).Pressure;

    %ocean levelset
    md.mask.ocean_levelset = md.results.TransientSolution(end).MaskOceanLevelset;

    %smb
    md.smb.mass_balance = md.results.TransientSolution(end).SmbMassBalance;

    % change the front
    pos=find(md.mesh.x>700000);  % originally 750000
	md.mask.ice_levelset(pos)=1; % positive for more no-ice


	%Update cluster and solve
	md.verbose=verbose('solution',true,'module',true,'convergence',false);
	md.cluster=cluster;
	md=solve(md,'sb');

	savemodel(org,md);
end% }}}
