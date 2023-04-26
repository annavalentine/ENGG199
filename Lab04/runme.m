steps=[10 11];

%Set-up cluster parameter
execution_dir='/jumbo/ice/ISSM/issm_anna/Lab04/Lab4Files/execution';
cluster=generic('np',5,'executionpath',execution_dir);

%Run Steps
org=organizer('repository','./Models','prefix','MISMIP.','steps',steps);

%Initialization
if perform(org,'Mesh_generation'),% {{{1  %1

	md=model();
	md=bamg(md,'domain','./Domain.exp','hmax',2000,'splitcorners',1);
	md.miscellaneous.name=['MISMIP'];

	savemodel(org,md);
end % }}}
if perform(org,'Parameterization'),% {{{1   %2 

	md=loadmodel(org,'Mesh_generation');

	md=setmask(md,'','');
	md=parameterize(md,'./Mismip.par');
	md=setflowequation(md,'SSA','all');

	savemodel(org,md);
end% }}}
if perform(org,'Transient_Steadystate'),% {{{1  %3

	md=loadmodel(org,'Parameterization');

	md.timestepping.time_step=1;
	md.timestepping.final_time=200000;
	md.settings.output_frequency=2000;
	md.stressbalance.maxiter=10;
	md.stressbalance.abstol=NaN;
	md.stressbalance.restol=1;
	md.verbose=verbose('convergence',false,'solution',true,'module',false);
	md.cluster=cluster;
	md=solve(md,'tr');

	savemodel(org,md);
end% }}}

%Experiments
if perform(org,'Experiment0'),% {{{1   %4

	md=loadmodel(org,'Transient_Steadystate');

	md=setflowequation(md,'SSA','all');
	%Initialize model from last solution
	md.initialization.vx=md.results.TransientSolution(end).Vx;
	md.initialization.vy=md.results.TransientSolution(end).Vy;
	md.initialization.vel=md.results.TransientSolution(end).Vel;
	md.geometry.thickness=md.results.TransientSolution(end).Thickness;
	md.geometry.base=md.results.TransientSolution(end).Base;
	md.geometry.surface=md.results.TransientSolution(end).Surface;
	md.mask.ocean_levelset=md.results.TransientSolution(end).MaskOceanLevelset;

    md.friction = friction(md);
    md.friction.coefficient=sqrt(3.160*10^6)*ones(md.mesh.numberofvertices,1); 
    %md.friction.coefficientcoulomb=sqrt(0.5)*ones(md.mesh.numberofvertices,1); 
    md.friction.p=3*ones(md.mesh.numberofelements,1);
    md.friction.q=zeros(md.mesh.numberofelements,1);

	md.timestepping.time_step=0.25;
	md.timestepping.final_time=100;
	md.settings.output_frequency=4;
	md.stressbalance.maxiter=10;
	md.stressbalance.restol=1;
	md.stressbalance.reltol=0.001;
	md.stressbalance.abstol=NaN;
	md.verbose=verbose('convergence',false,'solution',true,'module',false);
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','BasalforcingsFloatingiceMeltingRate','GroundedArea'};
	md.cluster=cluster;
	md.miscellaneous.name=['MISMIP_Exp0'];
	md=solve(md,'tr');

	gl_pos=zeros(length(md.results.TransientSolution),1);
	for i=1:length(md.results.TransientSolution),
		[elements,x,y,z,s,data]=SectionValues(md,md.results.TransientSolution(i).MaskOceanLevelset,'SectionX0.exp',[100 10]);
		gl_pos(i) = x(max(find(data<0)));
	end

	savemodel(org,md);
end% }}}
if perform(org,'Experiment1'),% {{{1   %5 

	md=loadmodel(org,'Transient_Steadystate');

	md=setflowequation(md,'SSA','all');
	md.initialization.vx=md.results.TransientSolution(end).Vx;
	md.initialization.vy=md.results.TransientSolution(end).Vy;
	md.initialization.vel=md.results.TransientSolution(end).Vel;
	md.geometry.thickness=md.results.TransientSolution(end).Thickness;
	md.geometry.base=md.results.TransientSolution(end).Base;
	md.geometry.surface=md.results.TransientSolution(end).Surface;
	md.mask.ocean_levelset=md.results.TransientSolution(end).MaskOceanLevelset;

    md.friction = friction(md);
    md.friction.coefficient=sqrt(3.160*10^6)*ones(md.mesh.numberofvertices,1); 
    %md.friction.coefficientcoulomb=sqrt(0.5)*ones(md.mesh.numberofvertices,1); 
    md.friction.p=3*ones(md.mesh.numberofelements,1);
    md.friction.q=zeros(md.mesh.numberofelements,1);

	md.basalforcings.meltrate_factor=0.2;
	md.timestepping.time_step=0.25;
	md.timestepping.final_time=100;
	md.settings.output_frequency=4;
	md.stressbalance.maxiter=100;
	md.stressbalance.restol=1;
	md.stressbalance.reltol=0.001;
	md.stressbalance.abstol=NaN;
	md.verbose=verbose('convergence',false,'solution',true,'module',false);
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','BasalforcingsFloatingiceMeltingRate','GroundedArea'};
	md.miscellaneous.name=['MISMIP_Exp1'];
	md.cluster=cluster;
	md=solve(md,'tr');

	gl_pos=zeros(length(md.results.TransientSolution),1);
	for i=1:length(md.results.TransientSolution),
		[elements,x,y,z,s,data]=SectionValues(md,md.results.TransientSolution(i).MaskOceanLevelset,'SectionX0.exp',[100 10]);
		gl_pos(i) = x(max(find(data<0)));
    end

	savemodel(org,md);
end% }}}

% set the model friction to weertman with parameters 
if any(steps==7)
    
end


if any(steps==8)
    %md=loadmodel('./Models/MISMIP.Experiment1');
    %md=loadmodel('./Models/MISMIP.Experiment0');

    volume=[];
    for i=1:101;
    volume=[volume md.results.TransientSolution(i).IceVolume];
    end

    for i=1:length(md.results.TransientSolution),
		[elements,x,y,z,s,data]=SectionValues(md,md.results.TransientSolution(i).MaskOceanLevelset,'SectionX0.exp',[100 10]);
		gl_pos(i) = x(max(find(data<0)));
	end

    figure;
    subplot(2, 1, 1)
    plot([1:101], gl_pos)
    title('Grounding Line');
    xlabel('years');


    subplot(2,1,2);
    plot([1:101],volume);
    %Title this plot Mean Velocity and add an x label of years
    title('Ice Volume');
    xlabel('years');

end
if any(steps==9)   % change sliding law to weertman 
    md.friction = friction(md);
    md.friction.coefficient=sqrt(3.160*10^6)*ones(md.mesh.numberofvertices,1); 
    %md.friction.coefficientcoulomb=sqrt(0.5)*ones(md.mesh.numberofvertices,1); 
    md.friction.p=3*ones(md.mesh.numberofelements,1);
    md.friction.q=zeros(md.mesh.numberofelements,1);
end 
if any(steps==10)
    % calculate stress-field for steady-state of each model. 
    % weertman
    %make my equation easier
    C = sqrt(3.160*10^6); %md.friction.coefficient;
    q = 0; %md.friction.q;
    p = 3; %md.friction.p;
    s = 1/p;
    r = q/p;
    rho_ice = md.materials.rho_ice;
    g = md.constants.g;
    rho_water = md.materials.rho_water;
    thickness = md.geometry.thickness;
    bed = md.geometry.bed;
    Neff=rho_ice*g*thickness+rho_water*g*bed; 
    %N_eff = md.results.TransientSolution.Pressure;
    U_b = md.results.TransientSolution.Vel;
    basalStress_W= C.^2 .*(Neff.^r).*(abs(U_b).^(s-1)).* U_b;
    %%%%%%%%%%%%%%%%%%% now for coulomb
    C_C = sqrt(0.5);  %coulomb coeff
    basalStress_C = -min(basalStress_W, (C_C.^2)*Neff);
    
end

if any(steps==11)
    figure;
    %subplot(2, 1, 1)
    plotmodel(md, 'data', basalStress_W, 'title', 'Weertman', 'axis', 'tight','data', basalStress_C, 'title', 'Coulomb', 'axis', 'tight')
   
    
end