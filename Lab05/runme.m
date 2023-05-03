steps = [7];
clustername = oshostname;

%Cluster parameters
%cluster=generic('name',oshostname(),'np',15);
%cluster.executionpath='/local/helene/ENGS199/execution';
execution_dir='/jumbo/ice/ISSM/issm_anna/Lab04/Lab4Files/execution';
cluster=generic('np',5,'executionpath',execution_dir);


org=organizer('repository',['./Models'],'prefix','Kjer','steps',steps); clear steps;

if perform(org,'Mesh'),% {{{  %1

	md=triangle(model,['Exp/Domain_TEMP.exp'],500);
	
	load('/local/ModelData/GreenlandVelocity/velocity.mat');
	velx = InterpFromGridToMesh(x_m,y_m,double(vx),md.mesh.x,md.mesh.y,NaN);
	vely = InterpFromGridToMesh(x_m,y_m,double(vy),md.mesh.x,md.mesh.y,NaN);
	vel  = sqrt(velx.^2+vely.^2);

	%refine mesh using surface velocities as metric
	md=bamg(md,'hmin',100,'hmax',1100,'field',vel,'err',5);
	[md.mesh.lat,md.mesh.long]  = xy2ll(md.mesh.x,md.mesh.y,+1,45,70);
	md.mesh.epsg=3413;

	savemodel(org,md);
end %}}}
if perform(org,'Param'),% {{{     % 2

	md=loadmodel(org,'Mesh');
	md=setflowequation(md,'SSA','all');
	
	md=setmask(md,'','');
	md=parameterize(md,'Kjer.par');

	savemodel(org,md);
end%}}}
if perform(org,'Inversion_drag'),% {{{   % 3

	md=loadmodel(org,'Param');

	%Control general
	md.inversion=m1qn3inversion(md.inversion);
	md.inversion.iscontrol=1;
	md.verbose=verbose('solution',false,'control',true);
	md.transient.amr_frequency = 0;

	%Cost functions
	md.inversion.cost_functions=[101 103 501];
	md.inversion.cost_functions_coefficients=zeros(md.mesh.numberofvertices,numel(md.inversion.cost_functions));
	md.inversion.cost_functions_coefficients(:,1)=5000;  %2000
	md.inversion.cost_functions_coefficients(:,2)=10;
	md.inversion.cost_functions_coefficients(:,3)=.2*50^-3; %-3
	pos=find(md.mask.ice_levelset>0);
	md.inversion.cost_functions_coefficients(pos,1:2)=0;
	
	%Controls
	md.inversion.control_parameters={'FrictionCoefficient'};
	md.inversion.maxsteps=100;
	md.inversion.maxiter =100;
	md.inversion.min_parameters=0.05*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters=400*ones(md.mesh.numberofvertices,1);
	md.inversion.control_scaling_factors=1;

	%Additional parameters
	md.stressbalance.restol=0.01;
	md.stressbalance.reltol=0.1;
	md.stressbalance.abstol=NaN;

	%Go solve
	md.cluster=cluster;
	md=solve(md,'sb');

	%Put results back into the model
	md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;
	md.initialization.vx=md.results.StressbalanceSolution.Vx;
	md.initialization.vy=md.results.StressbalanceSolution.Vy;

	savemodel(org,md);
end%}}}
if perform(org,'Transient_Prep'),% {{{  %4

	md=loadmodel(org,'Inversion_drag');

	md.initialization.pressure = zeros(md.mesh.numberofvertices,1); %FIXME
	md.masstransport.spcthickness = NaN(md.mesh.numberofvertices,1); %FIXME
	md.initialization.temperature = (273.15-5)*ones(md.mesh.numberofvertices,1); %FIXME
	md.transient.amr_frequency = 0;
	md.settings.output_frequency=5;
	
	% Set parameters
	md.inversion.iscontrol=0;
	md.timestepping.start_time = 2007;
	md.timestepping.time_step  = 0.02;
	md.timestepping.final_time = 2017; 

	md.transient.ismovingfront=1;
	md.transient.isthermal=0;
	md.transient.isstressbalance=1;
	md.transient.ismasstransport=1;
	md.transient.isgroundingline=1;
	md.groundingline.migration = 'SubelementMigration';
	
	% spclevelset
	pos=find(md.mask.ice_levelset<0); md.mask.ice_levelset(pos)=-1;
	pos=find(md.mask.ice_levelset>0); md.mask.ice_levelset(pos)=+1;
	md.levelset.spclevelset = NaN(md.mesh.numberofvertices,1);
	pos = find(md.mesh.vertexonboundary);
	md.levelset.spclevelset(pos) = md.mask.ice_levelset(pos);
	md.levelset.reinit_frequency = 1;
	md.levelset.stabilization = 1;

	disp('Assigning qsg');
	timestamps = 1992:1/12:2016-1/12; % from 1992 to 101
	md.results.qsg = zeros(md.mesh.numberofvertices+1,numel(timestamps));
	md.results.qsg(end,:) =  timestamps;
	data=load('./Kjer_Forcing.mat');
	md.results.qsg(1:end-1,:) = repmat(data.qsg,md.mesh.numberofvertices,1); % m^3/day

	disp('Assigning TF');
	md.results.TF = zeros(md.mesh.numberofvertices+1,numel(timestamps));
	md.results.TF(end,:) =  timestamps;
	%tricky part here... find deepest depths
	[X Y]=meshgrid(min(md.mesh.x)-10e3:150:max(md.mesh.x)+10e3,min(md.mesh.y)-10e3:150:max(md.mesh.y)+10e3);
	B = interpBedmachineGreenland(X,Y,'bed','cubic','/local/ModelData/BedMachineGreenland/BedMachineGreenland-2022-03-17.nc');
	D = zeros(size(X));
	I = zeros(size(X));
	FAR = zeros(size(X)); FAR(:,1) = 1;
	for i=0:numel(data.depth)-1,
		thrsld = -data.depth(end-i);
		A=(B<thrsld);
		CC=bwlabel(A);
		pos=find(FAR& B<thrsld );
		list = unique(CC(pos));
		pos = find(ismember(CC,list) & D==0);
		D(pos) = -thrsld;
		I(pos) = numel(data.depth) -i;
	end
	%subplot(1,2,1); imagesc(D); subplot(1,2,2); imagesc(-B); caxis([0 max(data.depths)])
	md.results.TFdepths = InterpFromGridToMesh(X(1,:)',Y(:,1),D,md.mesh.x,md.mesh.y,NaN);
	md.results.TFindex  = InterpFromGridToMesh(X(1,:)',Y(:,1),I,md.mesh.x,md.mesh.y,NaN);

	for j=1:md.mesh.numberofvertices,
		if md.results.TFindex(j)>0,
			indices = 1:ceil(min(md.results.TFindex(j),size(data.TF,1)));

			%Maximum temperature
			md.results.TF(j,:) = max(data.TF(indices,:),[],1);

			%Depth-averaged temperature
			z=100:5:max(data.depth(indices));
			if indices==1
				Tz=0;%data.TF(indices,:);
			else
				Tz = interp1(data.depth(indices),data.TF(indices,:),z);
			end
			md.results.TF(j,:) = mean(Tz);
		end
	end
	pos=find(isnan(md.results.TF));
	md.results.TF(pos)=0;

	savemodel(org,md);
end%}}}

if perform(org,'Transient_calvingvonmises'),% {{{  %5

	md=loadmodel(org,'Transient_Prep');
	md.groundingline.friction_interpolation='SubelementFriction1';
	md.groundingline.melt_interpolation='NoMeltOnPartiallyFloating';

	% calving parameters
	md.calving = calvingvonmises();
    %only use grounded! 
	md.calving.stress_threshold_groundedice=0.75*10^6*ones(md.mesh.numberofvertices,1); %og 2.9*10^6, %% this is pretty good
	md.calving.stress_threshold_floatingice=150*10^3*ones(md.mesh.numberofvertices,1);  %og 150*10^3

	%Fix ice rise
	pos=find(ContourToNodes(md.mesh.x,md.mesh.y,'./Exp/nocalving.exp',2));
	md.calving.stress_threshold_groundedice(pos)=0.5*10^6;  %was 5*10^6
	md.calving.stress_threshold_floatingice(pos)=2*10^4; %was 2*10^4
	
	% meltingrate
	timestamps = md.results.TF(end,:);
	md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices+1,numel(timestamps));
	md.frontalforcings.meltingrate(end,:) = timestamps;

	A = 2.25e-4;
	B = 0.1125;
	Alpha = 0.39;
	Beta  = 1.18;
	BED = -repmat(min(0,md.geometry.bed),[1 numel(timestamps)]);
	TF  = max(0,md.results.TF(1:end-1,:));
	qsg = max(0,md.results.qsg(1:end-1,:));
	md.frontalforcings.meltingrate(1:end-1,:) = (A*BED.*qsg.^Alpha + B).*TF.^Beta;
	md.frontalforcings.meltingrate(1:end-1,:) = md.frontalforcings.meltingrate(1:end-1,:)*365*0.7; %Conversion from m/day to m/year

	disp('Extending time series by repeating year 2015');
	md.frontalforcings.meltingrate = [md.frontalforcings.meltingrate repmat(md.frontalforcings.meltingrate(:,end-11:end),[1 md.timestepping.final_time-2016]) ];
	md.frontalforcings.meltingrate(end,:) = [timestamps 2016:1/12:md.timestepping.final_time-1/12];
	
	md.cluster = cluster;
	md.miscellaneous.name = ['Kjer_vonMises_calibration_' num2str(mode(md.calving.stress_threshold_floatingice)./10^3) '_' num2str(mode(md.calving.stress_threshold_groundedice)./10^3)];
	
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','MaskOceanLevelset','MaskIceLevelset','CalvingCalvingrate','SigmaVM'};

	md.settings.waitonlock = 1;
	md.verbose.solution = 1;
	md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
	md=solve(md,'Transient');
	savemodel(org,md);

end%}}}
if perform(org,'Transient_hab'),% {{{   % 6 

	md=loadmodel(org,'Transient_Prep');
	md.groundingline.friction_interpolation='SubelementFriction1';
	md.groundingline.melt_interpolation='NoMeltOnPartiallyFloating';

	% calving parameters
	md.calving = calvinghab();
	md.calving.flotation_fraction = 9.16*10^-2*ones(md.mesh.numberofvertices,1);;  %7.6*10^-2 smaller = closer

	pos=find(ContourToNodes(md.mesh.x,md.mesh.y,'./Exp/nocalving.exp',2));
	md.calving.flotation_fraction(pos)=0;

	% meltingrate
	timestamps = md.results.TF(end,:);
	md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices+1,numel(timestamps));
	md.frontalforcings.meltingrate(end,:) = timestamps;

	A = 2.25e-4;
	B = 0.1125;
	Alpha = 0.39;
	Beta  = 1.18;
	BED = -repmat(min(0,md.geometry.bed),[1 numel(timestamps)]);
	TF  = max(0,md.results.TF(1:end-1,:));
	qsg = max(0,md.results.qsg(1:end-1,:));
	md.frontalforcings.meltingrate(1:end-1,:) = (A*BED.*qsg.^Alpha + B).*TF.^Beta;
	md.frontalforcings.meltingrate(1:end-1,:) = md.frontalforcings.meltingrate(1:end-1,:).*365*0.7; %Conversion from m/day to m/year
	
	disp('Extending time series by repeating year 2015');
	md.frontalforcings.meltingrate = [md.frontalforcings.meltingrate repmat(md.frontalforcings.meltingrate(:,end-11:end),[1 md.timestepping.final_time-2016]) ];
	md.frontalforcings.meltingrate(end,:) = [timestamps 2016:1/12:md.timestepping.final_time-1/12];
	
	md.cluster = cluster;
	md.miscellaneous.name = ['Kjer_hab_calibration_' num2str(max(md.calving.flotation_fraction).*10^3)];
	md.verbose.solution = 1;
	
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','MaskOceanLevelset','MaskIceLevelset'};

	md.settings.waitonlock = 1;
	md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
	md=solve(md,'Transient');

	savemodel(org,md);
end%}}}

if perform(org,'Plot'),% {{{  %7
	
	%md=loadmodel(org,'Transient_hab');  %(step 6) 
	md=loadmodel(org,'Transient_calvingvonmises');  % step 5 

	cmap = demmap(50,-700,700);
	plotmodel(md,'data',md.geometry.bed,'colormap',cmap,'caxis',[-700 700],'figure',1);
	hold on
	title('Observed ice front');

	mapa=flipud(hot(11));

	for ii=0:9
		expdisp(['./Exp/Kjer_' num2str(2007+ii) '.exp'],'linestyle',mapa(ii+1,:));
	end
	
	text(.01,0.99,'(a)','Fontsize',18,'Fontweight','bold','unit','normalized','VerticalAlignment','top','HorizontalAlignment','left');
	set(gca,'visible','on'); 

	plotmodel(md,'data',md.geometry.bed,'figure',2,'colormap',cmap,'caxis',[-700 700]);
	title('Modeled ice front');

	h0=isoline(md,md.mask.ice_levelset,'value',0,'output','matrix');
	hold on; plot(h0(:,1),h0(:,2),'-','Color',mapa(1,:));
	for ii=1:9
		hii=isoline(md,md.results.TransientSolution(ii*10).MaskIceLevelset,'value',0,'output','matrix');
		hold on; plot(hii(:,1),hii(:,2),'-','Color',mapa(ii+1,:));
	end
	
	text(.01,0.99,'(b)','Fontsize',18,'Fontweight','bold','unit','normalized','VerticalAlignment','top','HorizontalAlignment','left');
	set(gca,'visible','on'); 
	
end%}}}	


