%% Misorientations for a given phase

% set phase of interest
phase = 'f';


% minimum and maximum angles (degrees) of misorientations to consider...
minAng = 2;     % minimum
maxAng = 10;    % maximum

% duplicate ebsd dataset for analysis
eLowAngle = ebsd;

% compute boundaries with minimum segmentation angle (defined above)
[gLowAngle,eLowAngle.grainId,eLowAngle.mis2mean] = calcGrains(eLowAngle('indexed'),'angle',minAng*degree,'tight');

% check if sample name variable exists (for plotting later)
if exist('sampleName')<1
sampleName = 'test';
end

%% Boundaries

% All boundaries from the grainset
gB = gLowAngle.boundary;

% Outer boundaries for the specified phase
gB_phase = gB(phase,phase);

% Inner boundaries for the specified phase
innerB = gLowAngle.innerBoundary(phase,phase);

% Concatenate all boundaries (inner and outer)
all_B = [gB_phase;innerB];


%% Compute the misorientations axross grain boundaries for specified phase

% Boundary misorientations
misor = all_B.misorientation;

% AXES of the boundary dis/misorientation
ax = (misor.axis);

% ANGLES of the boundary dis/misorientation
ang = misor.angle./degree;

% Select only "low-angle" boundaries (i.e., 2-10 degrees) using 'ang'
low_B = all_B(ang<=10 & ang>=2);

% Misorientations from the low-angle boundaries
lowAngMis = low_B.misorientation;


%% Compute the axes directions relative to the specimen reference frame
% (rather than in the crystal reference frame)

% Pairs of boundary orientations
oris = eLowAngle(low_B.ebsdId).orientations;

% Misorientation axes in specimen coordinates
specLowMisAx = axis(oris(:,1),oris(:,2));



%% Plots:

% plot axis distribtuion of low-angle neighboring misorientations
figure,
plotAxisDistribution(lowAngMis,'antipodal','smooth','figSize','small')
cb = mtexColorbar('Title','M.U.D.');
setColorRange([1 max(cb.Limits)]);
cb.Ticks = linspace(min(cb.Limits),max(cb.Ticks),3);
saveas(gcf,sprintf('%s_%s_lowAngleMisOrs_CRYSTAL_REF.png',sampleName,phase));


% plot axes in specimen reference frame
figure
plot(specLowMisAx,'antipodal','lower','smooth','halfwidth',10*degree,'figSize','small')
drawnow
hold on
cb = mtexColorbar('Title','M.U.D.');
setColorRange([1 max(cb.Limits)]);
cb.Ticks = [min(cb.Limits):.5:max(cb.Ticks)];
annotation('textbox',[0 .8 0 .2],'String',sprintf('n = %i',length(specLowMisAx)),'FitBoxToText','on');
saveas(gcf,sprintf('%s_%s_lowAngleMisOrs_SPEC_REF.png',sampleName,phase));
