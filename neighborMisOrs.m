function [misXTAL,misSPEC] = neighborMisOrs(ebsd)
%% neighbor misorientations per orientation

% compute all adjacent measurements
[~,~,I_FD] = spatialDec([ebsd.prop.x(:), ebsd.prop.y(:)],ebsd.unitCell,'unitCell');
A_D = I_FD.' * I_FD;

% order
n = 2;

A_D1 = A_D;
for i = 1:n-1  
  A_D = A_D + A_D*A_D1 + A_D1*A_D;
end
clear A_D1

% extract adjacent pairs
[Dl, Dr] = find(A_D);

% take only ordered pairs of same, indexed phase 
use = Dl > Dr & ebsd.phaseId(Dl) == ebsd.phaseId(Dr) & ebsd.isIndexed(Dl);
Dl = Dl(use); Dr = Dr(use);
phaseId = ebsd.phaseId(Dl);

% calculate misorientation angles and axes
misXTAL = struct();
misSPEC = struct();
axPhase = struct();
omegaPhase = struct();

% iterate all phases
for p=1:numel(ebsd.phaseMap)
  
  currentPhase = phaseId == p;
  if any(currentPhase)
    
    o_Dl = orientation(ebsd.rotations(Dl(currentPhase)),ebsd.CSList{p});
    o_Dr = orientation(ebsd.rotations(Dr(currentPhase)),ebsd.CSList{p});
    phase = ebsd(Dl(currentPhase)).mineral;
    
    % misorientations in crystal coordinates
    misXTAL.(phase) = o_Dl.\o_Dr;
    
    % misorientations in specimen coordinates
    misSPEC.(phase) = orientation('axis',axis(o_Dl,o_Dr),'angle',angle(o_Dl,o_Dr));
    
  end
end


end

