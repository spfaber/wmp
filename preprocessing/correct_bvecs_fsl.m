function ecc_bvecs(ecclog,bvecsfile)
% ecc_bvecs(ecclog,bvecsfile)
%
% S.Jbabdi 03/2009

% read bvecs file
bvecs = load(bvecsfile);
if(size(bvecs,2)==3 && size(bvecs,1)>3)
    bvecs = bvecs';
end

% read ecc-log file
unix(['cat ' ecclog ' | grep Final -A 4 | sed s/"Final result:"//g | sed s/--//g > /tmp/tmpecclog']);

mat = load('/tmp/tmpecclog');

% rotate bvecs
rotbvecs = zeros(size(bvecs));
for i = 1:size(bvecs,2)
    M = mat((i-1)*4+1:i*4,:);
    M = M(1:3,1:3);
    
    % extract rotation matrix
    [u,s,v] = svd(M*M');
    R = inv(u*sqrt(s)*v')*M;
    
    rotbvecs(:,i) = R*bvecs(:,i);
end

save([bvecsfile '_rot'],'rotbvecs','-ascii');


