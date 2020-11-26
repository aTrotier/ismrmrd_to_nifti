function [img,head] = ismrmrd_read_img(struct)
%This function read ismrmrd out.h5 image
%   input :
%           struct.path : string

if nargin < 1 || ~isfield(struct,'path')
    [file,path]=uigetfile('*.h5');
    struct.path = fullfile(path,file);
%error('Missing struct.path which points to the .h5 file');
end

S=hdf5info(struct.path);

if exist(struct.path, 'file')
    dset = ismrmrd.Dataset(struct.path, 'dataset');
else
    error(['File ' struct.path ' does not exist.  Please generate it.'])
end

%% old way
disp(dset.fid.identifier)

S=hdf5info(struct.path);

for i=1:length(S.GroupHierarchy(1).Groups(1).Groups)

attributes=S.GroupHierarchy(1).Groups(1).Groups(i).Datasets(1).Name;
dataset=S.GroupHierarchy(1).Groups(1).Groups(i).Datasets(2).Name;
header=S.GroupHierarchy(1).Groups(1).Groups(i).Datasets(3).Name;

head{i} = h5read(struct.path,header);
img{i}=squeeze(h5read(struct.path,dataset));
end

%% new way to read ?
% S=h5info(struct.path);
% 
% path=S.Groups(1).Groups(1).Name;
% dataset=S.Groups(1).Groups(1).Datasets(2).Name;
% header=S.Groups(1).Groups(1).Datasets(3).Name;
% 
% img=h5read(struct.path,[path '/' dataset]);

%%

end

