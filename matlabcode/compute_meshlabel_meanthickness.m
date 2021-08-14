function compute_meshlabel_meanthickness(MeshVTK, id, side, OUT_FILENAME1, OUT_FILENAME2)

%addpath('~longxie/SideProjects/Cortical_Surface_Label_Fusion/matlab_code/Basic Functions/matlab/')

%% generate label
% read VTK mesh
m = vtk_polydata_read(MeshVTK);

% find the PROB and thickness field
found_PROB = 0;
found_Thickness = 0;
for ii = 1:length(m.point_data)
   
    if strcmp(m.point_data(ii).name, 'PROB')
       
        mat = m.point_data(ii).data;
        found_PROB = 1;
        
    elseif strcmp(m.point_data(ii).name, 'Thickness')
        
        Thickness = m.point_data(ii).data;
        found_Thickness = 1;
        
    else
        continue;
    end
    
end

if found_PROB == 0 || found_Thickness == 0
    error('field PROB or thickness does not exist. \n');
else
    fprintf('Done loading file. \n')
end

% generate label
[~, idx] = sort(mat, 2);
num_label = size(mat, 2);
label = idx(:, end);
m = vtk_add_point_data(m, 'Label', label, 1);
vtk_polydata_write(MeshVTK, m);

%% generate mean thickness for each label
mean_thickness = zeros(1, num_label);
median_thickness = zeros(1, num_label);
for ii = 1:num_label
    mean_thickness(1,ii) = 2 * mean(Thickness(label == ii));
    median_thickness(1,ii) = 2 * median(Thickness(label == ii));
end

fprintf('Done computing mean thickness. \n')

%% save the result to txt file
if exist(OUT_FILENAME1, 'file')
    delete(OUT_FILENAME1)
end
if exist(OUT_FILENAME2, 'file')
    delete(OUT_FILENAME2)
end
tempmeantext = [id, ',', side];
tempmediantext = [id, ',', side];
for jj = 1:num_label
    tempmeantext = sprintf('%s,%1.10f', tempmeantext, mean_thickness(1, jj));
    tempmediantext = sprintf('%s,%1.10f', tempmediantext, median_thickness(1, jj));
end
system(sprintf('echo "%s" >> %s', tempmeantext, OUT_FILENAME1));
system(sprintf('echo "%s" >> %s', tempmediantext, OUT_FILENAME2));

fprintf('Done! \n')
