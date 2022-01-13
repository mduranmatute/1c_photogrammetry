Base_Loc = '..\Examples\Data_Files\MATLAB\Measurements\';
SaveB_Loc = '..\Examples\Data_Files\NETCDF\Measurements\';

list = ls([Base_Loc '*.mat']);
[nz,~] = size(list);

for iz = 1:nz
    Filename = [Base_Loc, list(iz,:)];
    DataFile = load(Filename);
    
    MesO = DataFile.Pmid;
    InX = DataFile.xq;
    InY = DataFile.yq;
    InZ = DataFile.vq;
    
    [N_dim,N_var] = size(MesO);
    [X_dim, Y_dim] = size(InZ);
    NL = 1:N_dim;
%===================================================== Name and create file
    Name_sep = regexp(list(iz,:), '\.', 'split');
    NCF_name = [Name_sep{1} '.nc'];
    ncid = netcdf.create([SaveB_Loc NCF_name],'NC_WRITE');
        
%======================================================== Define dimensions
    dimidNPnt = netcdf.defDim(ncid,'Num_Points',N_dim);        
    dimidInd = netcdf.defDim(ncid,'x,y,z',3);
    dimidX = netcdf.defDim(ncid,'X_int',X_dim);
    dimidY = netcdf.defDim(ncid,'Y_int',Y_dim);
    
%=============================== Define IDs for the the dimension variables
    NPoint_ID = netcdf.defVar(ncid,'Num_Points','NC_Int',[dimidNPnt]);
    netcdf.putAtt(ncid,NPoint_ID,'long_name',...
                               'Number of measured points');
    netcdf.putAtt(ncid,NPoint_ID,'units','[]');
    netcdf.putAtt(ncid,NPoint_ID,'description',['Number of ' ...
            'measurements obtained from the projected pattern.']);
        
    X_ID = netcdf.defVar(ncid,'X_int','double',[dimidX]);
    netcdf.putAtt(ncid,X_ID,'long_name','X position');
    netcdf.putAtt(ncid,X_ID,'units','mm');
    netcdf.putAtt(ncid,X_ID,'description',['Interpolated horizontal ' ...
            'position of the calculated surface height. The ' ...
            'interpolation is employed to place the measurement in a ' ...
            'uniform grid.']);
        
    Y_ID = netcdf.defVar(ncid,'Y_int','double',[dimidY]);
    netcdf.putAtt(ncid,Y_ID,'long_name','Y position');
    netcdf.putAtt(ncid,Y_ID,'units','mm');
    netcdf.putAtt(ncid,Y_ID,'description',['Interpolated vertical ' ...
            'position of the calculated surface height. The ' ...
            'interpolation is employed to place the measurement in a ' ...
            'uniform grid.']);
%============================================== Define ID for main varibles
%---------------------- Measurements obtained from the projected patterns
    Measurement_ID = netcdf.defVar(ncid,'Measurements','double',...
            [dimidNPnt dimidInd]);
    netcdf.putAtt(ncid,Measurement_ID,'long_name',['Measurement ' ...
            'points.']);
    netcdf.putAtt(ncid,Measurement_ID,'units','mm');
    netcdf.putAtt(ncid,Measurement_ID,'description',['Measurement ' ...
            'points obtained from the projected patterns. Each point ' ...
            'is given by its 3D spatial coordinate (x,y,z) such that: ' ...
            'x:(1:Num_Points,1), y:(1:Num_Points,2), z:(1:Num_Points,3)']);

%---------------------- Interpolated surface
    Surface_ID = netcdf.defVar(ncid,'Surface','double',...
            [dimidX dimidY]);
    netcdf.putAtt(ncid,Surface_ID,'long_name','Interpolated surface');
    netcdf.putAtt(ncid,Surface_ID,'units','mm');
    netcdf.putAtt(ncid,Surface_ID,'description',...
            'Surface obtained from the interpolation of the measurements');
        
%---------------------- Global 
    varid = netcdf.getConstant('GLOBAL');
    netcdf.putAtt(ncid,varid,'Description',['Data '...
            'used in the manuscript "High-resolution single-camera ' ...
            'photogrammetry: incorporation of refraction at a fluid ' ...
            'interface" by Gonzalez-Vera, A. S., Wilting, T. J. S., ' ...
            'Holten, A. P. C., van Heijst, G. J. F. & Duran-Matute, M. '...
            'Experiments in Fluids 61(1),3. Data was created and used ' ...
            'for experiments by the authors and described in the ' ...
            'manuscript.']);
    netcdf.putAtt(ncid,varid,'Authors',['Gonzalez-Vera, A. S., ' ...
            'Wilting, T. J. S., Holten, A. P. C., van Heijst, ' ...
            'G. J. F. & Duran-Matute, M.']);        
    netcdf.putAtt(ncid,varid,'email for correspondence',...
            'm.duran.matute@tue.nl');
    netcdf.putAtt(ncid,varid,'Institution',['Eindhoven University ' ...
            'of Technology, Department of Applied Physics']);
    netcdf.putAtt(ncid,varid,'Source','Laboratory Experiments');
    netcdf.putAtt(ncid,varid,'Creation date',datestr(now));
        
    netcdf.putAtt(ncid,varid,'Title',['Data used for the ' ...
            '"High-resolution single-camera photogrammetry: ' ...
            'incorporation of refraction at a fluid interface"']);
        
    % End defining NetCDF
    netcdf.endDef(ncid);
        
    %-------------------------- Store dimension variables
    netcdf.putVar(ncid,NPoint_ID,NL);
    netcdf.putVar(ncid,X_ID,InX(1,:));
    netcdf.putVar(ncid,Y_ID,InY(:,1));
        
    %-------------------------- Store measurements and interpolated surface
    netcdf.putVar(ncid,Measurement_ID ,MesO);
    netcdf.putVar(ncid,Surface_ID,InZ);
    
    %We're done, close the netcdf
    netcdf.close(ncid);
end