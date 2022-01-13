Base_Loc = '..\Examples\Data_Files\MATLAB\Patterns\';
SaveB_Loc = '..\Examples\Data_Files\NETCDF\Patterns\';

list = ls([Base_Loc '*.mat']);
[nz,~] = size(list);

for iz = 1:nz
    Filename = [Base_Loc, list(iz,:)];
    PatternFile = load(Filename);

    Pattern = PatternFile.DotsToProject(:,:,:) * 1;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Save to NETCDF files                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%========================================= Value of the variable dimensions
        [Y_dim,X_dim,N_dim] = size(Pattern);
        YL = 1:Y_dim;
        XL = 1:X_dim;
        NL = 1:N_dim;

%===================================================== Name and create file
        %cd(SaveB_Loc)
        Name_sep = regexp(list(iz,:), '\.', 'split');
        NCF_name = [Name_sep{1} '.nc'];
        ncid = netcdf.create([SaveB_Loc NCF_name],'NC_WRITE');
        
%======================================================== Define dimensions
        dimidNPat = netcdf.defDim(ncid,'Num_Pattern',N_dim);        
        dimidInd = netcdf.defDim(ncid,'SingleValue',1);
        dimidX = netcdf.defDim(ncid,'X_p',X_dim);
        dimidY = netcdf.defDim(ncid,'Y_p',Y_dim);
        
%=============================== Define IDs for the the dimension variables
        NPat_ID = netcdf.defVar(ncid,'Num_Pattern','NC_Byte',[dimidNPat]);
        netcdf.putAtt(ncid,NPat_ID,'long_name',...
                                   'Number of projected patterns');
        netcdf.putAtt(ncid,NPat_ID,'units','[]');
        netcdf.putAtt(ncid,NPat_ID,'description',['Number of pattern ' ...
            'images that will be projected.']);
        
        X_ID = netcdf.defVar(ncid,'X_p','NC_Short',[dimidX]);
        netcdf.putAtt(ncid,X_ID,'long_name','X resolution');
        netcdf.putAtt(ncid,X_ID,'units','Pixels');
        netcdf.putAtt(ncid,X_ID,'description',['Horizontal pixel ' ...
            'resolution of the images projected. Tied to the Maximum ' ...
            'resolution of the Projector, not the computer screen.']);
        
        Y_ID = netcdf.defVar(ncid,'Y_p','NC_Short',[dimidY]);
        netcdf.putAtt(ncid,Y_ID,'long_name','Y resolution');
        netcdf.putAtt(ncid,Y_ID,'units','Pixels');
        netcdf.putAtt(ncid,Y_ID,'description',['Vertical pixel ' ...
            'resolution of the images projected. Tied to the Maximum ' ...
            'resolution of the Projector, not the computer screen.']);
        
%=============================================== Define ID for main varible
        % Images of the projected patterns
        PatternImg_ID = netcdf.defVar(ncid,'Pattern_Images','NC_Byte',...
            [dimidY dimidX dimidNPat]);
        netcdf.putAtt(ncid,PatternImg_ID,'long_name',['Images of ' ...
            'patterns that will be projected.']);
        netcdf.putAtt(ncid,PatternImg_ID,'units','pixels');
        netcdf.putAtt(ncid,PatternImg_ID,'description',['Images of ' ...
            'the patterns that will be projected. In the MATLAB ' ...
            'script, this matirx is implemented as a logical data type.']);
%=========================================== Give a description of the file        
        % Global 
        varid = netcdf.getConstant('GLOBAL');
        netcdf.putAtt(ncid,varid,'Description',['Experimental data '...
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
%         netcdf.putAtt(ncid,varid,'Disclaimer',['This version of the ' ...
%             'data is made available for review purposes. Upon ' ...
%             'acceptance of the manuscript, the data will be made ' ...
%             'publicly available through the 4TU.Center for Research ' ...
%             'data at data.4tu.nl']);
        netcdf.putAtt(ncid,varid,'Source','Laboratory Experiments');
        netcdf.putAtt(ncid,varid,'Creation date',datestr(now));
        
        netcdf.putAtt(ncid,varid,'Title',['Data used for the ' ...
            '"High-resolution single-camera photogrammetry: ' ...
            'incorporation of refraction at a fluid interface"']);
        
        % End defining NetCDF
        netcdf.endDef(ncid);
        
        %-------------------------- Store dimension variables
        netcdf.putVar(ncid,NPat_ID,NL);
        netcdf.putVar(ncid,X_ID,XL);
        netcdf.putVar(ncid,Y_ID,YL);
        
        %-------------------------- Store patterns
        netcdf.putVar(ncid,PatternImg_ID,Pattern);                
        %We're done, close the netcdf
        netcdf.close(ncid);
% 
end


Patt = ncread([SaveB_Loc 'CalibrationPattern_Example.nc'],'Pattern_Images');
imgname=['C:\Users\agonzalez\surfdrive\Work\Photogrammetry\Examples\'...
    'Calibration\Photographs\Projector\z_7\Proj_z7_img_1_0016.png'];
Icam = imread(imgname);
Ipat = Patt(:,:,1);
figure()
 subplot(2,1,1); imagesc(Ipat);
        title('pattern dots');
        subplot(2,1,2); imagesc(Icam); 
        colormap(([1:255; 1:255; 1:255])'/256);
        title('image dots');
