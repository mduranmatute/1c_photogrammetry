Base_Loc= '..\Examples\Data_Files\MATLAB\Data_process_steps\Measurements\';
SaveB_Loc='..\Examples\Data_Files\NETCDF\Data_process_steps\Measurements\';

list = ls([Base_Loc '*.mat']);
[nz,~] = size(list);

for iz = 1:nz
%===================================================== Name and create file
    Name_sep = regexp(list(iz,:), '\.', 'split');
    NCF_name = [Name_sep{1} '.nc'];
    ncid = netcdf.create([SaveB_Loc NCF_name],'NC_WRITE');
    
    Step_Name = regexp(Name_sep{1}, '\_', 'split');
    Filename = [Base_Loc, list(iz,:)];
    DataFile = load(Filename); 
    
    if strcmp(Step_Name{1},'DotsInCam')
        CamDPos = DataFile.camdots;
        ImNames = DataFile.image_names;
        LocSep = regexp(ImNames{1}, '\\', 'split');
        LocNm = [];
        for nm = 1:length(LocSep)-1
            LocNm = [LocNm LocSep{nm} '\'];
        end
        LnNm = length(LocNm);
        
        NumF = length(CamDPos);
        LngC = zeros(1,NumF);
        for il = 1:NumF
            LngC(il) = length(CamDPos(il).x);
        end
        MxLP = max(LngC);
        Xpos = zeros(NumF,MxLP)*NaN;
        Ypos = Xpos;
        
        for il = 1:NumF
            for ij = 1:LngC(il)
                Xpos(il,ij) = CamDPos(il).x(ij);
                Ypos(il,ij) = CamDPos(il).y(ij);
            end
        end
%======================================================== Define dimensions
        dimidNPat = netcdf.defDim(ncid,'Num_Patterns',NumF);        
        dimidInd = netcdf.defDim(ncid,'Max_Ind',MxLP);     
%=============================== Define IDs for the the dimension variables
        NPattern_ID = netcdf.defVar(ncid,'N_Patterns','NC_Int',[dimidNPat]);
        netcdf.putAtt(ncid,NPattern_ID,'long_name','Number of patterns');
        netcdf.putAtt(ncid,NPattern_ID,'units','[]');
        netcdf.putAtt(ncid,NPattern_ID,'description',['Number of ' ...
                'projected patterns.']);
        
        MaxDots_ID = netcdf.defVar(ncid,'Max_Dots','NC_int',[dimidInd]);
        netcdf.putAtt(ncid,MaxDots_ID,'long_name',['Maximum number of ' ...
                                        'dots']);
        netcdf.putAtt(ncid,MaxDots_ID,'units','[]');
        netcdf.putAtt(ncid,MaxDots_ID,'description',['Maximum number ' ...
            'of dots detected in the photographed images. The maximum ' ...
            'is used to place all detected dots in a single matrix. ' ...
            'NaN values indicate non-existent dot positions.']);
%============================================== Define ID for main varibles
%---------------------- position of the detected dots in photographs
        DotsX_ID = netcdf.defVar(ncid,'X_pos','double',...
            [dimidNPat dimidInd]);
        netcdf.putAtt(ncid,DotsX_ID,'long_name',['X position of dots ' ...
            'points.']);
        netcdf.putAtt(ncid,DotsX_ID,'units','pixels');
        netcdf.putAtt(ncid,DotsX_ID,'description',['X position of the ' ...
            'dots detected in the photographed images.']);

        DotsY_ID = netcdf.defVar(ncid,'Y_pos','double',...
            [dimidNPat dimidInd]);
        netcdf.putAtt(ncid,DotsY_ID,'long_name',['Y position of dots ' ...
            'points.']);
        netcdf.putAtt(ncid,DotsY_ID,'units','pixels');
        netcdf.putAtt(ncid,DotsY_ID,'description',['Y position of the ' ...
            'dots detected in the photographed images.']);
%---------------------- Global 
        varid = netcdf.getConstant('GLOBAL');
        netcdf.putAtt(ncid,varid,'Description',['Data '...
            'used in the manuscript "High-resolution single-camera ' ...
            'photogrammetry: incorporation of refraction at a fluid ' ...
            'interface" by Gonzalez-Vera, A. S., Wilting, T. J. S., ' ...
            'Holten, A. P. C., van Heijst, G. J. F. & Duran-Matute, M. '...
            'Experiments in Fluids 61(1),3. Data was created and used ' ...
            'by the authors and described in the manuscript. The ' ...
            'image files processed are lcated in:' LocNm]);
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
        netcdf.putVar(ncid,NPattern_ID,1:NumF);
        netcdf.putVar(ncid,MaxDots_ID,1:MxLP);
        
    %-------------------------- Store measurements and interpolated surface
        netcdf.putVar(ncid,DotsX_ID,Xpos);
        netcdf.putVar(ncid,DotsY_ID,Ypos);
    
    %We're done, close the netcdf
        netcdf.close(ncid);
        
    elseif strcmp(Step_Name{1},'Projected')
        PatPos = DataFile.patdots;
        NumF = length(PatPos);
        LngC = zeros(1,NumF);
        for il = 1:NumF
            LngC(il) = length(PatPos(il).x);
        end
        MxLP = max(LngC);
        Xpos = zeros(NumF,MxLP)*NaN;
        Ypos = Xpos;
        
        for il = 1:NumF
            for ij = 1:LngC(il)
                Xpos(il,ij) = PatPos(il).x(ij);
                Ypos(il,ij) = PatPos(il).y(ij);
            end
        end
%======================================================== Define dimensions
        dimidNPat = netcdf.defDim(ncid,'Num_Patterns',NumF);        
        dimidInd = netcdf.defDim(ncid,'Max_Ind',MxLP);     
%=============================== Define IDs for the the dimension variables
        NPattern_ID = netcdf.defVar(ncid,'N_Patterns','NC_Int',[dimidNPat]);
        netcdf.putAtt(ncid,NPattern_ID,'long_name','Number of patterns');
        netcdf.putAtt(ncid,NPattern_ID,'units','[]');
        netcdf.putAtt(ncid,NPattern_ID,'description',['Number of ' ...
                'projected patterns.']);
        
        MaxDots_ID = netcdf.defVar(ncid,'Max_Dots','NC_int',[dimidInd]);
        netcdf.putAtt(ncid,MaxDots_ID,'long_name',['Maximum number of ' ...
                                        'dots']);
        netcdf.putAtt(ncid,MaxDots_ID,'units','[]');
        netcdf.putAtt(ncid,MaxDots_ID,'description',['Maximum number ' ...
            'of dots detected in the projected patterns.']);
%============================================== Define ID for main varibles
%---------------------- position of the detected dots in photographs
        DotsX_ID = netcdf.defVar(ncid,'X_pos','double',...
            [dimidNPat dimidInd]);
        netcdf.putAtt(ncid,DotsX_ID,'long_name',['X position of dots ' ...
            'points.']);
        netcdf.putAtt(ncid,DotsX_ID,'units','pixels');
        netcdf.putAtt(ncid,DotsX_ID,'description',['X position of the ' ...
            'dots detected in the projected patterns.']);

        DotsY_ID = netcdf.defVar(ncid,'Y_pos','double',...
            [dimidNPat dimidInd]);
        netcdf.putAtt(ncid,DotsY_ID,'long_name',['Y position of dots ' ...
            'points.']);
        netcdf.putAtt(ncid,DotsY_ID,'units','pixels');
        netcdf.putAtt(ncid,DotsY_ID,'description',['Y position of the ' ...
            'dots detected in the projected patterns.']);
%---------------------- Global 
        varid = netcdf.getConstant('GLOBAL');
        netcdf.putAtt(ncid,varid,'Description',['Data '...
            'used in the manuscript "High-resolution single-camera ' ...
            'photogrammetry: incorporation of refraction at a fluid ' ...
            'interface" by Gonzalez-Vera, A. S., Wilting, T. J. S., ' ...
            'Holten, A. P. C., van Heijst, G. J. F. & Duran-Matute, M. '...
            'Experiments in Fluids 61(1),3. Data was created and used ' ...
            'by the authors and described in the manuscript.']);
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
        netcdf.putVar(ncid,NPattern_ID,1:NumF);
        netcdf.putVar(ncid,MaxDots_ID,1:MxLP);
        
    %-------------------------- Store measurements and interpolated surface
        netcdf.putVar(ncid,DotsX_ID,Xpos);
        netcdf.putVar(ncid,DotsY_ID,Ypos);

    %We're done, close the netcdf
        netcdf.close(ncid); 
        
    elseif strcmp(Step_Name{1},'SortedDots')
        DotSPos = DataFile.image_dots;
        NumF = length(DotSPos);
        LngC = zeros(1,NumF);
        for il = 1:NumF
            LngC(il) = length(DotSPos(il).x);
        end
        MxLP = max(LngC);
        Xpos = zeros(NumF,MxLP)*NaN;
        Ypos = Xpos;
        
        for il = 1:NumF
            for ij = 1:LngC(il)
                Xpos(il,ij) = DotSPos(il).x(ij);
                Ypos(il,ij) = DotSPos(il).y(ij);
            end
        end
%======================================================== Define dimensions
        dimidNPat = netcdf.defDim(ncid,'Num_Patterns',NumF);        
        dimidInd = netcdf.defDim(ncid,'Max_Ind',MxLP);     
%=============================== Define IDs for the the dimension variables
        NPattern_ID = netcdf.defVar(ncid,'N_Patterns','NC_Int',[dimidNPat]);
        netcdf.putAtt(ncid,NPattern_ID,'long_name','Number of patterns');
        netcdf.putAtt(ncid,NPattern_ID,'units','[]');
        netcdf.putAtt(ncid,NPattern_ID,'description',['Number of ' ...
                'projected patterns.']);
        
        MaxDots_ID = netcdf.defVar(ncid,'Max_Dots','NC_int',[dimidInd]);
        netcdf.putAtt(ncid,MaxDots_ID,'long_name',['Maximum number of ' ...
                                        'sorted dots']);
        netcdf.putAtt(ncid,MaxDots_ID,'units','[]');
        netcdf.putAtt(ncid,MaxDots_ID,'description',['Maximum number ' ...
            'of dots sorted by comparing the projected dot patterns ' ...
            'to the photographed patterns.']);
%============================================== Define ID for main varibles
%---------------------- position of the detected dots in photographs
        DotsX_ID = netcdf.defVar(ncid,'X_pos','double',...
            [dimidNPat dimidInd]);
        netcdf.putAtt(ncid,DotsX_ID,'long_name',['X position of sorted' ...
            'points.']);
        netcdf.putAtt(ncid,DotsX_ID,'units','pixels');
        netcdf.putAtt(ncid,DotsX_ID,'description',['X position of the ' ...
            'sorted dots obtained by comparing projected patterns with' ...
            ' patterns captured by photographs.']);

        DotsY_ID = netcdf.defVar(ncid,'Y_pos','double',...
            [dimidNPat dimidInd]);
        netcdf.putAtt(ncid,DotsY_ID,'long_name',['Y position of the ' ...
            'sorted dots.']);
        netcdf.putAtt(ncid,DotsY_ID,'units','pixels');
        netcdf.putAtt(ncid,DotsY_ID,'description',['Y position of the ' ...
            'sorted dots obtained by comparing projected patterns with' ...
            ' patterns captured by photographs.']);
%---------------------- Global 
        varid = netcdf.getConstant('GLOBAL');
        netcdf.putAtt(ncid,varid,'Description',['Data '...
            'used in the manuscript "High-resolution single-camera ' ...
            'photogrammetry: incorporation of refraction at a fluid ' ...
            'interface" by Gonzalez-Vera, A. S., Wilting, T. J. S., ' ...
            'Holten, A. P. C., van Heijst, G. J. F. & Duran-Matute, M. '...
            'Experiments in Fluids 61(1),3. Data was created and used ' ...
            'by the authors and described in the manuscript.']);
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
        netcdf.putVar(ncid,NPattern_ID,1:NumF);
        netcdf.putVar(ncid,MaxDots_ID,1:MxLP);
        
    %-------------------------- Store measurements and interpolated surface
        netcdf.putVar(ncid,DotsX_ID,Xpos);
        netcdf.putVar(ncid,DotsY_ID,Ypos);

    %We're done, close the netcdf
        netcdf.close(ncid);
    end
end