function [] = createSfunctions(qp_W,qp_gradJ,qp_gradhT,qp_h)

% Get correct root
% temp.parts = regexp(cd,'\','split')';
% temp.idx = find(strcmp(temp.parts,'GroupWork'));
% root = fullfile(temp.parts{1:temp.idx-1});
folder = 'SfunctionGeneration';
pathcd = which('createSfunctions');
temp.idx = regexp(pathcd,folder);
rootPath = pathcd(1:temp.idx+numel(folder));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Script Settings                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sfunc_qp        = 1;       % Create s-Functions for derivative information
ndef_real_t     = 1;       % Define flag to define real_t = double
mex_getFunc     = 1;       % Compile mex-functions to read out work array sizes
%-------------------------------------------------------------------------%


% make sure that there are the necessary functions in the workspace
% % funcNames = {'qp_W','qp_gradJ','qp_gradhT','qp_h'};
% % 
% % for ii=1:length(funcNames)
% %     if ~(1==exist(funcNames{ii}))
% %         error('Error. \n%s not found in the workspace.',funcNames{ii})
% %     end
% % end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Build s-Functions                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
buildFolder = 'cSourceFiles';
if sfunc_qp
    
    % % Create building folder for sfun generation
      % remove old folder and generate new one
        if isfolder(buildFolder)
            try
                cd(buildFolder);
                delete('get*');
                delete('qp*');
                delete('Wrap*');
                cd ..
                rmdir(buildFolder);
            catch
                warning('createSfunctions: Couldn,t remove previously compiled functions.');
            end
        end
        mkdir(buildFolder);
        cd(buildFolder);
      % Copy needed c-files to current folder
        copyfile(fullfile(rootPath,'getArrSizes.c'),'.\');
        copyfile(fullfile(rootPath,'getSparsity.c'),'.\');
        copyfile(fullfile(rootPath,'Wrap_qp_W_orig.c'),'.\');
        copyfile(fullfile(rootPath,'Wrap_qp_gradJ_orig.c'),'.\');
        copyfile(fullfile(rootPath,'Wrap_qp_gradhT_orig.c'),'.\');
        copyfile(fullfile(rootPath,'Wrap_qp_h_orig.c'),'.\');
      % Copy needed header files to current folder
        copyfile(fullfile(rootPath,'Wrap_qp_W.h'),'.\');
        copyfile(fullfile(rootPath,'Wrap_qp_gradJ.h'),'.\');
        copyfile(fullfile(rootPath,'Wrap_qp_gradhT.h'),'.\');
        copyfile(fullfile(rootPath,'Wrap_qp_h.h'),'.\');
        
    % Code generation of CasADi functions
      qp_W.generate('qp_W',struct('with_header',true, 'with_export', false));
      qp_gradJ.generate('qp_gradJ',struct('with_header',true, 'with_export', false));
      qp_gradhT.generate('qp_gradhT',struct('with_header',true, 'with_export', false));
      qp_h.generate('qp_h',struct('with_header',true, 'with_export', false));
    %
    % Update mex
    if mex_getFunc
        mex getArrSizes.c qp_W.c qp_gradJ.c qp_gradhT.c qp_h.c   
        mex getSparsity.c qp_W.c qp_gradJ.c qp_gradhT.c qp_h.c 
    end
    
    % Get size of working arrays (substitute zero array size by one)
    arrSizes = getArrSizes(1, 0, 1, zeros(4,4));
    ind = arrSizes == 0;
    arrSizes(ind) = 1;
    
    % Get sparsity information of resulting matrices
    arrSparsity = getSparsity(1, 0, 1, zeros(3,10));

    % Write arraySizes and sparsity information to each wrapper
    wrapNames       = {'Wrap_qp_W.c', 'Wrap_qp_gradJ.c', ...
                       'Wrap_qp_gradhT.c', 'Wrap_qp_h.c'};
    wrapNamesOrig   = {'Wrap_qp_W_orig.c', 'Wrap_qp_gradJ_orig.c', ...
                       'Wrap_qp_gradhT_orig.c', 'Wrap_qp_h_orig.c'};

    buffSizes = [numel_out(qp_W);sqrt(numel_out(qp_W));...
                 numel_out(qp_h)*sqrt(numel_out(qp_W));numel_out(qp_h)];
    % Parse wrapper functions and replace with actual working array sizes
    for ii = 1:4
        A{1} = ['const double* argument[' num2str(arrSizes(1,ii)) ']; \n'];
        A{2} = ['double buff[' num2str(buffSizes(ii)) ']; \n'];
        A{3} = ['double* result[' num2str(arrSizes(2,ii)) ']; \n'];
        A{4} = ['int iw[' num2str(arrSizes(3,ii)) ']; \n'];
        A{5} = ['double w[' num2str(arrSizes(4,ii)) ']; \n'];
        A{6} = ['int n_row = ' num2str(arrSparsity(1,ii)) '; \n'];
        A{7} = ['int n_col = ' num2str(arrSparsity(2,ii)) '; \n'];
        A{8} = ['int nnz = ' num2str(arrSparsity(3,ii)) '; \n'];
        A{9} = ['int arr_off[' num2str(arrSparsity(2,ii)) ']; \n'];
        A{10} = ['int arr_nnz[' num2str(arrSparsity(3,ii)) ']; \n'];

        fid = fopen(wrapNamesOrig{ii},'r');
        jj = 1;
        kk = 1;
        tline = fgetl(fid);
        B{jj} = [tline '\n'];
        while ischar(tline)
            jj = jj+1;
            tline = fgetl(fid);

            B{jj} = [tline '\n'];

            if jj > 7 && jj <18
                B{jj} = A{kk};
                kk = kk + 1;
            end
        end
        fclose(fid);

        fid = fopen(wrapNames{ii},'wt');
        for jj = 1:length(B)-1
            fprintf(fid, B{jj});
        end
        fclose(fid);
        clear A B;
    end
    
    % change folder and delete old mex and c files
    if ~(exist('..\sFunctions','dir')==7)
        mkdir('../sFunctions/');
    end % if
    cd ../sFunctions/
    delete('sfun*');
    
    def                 = legacy_code('initialize');
    def.IncPaths        = {['../',buildFolder,'/']};
    def.SrcPaths        = {['../',buildFolder,'/']};
    
    % Hessian of Lagrangian
    def.SourceFiles     = {'qp_W.c', 'Wrap_qp_W.c'};
    def.HeaderFiles     = {'qp_W.h', 'Wrap_qp_W.h'};
    def.SFunctionName   = 'sfun_qp_W';
    def.OutputFcnSpec   = ['void Wrap_qp_W(double u1[], double u2[], double y1[' num2str(buffSizes(1)) '])'];
    legacy_code('sfcn_cmex_generate', def);
    % MSVC 2017 does not know real_t
    if ndef_real_t 
        legacy_code('compile', def, '-Dreal_t=double' );
    else
        legacy_code('compile', def);        
    end
%     legacy_code('slblock_generate', def);
    legacy_code('sfcn_tlc_generate', def);

    % Gradient of objective function
    def.SourceFiles     = {'qp_gradJ.c', 'Wrap_qp_gradJ.c'};
    def.HeaderFiles     = {'qp_gradJ.h', 'Wrap_qp_gradJ.h'};
    def.SFunctionName   = 'sfun_qp_gradJ';
    def.OutputFcnSpec   = ['void Wrap_qp_gradJ(double u1[], double u2[], double y1[' num2str(buffSizes(2)) '])'];
    legacy_code('sfcn_cmex_generate', def);
    % MSVC 2017 does not know real_t
    if ndef_real_t 
        legacy_code('compile', def, '-Dreal_t=double');
    else
        legacy_code('compile', def);        
    end
%     legacy_code('slblock_generate', def);
    legacy_code('sfcn_tlc_generate', def);
    
    % Gradient of boundary conditions
    def.SourceFiles     = {'qp_gradhT.c', 'Wrap_qp_gradhT.c'};
    def.HeaderFiles     = {'qp_gradhT.h', 'Wrap_qp_gradhT.h'};
    def.SFunctionName   = 'sfun_qp_gradhT';
    def.OutputFcnSpec   = ['void Wrap_qp_gradhT(double u1[], double u2[], double y1[' num2str(buffSizes(3)) '])'];
    legacy_code('sfcn_cmex_generate', def);
    % MSVC 2017 does not know real_t
    if ndef_real_t 
        legacy_code('compile', def, '-Dreal_t=double');
    else
        legacy_code('compile', def);        
    end
%     legacy_code('slblock_generate', def);
    legacy_code('sfcn_tlc_generate', def);
    
    % Boundary conditions
    def.SourceFiles     = {'qp_h.c', 'Wrap_qp_h.c'};
    def.HeaderFiles     = {'qp_h.h', 'Wrap_qp_h.h'};
    def.SFunctionName   = 'sfun_qp_h';
    def.OutputFcnSpec   = ['void Wrap_qp_h(double u1[], double u2[], double y1[' num2str(buffSizes(4)) '])'];
    legacy_code('sfcn_cmex_generate', def);
    % MSVC 2017 does not know real_t
    if ndef_real_t 
        legacy_code('compile', def, '-Dreal_t=double');
    else
        legacy_code('compile', def);        
    end
%     legacy_code('slblock_generate', def);
    legacy_code('sfcn_tlc_generate', def);
    
    cd ..
% %     % clear compiling folder and set path
% %     cd(buildFolder);
% %     delete('get*');
% %     delete('qp*');
% %     delete('Wrap*');
% %     cd ..
% %     rmdir(buildFolder);
end
end
