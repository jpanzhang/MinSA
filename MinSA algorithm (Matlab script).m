%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%！！！！ Developed by Jipan Zhang, jpanzhang@live.com ！！！！%%%%%%%%%%%
%%%%%%%%！！！Zhao's lab, Southwest univeristy, Chongqing, China ！！！！%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%               Part 1            %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%！！！！     Parameter configuration and image inputting    %%%%%%%%%%%%
clear all;  % Clear workspace
close all;
% Image file must be named as "Sample1_SHF2304PHF159.png" style
     % Sample1 is the your sample name
     % SHF2304 means the number of secondary hair follicles (SHF) is 2304
     % PHF2304 means the number of primary hair follicles (PHF) is 159
% Set your working directory, and put all clean images in it.     
fileFolder=fullfile('C:\Users\Admin\Desktop\MinSA analysis');  
dirOutput=dir(fullfile(fileFolder,'*.png'));
fileNames={dirOutput.name};
% Parameter configuration
Area_ori_img=50;     % Original image area (mm2)
fieldN=10;
fieldS=14;
pixSHF=528;
pixPHF=1372;
R1=0; G1=0; B1=255;      % The RGB value of blue dot, 0-0-255
R2=0; G2=255; B2=0;      % The RGB value of blue dot, 0-255-0
RepeatTimes=10;          % Simulation times for every fieldN〜fieldS
mkdir([fileFolder '\' 'THFD_final_output\'])


for kk=1:length(fileNames)                           % Loop level 1
    imName=fileNames(kk);                            % Image file name
    I=imread([fileFolder '\' cell2mat(imName)]);     % Load image
    I1=I;                                          
    I2=imrotate(I, 90);                             
    cf=[0.8 .^(linspace(2, fieldS*2, fieldS))]';     % The reduction ratio of fieldS
    display(['Image processing on sample' num2str(kk) ' (' num2str(kk) '/' num2str(length(fileNames)) ')']);
    display(['This procedure is very time-consuming, bacause '  num2str(fieldN*fieldS*RepeatTimes) ' sub-images will be extracted and calculated for every sample.']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%               Part 2            %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%   Image processing and MRE calculating    %%%%%%%%%%%%%%%
    for uu=1:2                                      % 1 means horizontal cut; 2 means longitudinal cut
        Img = eval([sprintf('I%d', uu)]);
        [mm,nn,~]=size(Img);                        % Image width and height
        % Get the truth value based on the file sample name
        a1=char(fileNames(kk));
        real_s = str2num(a1(strfind(char(fileNames(kk)), 'SHF')+3 : strfind(char(fileNames(kk)), 'PHF')-1));
        real_p = str2num(a1(strfind(char(fileNames(kk)), 'PHF')+3 : strfind(char(fileNames(kk)), '.png')-1));
        real_t = real_p + real_s;
        sampleName = a1(strfind(char(fileNames(kk)), 'Cropping_')+9 :strfind(char(fileNames(kk)), 'Cropping_')+12);
        % Build empty set
        MRE_dataset=[]; 
        
        % Load progress bar
        show = waitbar(0,'MinSA is running, please wait......'); 
        for  ii=1:fieldN                                 % Loop level 2
               for Num_techN=1:ii
                sub_Img=Img(:, (round(nn/ii)*(Num_techN-1)+1):min((round(nn/ii)*Num_techN), nn), :);  
                assignin('base',[sprintf('Img%d', Num_techN)],sub_Img);                              
            end
            for jj=1:fieldS                              % Loop level 3
               % calculating the  width and height of sub-image
                h = mm * 0.8 .^(linspace(1,fieldS,fieldS))';
                w = (nn * 0.8 .^(linspace(1,fieldS,fieldS))');
                hh= mm * sqrt(h(jj)*w(jj)/ (mm*nn/ii));    
                ww= nn/ii * sqrt(h(jj)*w(jj)/ (mm*nn/ii));
                if cf(jj)*ii >= 1                           
                    mre_nr_S_NiSj= NaN;
                    mre_nr_P_NiSj= NaN;
                    mre_nr_T_NiSj= NaN;
                elseif cf(jj)*ii < 1                        
                    S_NiSj=[];
                    P_NiSj=[];
                    for pp=1:RepeatTimes/2                        % Loop level 4
                        Ncolor1=[];
                        Ncolor2=[];
                        for Num_techN2=1:ii                      
                            sub_Img=eval([sprintf('Img%d', Num_techN2)]);
                            size_sub_Img = size(sub_Img);
                            f=rand(1,2);                     
                            rect=[ceil((size_sub_Img(2)-ww)*f(2)) ceil((size_sub_Img(1)-hh)*f(1)) ceil(ww)  ceil(hh)];  
                            gds_Img=imcrop(sub_Img, rect);
                            [row1, col1] = ind2sub(size(gds_Img), find(gds_Img(:,:,1)==R1 & gds_Img(:,:,2)==G1 & gds_Img(:,:,3)==B1));     
                            [row2, col2] = ind2sub(size(gds_Img), find(gds_Img(:,:,1)==R2 & gds_Img(:,:,2)==G2 & gds_Img(:,:,3)==B2));      
                            Ncolor1 = [Ncolor1, length(row1)];                
                            Ncolor2 = [Ncolor2, length(row2)];               
                            % imwrite(gds_Img,[fileFolder '\outputing_Img\Figure_' 'N' num2str(ii) '_S' num2str(jj) '_order' num2str(pp)...
                            %    '_SHF' num2str(length(row1)/pixSHF/cf(jj))  'PHF' num2str(length(row2)/pixPHF/cf(jj)) '.png']);       % Save sub-images (optional)
                        end
                        S_NiSj(pp)=mean(Ncolor1);               % SHF
                        P_NiSj(pp)=mean(Ncolor2);               % PHF
                    end
                    % Show progress bar
                    waitbar((ii+(uu-1)*fieldN)/(fieldN*2), show)     
                    
                    cf_s=mean((S_NiSj/pixSHF)/cf(jj))/real_s;    
                    cf_p=mean((P_NiSj/pixPHF)/cf(jj))/real_p;      
                    cf_t=mean(((S_NiSj/pixSHF)+(P_NiSj/pixPHF))/cf(jj))/(real_t);    
                    % nr_S_NiSj and nr_P_NiSj
                    nr_S_NiSj=(S_NiSj/pixSHF)/cf(jj)/cf_s;
                    nr_P_NiSj=(P_NiSj/pixPHF)/cf(jj)/cf_p;
                    nr_T_NiSj=(((S_NiSj/pixSHF)+(P_NiSj/pixPHF))/cf(jj))/cf_t;
                    % MRE, calculation of the mean of 600 RE values (MRE) at every fieldN〜fieldS
                    mre_nr_T_NiSj=mean(abs((nr_T_NiSj-real_t)/real_t));  
                end
                % MRE matrix generating
                MRE_dataset=[MRE_dataset,mre_nr_T_NiSj];  
            end
        end
        if uu==1
            MRE_dataset_1=MRE_dataset;
        elseif uu==2
            MRE_dataset_2=MRE_dataset;
        end
    delete(show);
    end
    orig_dataset = (MRE_dataset_1 + MRE_dataset_2)/2;
    assignin('base', [sprintf('orig_dataset%d', kk)], orig_dataset)
   
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%                        Part 3                       %%%%%%%%%%
    %%%%%%%  Curve fitting of MRE dataset (for every sample)    %%%%%%%%%% 
        dataset = reshape(orig_dataset, fieldS, fieldN);
        all_SA5 =[];
        all_SA10 =[];
        for ff=1:fieldN
            x = Area_ori_img *cf;                          % Area series of different sub-images
            y = dataset(:, ff);
            xx=log(x); yy=log(y);                          % Data is linearized
            xx = [ones(length(xx),1), xx ];                
            [coeff,bint,r,rint,stats] = regress(yy,xx);    % Linear regression
            coeff1 = exp(1)^(coeff(1));
            coeff2 = coeff(2);
            SA5 = (0.05/coeff1)^(1/coeff2);     % 5% MRE
            SA10 = (0.10/coeff1)^(1/coeff2);    % 10% MRE
            all_SA5 =[all_SA5, SA5];
            all_SA10 =[all_SA10, SA10];
        end
        dataset=[dataset; all_SA5; all_SA10];
        dataset=[[x; NaN; NaN] dataset];
        dataset=array2table(dataset, 'VariableNames', {'Area'; 'N1'; 'N2'; 'N3';	'N4';	'N5';	'N6';	'N7';	'N8';	'N9';  'N10'}, ...
            'RowNames',{'S1'; 'S2';	'S3';	'S4';	'S5';	'S6';	'S7';	'S8';	'S9';	'S10'; 'S11'; 'S12'; 'S13'; 'S14'; 'fieldS_MRE5%'; 'fieldS_MRE10%'});
        writetable(dataset, [fileFolder '\THFD_final_output\' 'MRE matrix_' char(sampleName) '.csv'], 'WriteRowNames', true);
end       
        

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%                        Part 4                       %%%%%%%%%%
    %%%%%%%    Curve fitting of MRE dataset (for all samples)    %%%%%%%%%% 
    all_dataset = [];
    for kk=1:length(fileNames)  
       orig_dataset =  evalin('base', [sprintf('orig_dataset%d', kk)]);
       all_dataset = [all_dataset; orig_dataset];
    end
        all_dataset = mean(all_dataset, 1);
        all_dataset = reshape(all_dataset, fieldS, fieldN);   
        all_parm=[]
for ff=1:fieldN
    x = Area_ori_img *cf * [0.8 .^(linspace(2, fieldS*2, fieldS))]';       % Area series of different sub-images
    y = all_dataset(:, ff);
    xx=log(x); yy=log(y);                          % Data is linearized
    xx = [ones(length(xx),1), xx ];                
    [coeff,bint,r,rint,stats] = regress(yy,xx);    % Linear regression
    R2=stats(1);
    pv=stats(3);
    coeff1 = exp(1)^(coeff(1));
    coeff2 = coeff(2);
    SA5 = (0.05/coeff1)^(1/coeff2);     % 5% MRE
    SA10 = (0.10/coeff1)^(1/coeff2);    % 10% MRE
    SA5_N = SA5*ff;     % 5% MRE
    SA10_N = SA10*ff;    % 10% MRE
    % parameters pooled
    single = [R2; coeff1; coeff2; pv; SA5; SA10; SA5_N; SA10_N];
    all_parm(:,ff)=single;
end
dataset2 =[all_dataset; all_parm];
dataset2 =[[x; NaN; NaN; NaN; NaN; NaN; NaN; NaN; NaN] dataset2];
dataset2 =array2table(dataset2, 'VariableNames', {'SA'; 'N1'; 'N2';	'N3';	'N4';	'N5';	'N6';	'N7';	'N8';	'N9';  'N10'}, ...
    'RowNames',{'S1'; 'S2';	'S3';	'S4';	'S5';	'S6';	'S7';	'S8';	'S9';	'S10'; 'S11'; 'S12'; 'S13'; 'S14';...
    'R2'; 'coeff1(a)'; 'coeff2(b)'; 'P_value'; 'fieldS_MRE5%'; 'fieldS_MRE10%'; 'fieldS*fieldN_MRE5%'; 'fieldS*fieldN_MRE10%'});
writetable(dataset2, [fileFolder '\THFD_final_output\' 'MRE matrix (for all samples)' '.csv'], 'WriteRowNames', true);
dataset2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%                        End                          %%%%%%%%%%%%%%%%
% For further details, see our research paper:
% Effect of number and area of view-fields on the measurement accuracy of hair follicle density in goats (Capra hircus),Small Ruminant Research,2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
