% https://github.com/SarahMorgan/Morphometric_Similarity_SZ/blob/master/Gene_analyses.md
% This code was written by Dr Petra Vértes and is taken from Whitaker and Vértes, PNAS 2016, please cite that paper if you this code in your own work.
clear all;
close all;

addpath(genpath('spatial_permutation_test'));
% 

gene_matrix = readmatrix("abagen_15631_gene_expression.csv");
genes=importdata('15631_genes.csv'); 
genes=convertCharsToStrings(genes);
geneindex=1:length(genes);

coor = importdata('308_regions_coordinates.txt');
coor_lh = coor(1:152,:);

for id=1:3
    t_vals = readmatrix(strcat('t_values_asd', num2str(id), '_control.csv'));

    t_vals = t_vals(1:152,1);
    
    X = gene_matrix;
    Y = t_vals;
    
    % z-score:
    X=zscore(X); % by column
    Y=zscore(Y);
%     [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,10,'CV',10);
%     plot(1:10,cumsum(100*PLSPctVar(2,:)),'-bo');
%     xlabel('Number of PLS components');
%     ylabel('Percent Variance Explained in Y');

    %perform full PLS and plot variance in Y explained by top 15 components
    %typically top 2 or 3 components will explain a large part of the variance
    %(hopefully!)
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y);
    dim=15;
    plot(1:dim,cumsum(100*PCTVAR(2,1:dim)),'-o','LineWidth',1.5,'Color',[140/255,0,0]);
    set(gca,'Fontsize',14)
    xlabel('Number of PLS components','FontSize',14);
    ylabel('Percent Variance Explained in Y','FontSize',14);
    grid on

%     dim=6;
%     [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim); 

    %%% plot correlation of PLS component 1 with t-statistic:
    figure
    plot(XS(:,1),Y,'r.')

%     tbl = table(XS(:,1),Y);
%     mdl = fitlm(tbl);
%     plot(mdl)

    [R_1,p_1]=corrcoef(XS(:,1),Y);
    
    xRange=xlim;
    yRange=ylim;
    % spatial permutation test
    perm_id_1 = rotate_parcellation_fixed_lh(coor_lh, 10000);
    p_perm = perm_sphere_p(X(:,1),Y,perm_id_1);

    xlabel('PLS1 scores','FontSize',14);
    ylabel('subtypes-TD t values','FontSize',14);
    str = ['r= ', num2str(round(R_1(1,2),3)),' p= ', num2str(round(p_1(1,2),3)), ', p_perm= ', num2str(round(p_perm,3))];
    text(0.9*xRange(1),0.9*yRange(2),str);
    grid on
    
    % permutation testing to assess significance of PLS result as a function of
    % the number of components (dim) included:

    rep=1000;
    for dim=1:3
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
    temp=cumsum(100*PCTVAR(2,1:dim)); % accumulated percentage
    Rsquared = temp(dim);
        for j=1:rep
            %j
            order=randperm(size(Y,1));
            Yp=Y(order,:);
            

            [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Yp,dim);

            temp=cumsum(100*PCTVAR(2,1:dim));
            Rsq(j) = temp(dim);
        end
    dim;
    R(dim)=Rsquared;
    p(dim)=length(find(Rsq>=Rsquared))/rep;
    end
    figure
    plot(1:dim, p,'ok','MarkerSize',8,'MarkerFaceColor','r');
    xlabel('Number of PLS components','FontSize',14);
    ylabel('p-value','FontSize',14);
    grid on
    
    R
    p

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Bootstrap to get the gene list
    
    %number of bootstrap iterations:
    bootnum=1000;
    
    % Do PLS in 3 dimensions (with 3 components):
    dim=3;
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
    
    %store regions' IDs and weights in descending order of weight for both components:

    [R1,p1]=corr([XS(:,1),XS(:,2), XS(:,3)],Y);

    %align PLS components with desired direction for interpretability 
    if R1(1,1)<0  %this is specific to the data shape we were using - will need ammending
        stats.W(:,1)=-1*stats.W(:,1);
        XS(:,1)=-1*XS(:,1);
    end
    if R1(2,1)<0 %this is specific to the data shape we were using - will need ammending
        stats.W(:,2)=-1*stats.W(:,2);
        XS(:,2)=-1*XS(:,2);
    end
    if R1(3,1)<0 %this is specific to the data shape we were using - will need ammending
        stats.W(:,3)=-1*stats.W(:,3);
        XS(:,3)=-1*XS(:,3);
    end

    [PLS1w,x1] = sort(stats.W(:,1),'descend');
    PLS1ids=genes(x1);
    geneindex1=geneindex(x1);
    [PLS2w,x2] = sort(stats.W(:,2),'descend');
    PLS2ids=genes(x2);
    geneindex2=geneindex(x2);
    [PLS3w,x3] = sort(stats.W(:,3),'descend');
    PLS3ids=genes(x3);
    geneindex3=geneindex(x3);

    %print out results
    csvwrite(strcat('PLS genes results/subtype_', string(id), '_PLS1_ROIscores.csv'),XS(:,1));
    csvwrite(strcat('PLS genes results/subtype_', string(id), '_PLS2_ROIscores.csv'),XS(:,2));
    csvwrite(strcat('PLS genes results/subtype_', string(id), '_PLS3_ROIscores.csv'),XS(:,3));

    %define variables for storing the (ordered) weights from all bootstrap runs
    PLS1weights=[];
    PLS2weights=[];
    PLS3weights=[];
    
    %start bootstrap
    for i=1:bootnum
        i;
        myresample = randsample(size(X,1),size(X,1),1);
        res(i,:)=myresample; %store resampling out of interest
        Xr=X(myresample,:); % define X for resampled subjects
        Yr=Y(myresample,:); % define X for resampled subjects
        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data
          
        temp=stats.W(:,1);%extract PLS1 weights
        newW=temp(x1); %order the newly obtained weights the same way as initial PLS 
        if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
            newW=-1*newW;
        end
        PLS1weights=[PLS1weights,newW];%store (ordered) weights from this bootstrap run
        
        temp=stats.W(:,2);%extract PLS2 weights
        newW=temp(x2); %order the newly obtained weights the same way as initial PLS 
        if corr(PLS2w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
            newW=-1*newW;
        end
        PLS2weights=[PLS2weights,newW]; %store (ordered) weights from this bootstrap run   

        temp=stats.W(:,3);%extract PLS3 weights
        newW=temp(x3); %order the newly obtained weights the same way as initial PLS 
        if corr(PLS3w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
            newW=-1*newW;
        end
        PLS3weights=[PLS3weights,newW]; %store (ordered) weights from this bootstrap run 
    end

    %get standard deviation of weights from bootstrap runs
    PLS1sw=std(PLS1weights');
    PLS2sw=std(PLS2weights');
    PLS3sw=std(PLS3weights');
    
    %get bootstrap weights
    temp1=PLS1w./PLS1sw';
    temp2=PLS2w./PLS2sw';
    temp3=PLS3w./PLS3sw';

    %order bootstrap weights (Z) and names of regions
    [Z1 ind1]=sort(temp1,'descend');
    PLS1=PLS1ids(ind1);
    geneindex1=geneindex1(ind1);
    [Z2 ind2]=sort(temp2,'descend');
    PLS2=PLS2ids(ind2);
    geneindex2=geneindex2(ind2);
    [Z3 ind3]=sort(temp3,'descend');
    PLS3=PLS3ids(ind3);
    geneindex3=geneindex3(ind3);

    %print out results
    % later use first column of these csv files for pasting into GOrilla (for
    % bootstrapped ordered list of genes) 
    fid1 = fopen(strcat('PLS genes results/subtype_', string(id), '_PLS1_geneWeights.csv'),'w');
    for i=1:length(genes)
      fprintf(fid1,'%s, %d, %f\n', PLS1{i},geneindex1(i), Z1(i));
    end
    fclose(fid1);

    fid2 = fopen(strcat('PLS genes results/subtype_', string(id), '_PLS2_geneWeights.csv'),'w');
    for i=1:length(genes)
      fprintf(fid2,'%s, %d, %f\n', PLS2{i},geneindex2(i), Z2(i));
    end
    fclose(fid2);

    fid3 = fopen(strcat('PLS genes results/subtype_', string(id), '_PLS3_geneWeights.csv'),'w');
    for i=1:length(genes)
      fprintf(fid3,'%s, %d, %f\n', PLS3{i},geneindex3(i), Z3(i));
    end
    fclose(fid3);   
end
