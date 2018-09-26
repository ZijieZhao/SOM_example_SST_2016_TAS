%% SOM tutorial/example in climate science - A SOM application on SST in 2016 over Tasman Sea.
% copy @ Zijie Zhao, University of Tasmania/University of Melbourne,
% zijiezhaomj@gmail.com

%% load the NOAA OI SST data for 2016
load('sst_2016');
sst_2016=double(sst_2016);
x=double(x);
y=double(y);

sst_2016(abs(sst_2016)>100)=nan;

%% arrange the data from 3D to 2D i.e. time * loc
sst_2d=[];
for i=1:size(sst_2016,3);
    data_here=sst_2016(:,:,i);
    sst_2d=[sst_2d;(data_here(:))'];
end

%% remove mean and std and store them for future application
[som_data,mu,sd]=zscore(sst_2d);

%% find the NAN in data and remove it. 
% Remember that we need to store the location of NAN since we need to know
% it when drawing the plot.

land_log=nansum(isnan(som_data))==size(som_data,1);
land_loc=find(nansum(isnan(som_data))==size(som_data,1));

som_data=som_data(:,~land_log);
sd=sd(~land_log);
mu=mu(~land_log);

%% determine the suitable SOM size for current data

% generate a number of SOM from size (1,2) to size (6,6).Actually, I
% normally do it from size (1,1) to size (10,10) in my own research, but it
% actually takes a lot of time to run so we just do it with a relatively
% small size in this small script.
for i=1:6
    for j=1:6
        
        if i*j~=1
        
        som=som_make(som_data,'munits',i*j,'msize',[i,j]);
        som_here=['som_' num2str(i) '_' num2str(j)];
        save(som_here,'som');
        clear som
        
        end
    end
end

% After generating and storing all SOMs, we need to determine which SOM
% size is the best one for current dataset.

% Here we use a correlation - based tool, which has been used in Zhao et
% al. (2018; in prep.). In plain language, As more nodes were included, 
% the generated patterns were reconstructed into a dataset with the same 
% size of original data by duplicating each patterns based on its allocated 
% temporal datasets and the correlation between these two datasets are calculated. 
% The final map size of the SOM was determined as the size in which the correlation 
% tends to be constant.

label_full=[];
cor_full=[];

for i=1:6
    for j=1:6
        if i*j~=1
        file_here=['som_' num2str(i) '_' num2str(j)];
        load(file_here)
        code=som.codebook;
        label=som_bmus(som,som_data);
        
        rep_data=code(label,:);
        
        cor_here=corr([rep_data(:) som_data(:)]);
        cor_here=cor_here(1,2);
        
        label_full=[label_full;[i j]];
        cor_full=[cor_full;cor_here];
        end
    end
end

[cor_full,i]=sort(cor_full);
label_full=label_full(i,:);

% draw the first difference of correlation in different map size.

figure('pos',[10 10 1000 1000]);

plot(1:length(diff(cor_full)),diff(cor_full),'linewidth',2);
xlim([1 35]);

label_text={};

for i=1:size(label_full,1);
    label_text{i}=['(' num2str(label_full(i,1)) ',' num2str(label_full(i,2)) ')'];
end

set(gca,'xtick',[5 8 10 15 19 25 30],'xticklabels',label_text([5 8 10 15 19 25 30]+1),'fontsize',16);

hold on

plot(8*ones(1000,1),linspace(0,0.035,1000),'r--');

% From this plot, we could see that, since size (4,1), the first difference
% of correlation in different map size tends to be constant and stationary,
% without so many gurgitation. So here we tend to choose size (4,1) as our
% determined map size. 

%% after determining that 4*1 is good, that's see the pattern presented in each node


load('som_4_1.mat');
label=som_bmus(som,som_data);

label_prec=[];

for i=1:4
    label_prec=[label_prec;nansum(label==i)./length(label)];
end

patterns=som.codebook.*sd+mu;

full_patterns=NaN(4,2898);
full_patterns(:,~land_log)=patterns;

index_plot=[1 2 3 4];
full_nodes=[1 1;2 1;3 1;4 1];

figure('pos',[10 10 1500 1500]);
h=tight_subplot(2,2,[0.05 0.01],[0.05 0.05],[0.05 0.07]);

m_proj('miller','lon',double([nanmin(x(:)) nanmax(x(:))]),'lat',double([nanmin(y(:)) nanmax(y(:))]));

for i=1:4
    data_here=full_patterns(i,:);
    data_here=reshape(data_here,69,42);
    
    
    axes(h(index_plot(i)));
    
    m_contourf(double(x),double(y),(data_here)',10:0.01:22,'linestyle','none');
    
    m_coast();
    if index_plot(i)~=3
        m_grid('xtick',[],'ytick',[],'fontsize',14);
    else
        m_grid('fontsize',14,'linestyle','none');
    end
    
    colormap(jet);
    
    text_here=['Node (' num2str(full_nodes(i,1)) ',' num2str(full_nodes(i,2)) '):' num2str(round(label_prec(i)*100,2)) '%'];
    title(text_here,'fontsize',16,'fontweight','bold');
    caxis([10 22]);
end

% From this plot, it could be determined that from Node (1,1) to Node
% (4,1), the SST over Tasman Sea tends to be higher. These four SOM
% patterns generally represent the major characteristics of SST in 2016
% over Tasman Sea. 

%% occurence time series
occurence_matrix=zeros(12,4);

date_full=datevec(datenum(2016,1,1):datenum(2016,12,31));

date_each_month={};

for i=1:12
    full=1:366;
    date_each_month{i}=full(date_full(:,2)==i);
end

for i=1:4
    label_here=find(label==i);
    
    temp_dist=[];
    
    for j=1:12
        add_here=nansum(ismember(label_here,date_each_month{j}));
        temp_dist=[temp_dist;add_here];
    end
    
    occurence_matrix(:,i)=temp_dist;
end

occurence_matrix=occurence_matrix./nansum(occurence_matrix,2);

occurence_table=table(occurence_matrix(:,1),occurence_matrix(:,2),occurence_matrix(:,3),occurence_matrix(:,4),'variablenames',...
    {'Node1','Node2','Node3','Node4'},'rownames',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});

occurence_table

% From this occurence_table, you can see that there is temporal variability
% among the SOM nodes. SST in Jul - Nov are mostly labelled into Node(1,1),
% while SST in Jun and Dec is labelled into Node(2,1), May is mostly
% labelled into Node(3,1) while Jan - Apr are labelled into Node(4,1).

% It means that, these determined periods could be potentially associated
% with the SST patterns shown in last figure. E.G. In Tasman Sea, during
% Jan to Apr in 2016, we expect to see the SST similar to Node(4,1).

% You may say: "wow, that's nice. But why don't we just use seasonally
% temporal average i.e. four average in four seasons, to represent the
% patterns in this year. It could be good too."

% That's good, but only when your dataset is relatively small. Considering
% that you have a dataset with 60 years, only using seasonal/monthly
% temporal average could not reveal some important characteristics and
% variabilities, such as ENSO patterns in SST. See Johnson, 2013.

%% Comments

% This is just a short note about how to simply apply SOM to climate
% science, but actually SOM has more things. Some further questions
% are,e.g.,
% how to apply SOM to discrete time series, how to explain SOM patterns,
% how to connect SOM based on one data with another data, what's the
% difference between SOM, PCA, EOF and other neural networks....A lot of
% things worth learning and I am free to chat and discuss. 




















        
        


    