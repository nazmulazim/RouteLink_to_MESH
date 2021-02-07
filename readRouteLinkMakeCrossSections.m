% AA = load ('D:\Project_Works\JTTI\NHD_Network\NHD_Geometry.csv');
clear all;
folderName = 'D:\Project_Works\JTTI\Florence_NC\Model\';
% folderName = 'D:\Project_Works\JTTI\ARBNM\Model\';
fileName = 'NetCDF_Florence.csv';
% fileName = 'Purged_Florence.csv';
% fileName = 'Route_link_Mat_A.csv';
% T = readtable('D:\Project_Works\JTTI\NHD_Network\NHD_Geometry.csv','PreserveVariableNames',true);
T = readtable([folderName fileName],'PreserveVariableNames',true);
A = table2array(T);
header = T.Properties.VariableNames;
clear T;

indexTo = find(contains(header,'to'));
indexLink = find(contains(header,'link'));
indexFrom = find(contains(header,'from_'));
indexBw = find(contains(header,'BtmWdth'));
indexDx = find(contains(header,'Length'));
indexTw = find(contains(header,'TopWdthCC'))-1;
indexTwCC = find(contains(header,'TopWdthCC'));
indexZ = find(contains(header,'alt'));
indexChSlp = find(contains(header,'ChSlp'));
indexMann = find(contains(header,'nCC'))-1;
indexSo = find(contains(header,'So'));


% lastNodeLink = 4185691;
% lastNodeLink = 9731454;
% lastNodeLink = 9733144;

lastNodeLink = 166737675; % for Florence

% lastNodeLink = 15208365; % for Amity River


% populating 'from' column
for i = 1:length(A(:,2))
    linkID = A(i,indexLink);
    fromValueIndex = A(:,indexTo)==linkID;
    
    if sum(fromValueIndex) == 0
        A(i,indexFrom)=0;
    else
        fromValue = A(A(:,indexTo)==linkID,indexLink);
        A(i,indexFrom)=fromValue(1);
    end
end

B = A;
% deleting all nodes that are not part of the network
for i = 1:length(A(:,2))
    fromID = A(i,indexFrom);
    if fromID == 0
        aa = i;
        bb = 1;
        clear indexCollect;
        nextLinkValue = A(i,indexLink);
        while nextLinkValue ~= lastNodeLink
            toValue = A(aa,indexTo);
            nextLinkValue = toValue;
            if sum(aa) == 0
                % the channel is not within the network
                B(indexCollect,:) = NaN;
                break;
            else
                indexCollect(bb) = aa;
                bb = bb + 1;
            end
            aa = find(A(:,indexLink)==toValue);
        end
    end
end

clear linkID aa bb fromID fromValue nextLinkValue toValue;

C=B;
C(isnan(C(:,1)),:)=[];
% The array C has got rid off all the links and nodes whose last node is
% not the given 'lastNodeLink'

[GC,GR] = groupcounts(C(:,indexTo));
% the GR element that has corresponding GC value more than 1 is the
% collection point of more than one channel


maxLinks = length(C(:,2));
maxNcomp = length(C(:,2));

nx1 = ones(1,maxNcomp)*NaN;
channelSequence = ones(1,maxNcomp)*NaN;

allLinks = ones(maxNcomp,maxLinks)*NaN;
allLinksTo = ones(maxNcomp,maxLinks)*NaN;

bwAll = ones(maxNcomp,maxLinks)*NaN;
TwAll = ones(maxNcomp,maxLinks)*NaN;
TwCCAll = ones(maxNcomp,maxLinks)*NaN;
ChSlpAll = ones(maxNcomp,maxLinks)*NaN;
dxAll = ones(maxNcomp,maxLinks)*NaN;
zAll = ones(maxNcomp,maxLinks)*NaN;
skAll = ones(maxNcomp,maxLinks)*NaN;

ndep = zeros(1,maxLinks);



j = 1;
for i = 1:length(C(:,2))
    fromID = C(i,indexFrom);
    fprintf([num2str(i) '\n']);
    linkID = C(i,indexLink);
    if fromID == 0 || sum(GC(GR==linkID)) > 1
        if fromID == 0
            ndep(j) = 0;
        elseif sum(GC(GR==linkID)) > 1
            ndep(j) = sum(GC(GR==linkID));
        end
        
        xx = 1;
        aa = i;
        bb = 1;
        while xx > 0
            allLinks(bb,j)=C(aa,indexLink);
            allLinksTo(bb,j)=C(aa,indexTo);
            dxAll(bb,j)=C(aa,indexDx);
            zAll(bb,j)=C(aa,indexZ);
            bwAll(bb,j)=C(aa,indexBw);
            TwAll(bb,j)=C(aa,indexTw);
            TwCCAll(bb,j)=C(aa,indexTwCC);
            ChSlpAll(bb,j) = C(aa,indexChSlp);
            skAll(bb,j) = 1/C(aa,indexMann);
        
            toValue = C(aa,indexTo);
            nextLinkValue = toValue;
            
            if nextLinkValue == lastNodeLink
                xx = 0;
                bb = bb + 1;
                
                channelSequence(j) = maxLinks;
                
                allLinks(bb,j)=C(aa,indexLink);
                allLinksTo(bb,j)=C(aa,indexTo);
                bwAll(bb,j)=C(aa,indexBw);
                TwAll(bb,j)=C(aa,indexTw);
                TwCCAll(bb,j)=C(aa,indexTwCC);
                ChSlpAll(bb,j) = C(aa,indexChSlp);
                skAll(bb,j) = 1/C(aa,indexMann);
                
                zAll(bb,j)=C(aa,indexZ)-C(aa,indexSo)*C(aa,indexDx);
                
                nx1(j) = bb;                
                j = j + 1;
                
            elseif GC(GR==nextLinkValue) > 1
                xx = 0;
                bb = bb + 1;
                
                allLinks(bb,j)=C(aa,indexLink);
                allLinksTo(bb,j)=C(aa,indexTo);
                bwAll(bb,j)=C(aa,indexBw);
                TwAll(bb,j)=C(aa,indexTw);
                TwCCAll(bb,j)=C(aa,indexTwCC);
                ChSlpAll(bb,j) = C(aa,indexChSlp);
                skAll(bb,j) = 1/C(aa,indexMann);
                
                nx1(j) = bb;
                aa = find(C(:,indexLink)==toValue);
                zAll(bb,j)=C(aa,indexZ);
                
                j = j + 1;
            else
                bb = bb+1;
                aa = find(C(:,indexLink)==toValue);
            end
        end
    end
end

nlinks = j-1;


% reducing all the array sizes
nx1(nlinks+1:end)=[];
ndep(nlinks+1:end)=[];
channelSequence(nlinks+1:end)=[];

allLinks(:,nlinks+1:end)=[];
allLinks(max(nx1)+1:end,:)=[];

allLinksTo(:,nlinks+1:end)=[];
allLinksTo(max(nx1)+1:end,:)=[];

dxAll(:,nlinks+1:end)=[];
dxAll(max(nx1)+1:end,:)=[];

zAll(:,nlinks+1:end)=[];
zAll(max(nx1)+1:end,:)=[];

bwAll(:,nlinks+1:end)=[];
bwAll(max(nx1)+1:end,:)=[];

TwAll(:,nlinks+1:end)=[];
TwAll(max(nx1)+1:end,:)=[];

TwCCAll(:,nlinks+1:end)=[];
TwCCAll(max(nx1)+1:end,:)=[];

ChSlpAll(:,nlinks+1:end)=[];
ChSlpAll(max(nx1)+1:end,:)=[];

skAll(:,nlinks+1:end)=[];
skAll(max(nx1)+1:end,:)=[];

%% Channel Sequence:
% At this level, we need to find the computation sequence
% of our channel. In MESH, the upstream channels need to be solved first.
% The sequence has to be in such an order so that all the upstream channels
% are solved before solving a channel 'j'
count = nlinks;
indexAA = find(channelSequence==maxLinks);
channelSequence(indexAA) = count;

uplinks=zeros(max(ndep),nlinks);
dslinks=ones(1,nlinks)*NaN;
ndepNew=ones(size(dslinks))*NaN;
boundaryConditionNew = zeros(2,nlinks);

previousChannelLastNode = allLinks(1,indexAA);
% count = count-1;
unsolvedLastNode=[];
aa = find(C(:,indexTo)==lastNodeLink); 
currentChannelLastNode = C(aa,indexLink);
unsolvedLastNode = [unsolvedLastNode currentChannelLastNode];
while count>0
        currentChannelLastNodeLocation = find(allLinks==unsolvedLastNode(1));
        pp = floor(min(currentChannelLastNodeLocation)/max(nx1))+1;
        if channelSequence(pp)==nlinks
            dslinks(pp)=0;
        end
        channelSequence(pp)=count;      
        previousChannelLastNodeLocation=find(allLinksTo==allLinks(1,pp));
        previousChannelLastNodeLocation=previousChannelLastNodeLocation(1:2:end);
        qq = floor(previousChannelLastNodeLocation/max(nx1))+1;
        dslinks(qq) = count;
        if size(qq)>0
            uplinks(1:size(qq),count)=qq;
        end
        count=count-1;
        currentChannelLastNode=allLinks(previousChannelLastNodeLocation);
        unsolvedLastNode = [unsolvedLastNode currentChannelLastNode'];
        unsolvedLastNode(1)=[];
end

uplinksNew = uplinks*0;
dslinksNew = ones(size(dslinks))*0;

for j = 1:nlinks
    qq=channelSequence==j;
    dslinksNew(j) = dslinks(qq);
    for i=1:length(uplinks(:,j))
        if uplinks(i,j)>0
            uplinksNew(i,j) = channelSequence(uplinks(i,j));
        end
    end
end
dslinksNew(nlinks) = 0;
clear uplinks dslinks qq pp;
        
nx1New = ones(1,nlinks)*NaN;
nx2New = ones(1,nlinks)*NaN;



%% writing cross sections, dx file, skk file
saveFolder = [folderName 'Geometry_RouteLink\'];
nx2=nx1; % keeping a copy of old variables
for j= 1:nlinks
    mkdir ([saveFolder '\CS' num2str(channelSequence(j))]);
    % interpolate new cross sections if there are only two sections in a river reach
    if nx1(j)==2
        
        % Changing all branches that has only 1 link will be converted to n links; i.e. n+1 nodes
        n = 3; 
        
        Bw = bwAll(1,j);
        Tw = TwAll(1,j);
        z  = 1./ChSlpAll(1,j);
        TwCC = TwCCAll(1,j);
        bed = zAll(1,j);
        bfd =  (Tw - Bw)/(2.0*z);
    
        x1 = [-0.5*TwCC -0.5*TwCC -0.5*Tw -0.5*Bw 0.5*Bw 0.5*Tw 0.5*TwCC 0.5*TwCC] +0.5*TwCC;
        y1 = [4*bfd     bfd       bfd     0       0      bfd    bfd      4*bfd   ]+bed;
        
        Bw = bwAll(2,j);
        Tw = TwAll(2,j);
        z  = 1./ChSlpAll(2,j);
        TwCC = TwCCAll(2,j);
        bed = zAll(2,j);
        bfd =  (Tw - Bw)/(2.0*z);
    
        x2 = [-0.5*TwCC -0.5*TwCC -0.5*Tw -0.5*Bw 0.5*Bw 0.5*Tw 0.5*TwCC 0.5*TwCC] +0.5*TwCC;
        y2 = [4*bfd     bfd       bfd     0       0      bfd    bfd      4*bfd   ]+bed;
        
        [P,Q]=interpolateXsections(x1,y1,x2,y2,n);
        
        nx1(j) = n+1;
        
        for i = 1: nx1(j)
            fileNameNew = [saveFolder 'CS' num2str(channelSequence(j)) '\Test_' sprintf('%04d', i) ];

            newX = P(:,i);
            newY = Q(:,i);

            fid = fopen(fileNameNew,'w');
            if i==1
                fprintf(fid,'%s\t%s\t%s\n','x','y',['link ' num2str(allLinks(1,j))]);
            elseif i == nx1(j)
                fprintf(fid,'%s\t%s\t%s\n','x','y',['link ' num2str(allLinks(2,j))]);
            else
                fprintf(fid,'%s\t%s\t%s\n','x','y',['link interpol ' num2str(allLinks(1,j)) ' and ' num2str(allLinks(2,j))]);
            end
            for k = 1:length(newX)
               fprintf(fid,'%12.3f\t%12.3f\n',newX(k),newY(k));
            end
            fclose(fid);
        end
        
        fileNameNew = [saveFolder 'dx' sprintf('%04d', channelSequence(j)) ];    
        fid = fopen(fileNameNew,'w');
        dxAll(1:n,j)=dxAll(1,j)/n;
        for k = 1:n+1
            if isnan(dxAll(k,j))
                dxAll(k,j)=0;
            end
            fprintf(fid,'%12.3f\t%12.3f\n',sum(dxAll(1:k,j))-dxAll(k,j),dxAll(k,j));
        end
        fclose(fid);
        
        fileNameNew = [saveFolder 'skk' sprintf('%04d', channelSequence(j)) ];    
        fid = fopen(fileNameNew,'w');
        skAll(1:n+1,j)=skAll(1,j);
        for k = 1:nx1(j)
            fprintf(fid,'%12.3f\t%12.3f\n',sum(dxAll(1:k,j))-dxAll(k,j),skAll(k,j));
        end
        fclose(fid);
        
    else
    
        for i = 1: nx1(j)
            fileNameNew = [saveFolder 'CS' num2str(channelSequence(j)) '\Test_' sprintf('%04d', i) ];

            Bw = bwAll(i,j);
            Tw = TwAll(i,j);
            z  = 1./ChSlpAll(i,j);
            TwCC = TwCCAll(i,j);
            bed = zAll(i,j);

            bfd =  (Tw - Bw)/(2.0*z);

            newX = [-0.5*TwCC -0.5*TwCC -0.5*Tw -0.5*Bw 0.5*Bw 0.5*Tw 0.5*TwCC 0.5*TwCC] +0.5*TwCC;
            newY = [4*bfd     bfd       bfd     0       0      bfd    bfd      4*bfd   ];

            fid = fopen(fileNameNew,'w');
            fprintf(fid,'%s\t%s\t%s\n','x','y',['link ' num2str(allLinks(i,j))]);
            for k = 1:length(newX)
               fprintf(fid,'%12.3f\t%12.3f\n',newX(k),newY(k)+bed);
            end
            fclose(fid);
        end
    
        fileNameNew = [saveFolder 'dx' sprintf('%04d', channelSequence(j)) ];    
        fid = fopen(fileNameNew,'w');
        for k = 1:nx1(j)
            if isnan(dxAll(k,j))
                dxAll(k,j)=0;
            end
            fprintf(fid,'%12.3f\t%12.3f\n',sum(dxAll(1:k,j))-dxAll(k,j),dxAll(k,j));
        end
        fclose(fid);

        fileNameNew = [saveFolder 'skk' sprintf('%04d', channelSequence(j)) ];    
        fid = fopen(fileNameNew,'w');
        for k = 1:nx1(j)
            fprintf(fid,'%12.3f\t%12.3f\n',sum(dxAll(1:k,j))-dxAll(k,j),skAll(k,j));
        end
        fclose(fid);
    end

end
%% correcting the sequence in some matiics
for j=1:nlinks
    qq=channelSequence==j;
    ndepNew(j) = ndep(qq);    
    nx1New(j) = nx1(qq);
    nx2New(j) = nx2(qq);
    if ndepNew(j) == 0
        boundaryConditionNew(1,j)=2;
    end
    if dslinksNew(j) == 0
        boundaryConditionNew(2,j)=1;
    end
end

clear nx1 xx bb bed bfd Bw count currentChannelLastNode currentChannelLastNodeLocation ...
    fromID fromValueIndex indexAA j k previousChannelLastNode previousChannelLastNodeLocation qq ...
    toValue Tw TwCC unsolvedLastNode xx z ndep;

%% write network file
fileNameNew = [folderName '\Geometry_purged\Network_file_' num2str(nlinks) ];
fid = fopen(fileNameNew,'w');
fprintf(fid,'%s\n','No of upstream channels at each river reach');
for j = 1:nlinks
    fprintf(fid,'%d\t%d\n',j, ndepNew(j));
end

fprintf(fid,'%s\n','which links are at the U/S of the current link j');
for j = 1:nlinks
    if sum(uplinksNew(:,j))>0
        string = num2str(uplinksNew(:,j)');
    else
        string = '0';
    end
    fprintf(fid,'%d\t%s\n',j,string);
end

fprintf(fid,'%s\n','which links are at the D/S of the current link j');
for j = 1:nlinks
    fprintf(fid,'%d\t%d\n',j, dslinksNew(j));
end

fprintf(fid,'%s\n','U/S and D/S boundary conditions of the link j: First column = 2	means Q boundary at u/s; 2nd column = 1 means WL boundary at d/s. Value = 0 means internal boundary');
for j = 1:nlinks
    fprintf(fid,'%d\t%s\n',j,num2str(boundaryConditionNew(:,j)'));
end
fclose(fid);

%% write input file

initialWaterDepth = ones(1,nlinks)*0.05;
initialDischarge = ones(1,nlinks)*0.01;

fileNameNew = [folderName '\Geometry_purged\input_file_' num2str(nlinks) ];
fid = fopen(fileNameNew,'w');
fprintf(fid,'%s\n','10.0 =: dtini           Time step dt (in seconds)');
fprintf(fid,'%s\n','25.0 =: dxini           Spatial step dx								--- parameter not used');
fprintf(fid,'%s\n','0.0 =: t0               starting time (in hours) // 0  // 78888 -- Time corresponds to the time origin in hours');
fprintf(fid,'%s\n','100 =: tfin             Final time (in hours)');
fprintf(fid,'%d\t%s\n',nlinks,'=: nlinks    number of river reaches in the model');
fprintf(fid,'%s\t%s\n',num2str(nx1New),'=: ncomp        number of nodes in all river reaches');
fprintf(fid,'%s\n','1.0 =: phi              source term treatment (0:explicit, 1:implicit)');
fprintf(fid,'%s\n','1.0 =: theta         ?');
fprintf(fid,'%s\n','1.0 =: thetas        ?');
fprintf(fid,'%s\n','1.0 =: thesinv       ?');
fprintf(fid,'%s\n','0. =: alfa2         emp parameter for artificial diffusion');
fprintf(fid,'%s\n','0. =: alfa4         maximum value for artificial diffusion');
fprintf(fid,'%s\n','1.0 =: f             SI Units = 1');
fprintf(fid,'%s\n','80 =: skk           1/Mannings n roughness');
fprintf(fid,'%s\t%s\n',num2str(initialWaterDepth),'=: yy            initial value of water depth (m)');
fprintf(fid,'%s\t%s\n',num2str(initialDischarge),'=: qq        initial value of uniform flow condition (m3s-1)');
fprintf(fid,'%s\n','0.8 =: cfl          max courant number for optimal optimization of dt');
fprintf(fid,'%s\n','0 =: ots           optimize dt (1:apply opt, !=1:user specify)		--- parameter not used');
fprintf(fid,'%s\n','0.0 =: yw            weir height									--- parameter not used');
fprintf(fid,'%s\n','20.0 =: bw            weir width									--- parameter not used');
fprintf(fid,'%s\n','1.1 =: w             weir coefficient (unused)						--- parameter not used');
fprintf(fid,'%s\n','1 =: option        DS imposed condition (1:y, 2:q, 3:rating curve)	--- parameter not used');
fprintf(fid,'%s\n','0.1 =: yn            water level of DS end							--- parameter not used');
fprintf(fid,'%s\n','0.0085 =: qn            water discharge of DS end					--- parameter not used');
fprintf(fid,'%s\n','700 =: igate         spatial index of gate location (unused)		--- parameter not used');
for j=1:nlinks
    fprintf(fid,'%s\t%s\n',[saveFolder 'CS' num2str(j) '\Test_'], ['=: xSection_path' num2str(j)]);
end
for j=1:nlinks
    fprintf(fid,'%s\t%s\n',[saveFolder 'skk' num2str(j,'%04d')], ['=: Mannings-stricklers path' num2str(j)]);
end
fprintf(fid,'%d\t%s\n',sum(boundaryConditionNew(1,:)==2),'=: 		no of upstream boundaries');
for j=1:nlinks
    if boundaryConditionNew(1,j)==2
        %fprintf(fid,'%s\t%s\n',[saveFolder 'Boundaries\up_Q_' num2str(j)],['=: upstream_path ' num2str(j)]);
        fprintf(fid,'%s\t%s\n',[saveFolder 'Boundaries\up_Q'],['=: upstream_path ' num2str(j)]);
    end
end
fprintf(fid,'%d\t%s\n',sum(boundaryConditionNew(2,:)==1),'=: 		no of downstream boundaries');
for j=1:nlinks
    if boundaryConditionNew(2,j)==1
        %fprintf(fid,'%s\t%s\n',[saveFolder 'Boundaries\dn_wl_' num2str(j)],['=: downstream_path ' num2str(j)]);
        fprintf(fid,'%s\t%s\n',[saveFolder 'Boundaries\dn_wl'],['=: downstream_path ' num2str(j)]);
    end
end
for j=1:nlinks
    fprintf(fid,'%s\t%s\n',[saveFolder 'Q_SK_Tables\'],['=: Q-SK Table path ' num2str(j)]);
end
for j=1:nlinks
    fprintf(fid,'%s\t%s\n',[saveFolder 'dx' num2str(j,'%04d')],['=: dx_path ' num2str(j)]);
end
for j=1:nlinks
    fprintf(fid,'%s\t%s\n',[saveFolder 'Lateral_' num2str(j,'%04d') '\'],['=: Lateral flow path ' num2str(j)]);
end
fprintf(fid,'%s\t%s\n',[saveFolder 'output\'],'=: output_path');
fprintf(fid,'%s\n','0 =: option_dsbc         2 is drawdown to critical, 1 is normal, and 0 is constant water level downstream');
fprintf(fid,'%s\n','1000 =: maxTableLength         maximum number of data in each cross section file');
fprintf(fid,'%s\n','201 =: nel         number of line in each cross section attribute table');
fprintf(fid,'%s\n','10.0 =: timesDepth		multipliyer of depths that will be allowed to be flooded at a section');
fprintf(fid,'%s\t%s\n',[saveFolder 'input\'],'=: other input path');
fprintf(fid,'%s\n','5000 =: boundaryFileMaxEntry 	Max No of entry in the input boundary file');
fprintf(fid,'%s\n','300 =: Result Saving interval (in seconds)');
% string = num2str(zeros(1,nlinks));
string = num2str(nx2New-1);
fprintf(fid,'%s\t%s\n',string,'=: No of lateral flow inputs to the system');
for j=1:nlinks
    string = num2str(1:nx2New(j)-1);
    fprintf(fid,'%s\t%s\n',string,['=: all the first nodes where a lateral flow starts at reach ' num2str(j)]);
end
for j=1:nlinks
    string = num2str(ones(1,nx2New(j)-1));
    fprintf(fid,'%s\t%s\n',string,['=: Lateral flow type for reach ' num2str(j) '=: Type = 1 = time series; Type 2 = flow as a function of upstream flow']);
end
for j=1:nlinks
    string = num2str(ones(1,nx2New(j)-1));
    fprintf(fid,'%s\t%s\n',string,['=: no of x-secs at the downstream that the lateral flow is applied for reach ' num2str(j)]);
end
string = num2str(zeros(1,nlinks));
fprintf(fid,'%s\t%s\n',string,'=: No of Q-Sk multiplier table');
for j=1:nlinks
    fprintf(fid,'%s\t%s\n','',[':= all the starting nodes under each table for reach ' num2str(j)]);
    fprintf(fid,'%s\t%s\n','',[':= all the ending nodes under each table for reach ' num2str(j)]);
end
fprintf(fid,'%s\t%s\n',[folderName 'Geometry_purged\Network_file_' num2str(nlinks) ],'=: network connectivity file');


%% writing bwAll, TwAll, TwCCAll, ChSlpAll
bwAllNew = ones(max(nx1New),nlinks)*NaN;
TwAllNew = ones(max(nx1New),nlinks)*NaN;
TwCCAllNew = ones(max(nx1New),nlinks)*NaN;
ChSlpAllNew = ones(max(nx1New),nlinks)*NaN;
skAllNew = ones(max(nx1New),nlinks)*NaN;
for j=1:nlinks
    qq=channelSequence==j;
    bwAllNew(:,j) = bwAll(:,qq);    
    TwAllNew(:,j) = TwAll(:,qq);  
    TwCCAllNew(:,j) = TwCCAll(:,qq);  
    ChSlpAllNew(:,j) = ChSlpAll(:,qq);
    skAllNew(:,j) = skAll(:,qq);
    
    if nx2New(j) == 2
        bwAllNew(1:n+1,j)=bwAllNew(1,j);
        TwAllNew(1:n+1,j)=TwAllNew(1,j);
        TwCCAllNew(1:n+1,j)=TwCCAllNew(1,j);
        ChSlpAllNew(1:n+1,j)=ChSlpAllNew(1,j);
        skAllNew(1:n+1,j)=skAllNew(1,j);
    end
end

fileNameNew = [folderName '\Geometry_RouteLink\bwAllNew' ];
fid = fopen(fileNameNew,'w');
for j = 1:nlinks
    string = num2str(bwAllNew(1:nx1New(j),j)');
    fprintf(fid,'%s\n',string );
end
fclose(fid);

fileNameNew = [folderName '\Geometry_RouteLink\TwAllNew' ];
fid = fopen(fileNameNew,'w');
for j = 1:nlinks
    string = num2str(TwAllNew(1:nx1New(j),j)');
    fprintf(fid,'%s\n',string );
end
fclose(fid);

fileNameNew = [folderName '\Geometry_RouteLink\TwCCAllNew' ];
fid = fopen(fileNameNew,'w');
for j = 1:nlinks
    string = num2str(TwCCAllNew(1:nx1New(j),j)');
    fprintf(fid,'%s\n',string );
end
fclose(fid);

fileNameNew = [folderName '\Geometry_RouteLink\ChSlpAllNew' ];
fid = fopen(fileNameNew,'w');
for j = 1:nlinks
    string = num2str(ChSlpAllNew(1:nx1New(j),j)');
    fprintf(fid,'%s\n',string );
end
fclose(fid);

fileNameNew = [folderName '\Geometry_RouteLink\skAllNew' ];
fid = fopen(fileNameNew,'w');
for j = 1:nlinks
    string = num2str(skAllNew(1:nx1New(j),j)');
    fprintf(fid,'%s\n',string );
end
fclose(fid);