clear all
load traceRealCellLocations.mat

numStates2D=20;   %excluding state 0
centerCoordinate=[37.762 -122.43];

leftUserEachSlot = zeros(length(cellOfUsers(:,1)), 1);
totalUserEachSlot = zeros(length(cellOfUsers(:,1)), 1);

for user = 1:length(cellOfUsers(1,:))
    prevCell = 0;
    for slot = length(cellOfUsers(:,1)):-1:1
        if cellOfUsers(slot, user) ~= 0
            totalUserEachSlot(slot) = totalUserEachSlot(slot)+1;
        end
        if (prevCell ~= 0) && (slot < length(cellOfUsers(:,1))) && (cellOfUsers(slot, user) ~= prevCell)
            leftUserEachSlot(slot+1) = leftUserEachSlot(slot+1) + 1;
        end
        prevCell = cellOfUsers(slot, user);
    end
end



availResourceTransFactorStore=[];
availResourceMigrationFactorStore=[];
maxUserEachCloudStore=[];
numCellsWithCloudStore=[];

avrCostNeverMigrateStore=[];
avrCostAlwaysMigrateStore=[];
avrCostMyopicStore=[];
avrCostThPolicyStore=[];
stdCostNeverMigrateStore=[];
stdCostAlwaysMigrateStore=[];
stdCostMyopicStore=[];
stdCostThPolicyStore=[];

avrCostGainOverNeverMigrateStore=[];
avrCostGainOverAlwaysMigrateStore=[];
avrCostGainOverMyopicStore=[];
stdCostGainOverNeverMigrateStore=[];
stdCostGainOverAlwaysMigrateStore=[];
stdCostGainOverMyopicStore=[];
maxCostGainOverNeverMigrateStore=[];
maxCostGainOverAlwaysMigrateStore=[];
maxCostGainOverMyopicStore=[];
minCostGainOverNeverMigrateStore=[];
minCostGainOverAlwaysMigrateStore=[];
minCostGainOverMyopicStore=[];

avrFirstMigrationStateThPolicyStore=[];
stdFirstMigrationStateThPolicyStore=[];
maxFirstMigrationStateThPolicyStore=[];
minFirstMigrationStateThPolicyStore=[];

for availResourceTransFactor=1.5
    for availResourceMigrationFactor=1.5
%         for maxUserEachCloud=[10,20,30,50,70,100,150,200]
        for maxUserEachCloud=50
%             for numCellsWithCloud=[20,30,50,70,100,150,200]
            for numCellsWithCloud=100


                %Init for algorithm
                gamma=0.9;  %discount factor
                powerFactor=0.8;


                Use2D=1;
                notUseValueIteration=1;
                notUsePolicyIteration=1;
                notFindActualValues=1;


                proportionalFactorMigrateEachTimeslot=zeros(length(numberOfUsersInCell(:,1)),1);
                proportionalFactorTransEachTimeslot=zeros(length(numberOfUsersInCell(:,1)),1);
                constFactorMigrateEachTimeslot=zeros(length(numberOfUsersInCell(:,1)),1);
                constFactorTransEachTimeslot=zeros(length(numberOfUsersInCell(:,1)),1);  

                actionsMyopicEachTimeslot=zeros(length(numberOfUsersInCell(:,1)),numStates2DTotal);
                actionsThPolicyEachTimeslot=zeros(length(numberOfUsersInCell(:,1)),numStates2D+1);

                PEachSlot=zeros(numStates2D+1,numStates2D+1,length(numberOfUsersInCell(:,1)));
                valuesThPolicyEachSlot=zeros(numStates2D+1,length(numberOfUsersInCell(:,1)));

                for timeslot=length(numberOfUsersInCell(:,1)):-1:1   %timeslot in data is saved in reversed order
                    %1:length(numberOfUsersInCell(:,1))

                    P2D = 1 - leftUserEachSlot(timeslot)/totalUserEachSlot(timeslot);

                    %recalculate costs
                    maxUsers=max(totalUsers);
                    tmpAvailResourceTrans=availResourceTransFactor*maxUsers;
                    tmpAvailResourceMigrate=availResourceMigrationFactor*maxUsers;
                    
                    %non-constant cost
                    proportionalFactorMigrate=-1/(1-totalUsers(timeslot)/(tmpAvailResourceTrans));
                    proportionalFactorTrans=-1/(1-totalUsers(timeslot)/(tmpAvailResourceTrans));
                    constFactorMigrate=1*1/(1-totalUsers(timeslot)/(tmpAvailResourceMigrate))-proportionalFactorMigrate;
                    constFactorTrans=0-proportionalFactorTrans;    

                    %constant costs
                    % proportionalFactorMigrate=0;
                    % proportionalFactorTrans=0;
                    % constFactorMigrate=1*1/(1-totalUsers(timeslot)/(tmpAvailResourceMigrate));
                    % constFactorTrans=1/(1-totalUsers(timeslot)/(tmpAvailResourceTrans));    
                    
                    disp(['** get actions -- timeslot: ', num2str(timeslot), ' availResourceTransFactor=', num2str(availResourceTransFactor),...
                        ' availResourceMigrationFactor=', num2str(availResourceMigrationFactor),...
                        ' maxUserEachCloud=', num2str(maxUserEachCloud),...
                        ' numCellsWithCloud=', num2str(numCellsWithCloud)])

                    supressOutput=1;
                    run([pwd, '/algorithms.m'])
                    if algReturn==1
                        return
                    end

                    PEachSlot(:,:,timeslot)=P;
                    valuesThPolicyEachSlot(:,timeslot)=valuesThPolicy;

                    %store data for current timeslot
                    proportionalFactorMigrateEachTimeslot(timeslot)=proportionalFactorMigrate;
                    proportionalFactorTransEachTimeslot(timeslot)=proportionalFactorTrans;
                    constFactorMigrateEachTimeslot(timeslot)=constFactorMigrate;
                    constFactorTransEachTimeslot(timeslot)=constFactorTrans;

                    actionsMyopicEachTimeslot(timeslot,:)=actionsMyopic;    %No transpose here due to specification in algorithms.m
                    actionsThPolicyEachTimeslot(timeslot,:)=actionsThPolicy';    

                end


                %% find actual cost from traces

                %initialization
                tmpStep=length(coordinatesCells2D(:,1))/numCellsWithCloud;

                cellsWithCloud=zeros(length(coordinatesCells2D(:,1)),1); %if cell has cloud, value is 1, otherwise 0
                i=tmpStep;
                while(round(i)<=length(coordinatesCells2D(:,1)))
                    cellsWithCloud(round(i))=1;
                    i=i+tmpStep;
                end
                cellsWithCloudIndexes=find(cellsWithCloud==1);

                serviceLocationNeverMigrate=zeros(1,length(cellOfUsers(1,:))); %if zero, service not placed; if vehicle shut off, then set service location to zero
                serviceLocationAlwaysMigrate=zeros(1,length(cellOfUsers(1,:))); %if zero, service not placed; if vehicle shut off, then set service location to zero
                serviceLocationMyopic=zeros(1,length(cellOfUsers(1,:))); %if zero, service not placed; if vehicle shut off, then set service location to zero
                serviceLocationThPolicy=zeros(1,length(cellOfUsers(1,:))); %if zero, service not placed; if vehicle shut off, then set service location to zero

                traceOneTimeslotCostNeverMigrate=zeros(length(numberOfUsersInCell(:,1)), length(cellOfUsers(1,:)));
                traceOneTimeslotCostAlwaysMigrate=zeros(length(numberOfUsersInCell(:,1)), length(cellOfUsers(1,:)));
                traceOneTimeslotCostMyopic=zeros(length(numberOfUsersInCell(:,1)), length(cellOfUsers(1,:)));
                traceOneTimeslotCostThPolicy=zeros(length(numberOfUsersInCell(:,1)), length(cellOfUsers(1,:)));

                numUsersEachCloudNeverMigrate=zeros(length(numberOfUsersInCell(:,1)),length(coordinatesCells2D(:,1)));
                numUsersEachCloudAlwaysMigrate=zeros(length(numberOfUsersInCell(:,1)),length(coordinatesCells2D(:,1)));
                numUsersEachCloudMyopic=zeros(length(numberOfUsersInCell(:,1)),length(coordinatesCells2D(:,1)));
                numUsersEachCloudThPolicy=zeros(length(numberOfUsersInCell(:,1)),length(coordinatesCells2D(:,1)));

                for timeslot=length(numberOfUsersInCell(:,1)):-1:1
                    disp(['** get costs -- timeslot: ', num2str(timeslot), ' availResourceTransFactor=', num2str(availResourceTransFactor),...
                        ' availResourceMigrationFactor=', num2str(availResourceMigrationFactor),...
                        ' maxUserEachCloud=', num2str(maxUserEachCloud),...
                        ' numCellsWithCloud=', num2str(numCellsWithCloud)])

                    %% find migration locations
                    tmpMigrateToCellNeverMigrate=zeros(length(cellOfUsers(1,:)),1);
                    tmpMigrateToCellAlwaysMigrate=zeros(length(cellOfUsers(1,:)),1);
                    tmpMigrateToCellMyopic=zeros(length(cellOfUsers(1,:)),1);
                    tmpMigrateToCellThPolicy=zeros(length(cellOfUsers(1,:)),1);

                    tmpMigrateToCellNeverMigrateOrig=zeros(length(cellOfUsers(1,:)),1);
                    tmpMigrateToCellAlwaysMigrateOrig=zeros(length(cellOfUsers(1,:)),1);
                    tmpMigrateToCellMyopicOrig=zeros(length(cellOfUsers(1,:)),1);
                    tmpMigrateToCellThPolicyOrig=zeros(length(cellOfUsers(1,:)),1);

                    for i=1:length(cellOfUsers(1,:))
                        if cellOfUsers(timeslot,i)~=0

                            %% never migrate except at border
                            tmpPrevLocationIndex=serviceLocationNeverMigrate(i);
                            if serviceLocationNeverMigrate(i)==0  %service initialization
                                tmpMigrateToCell=cellOfUsers(timeslot,i);

                                tmpMigrateToCellNeverMigrateOrig(i)=tmpMigrateToCell;
                                %** block for finding closest cell with cloud (initialization) **
                                tmpDiff=coordinatesCells2D(cellsWithCloudIndexes,:)-ones(length(cellsWithCloudIndexes),1)*coordinatesCells2D(tmpMigrateToCell,:);
                                tmpDist=tmpDiff(:,1).^2+tmpDiff(:,2).^2;
                                [~, tmpDistMinIndex]=min(tmpDist);
                                tmpMigrateToCell=cellsWithCloudIndexes(tmpDistMinIndex);
                                %****

                            else
                                tmpNorm1=norm(coordinatesCells2D(cellOfUsers(timeslot,i),:)-coordinatesCells2D(serviceLocationNeverMigrate(i),:))/cellDist;
                                if abs(tmpNorm1-round(tmpNorm1))<0.00000001
                                    tmpNorm1=round(tmpNorm1);
                                else
                                    tmpNorm1=ceil(tmpNorm1);
                                end

                                if tmpNorm1>=numStates2D
                                    tmpMigrateToCell=cellOfUsers(timeslot,i);
                                else
                                    tmpMigrateToCell=serviceLocationNeverMigrate(i);
                                end

                                tmpMigrateToCellNeverMigrateOrig(i)=tmpMigrateToCell;
                                %** block for finding closest cell with cloud (migration, SPECIFIC FOR NEVER/ALWAYS MIGRATE) **
                                tmpDiff=coordinatesCells2D(cellsWithCloudIndexes,:)-ones(length(cellsWithCloudIndexes),1)*coordinatesCells2D(tmpMigrateToCell,:);
                                tmpDist=tmpDiff(:,1).^2+tmpDiff(:,2).^2;
                                [~, tmpDistMinIndex]=min(tmpDist);
                                tmpMigrateToCell=cellsWithCloudIndexes(tmpDistMinIndex);
                                %****

                            end
                            tmpMigrateToCellNeverMigrate(i)=tmpMigrateToCell;





                            %% always migrate
                            tmpPrevLocationIndex=serviceLocationAlwaysMigrate(i);
                            if serviceLocationAlwaysMigrate(i)==0
                                tmpMigrateToCell=cellOfUsers(timeslot,i);

                                tmpMigrateToCellAlwaysMigrateOrig(i)=tmpMigrateToCell;
                                %** block for finding closest cell with cloud (initialization) **
                                tmpDiff=coordinatesCells2D(cellsWithCloudIndexes,:)-ones(length(cellsWithCloudIndexes),1)*coordinatesCells2D(tmpMigrateToCell,:);
                                tmpDist=tmpDiff(:,1).^2+tmpDiff(:,2).^2;
                                [~, tmpDistMinIndex]=min(tmpDist);
                                tmpMigrateToCell=cellsWithCloudIndexes(tmpDistMinIndex);
                                %****

                            else
                                tmpMigrateToCell=cellOfUsers(timeslot,i);

                                tmpMigrateToCellAlwaysMigrateOrig(i)=tmpMigrateToCell;
                                %** block for finding closest cell with cloud (migration, SPECIFIC FOR NEVER/ALWAYS MIGRATE) **
                                tmpDiff=coordinatesCells2D(cellsWithCloudIndexes,:)-ones(length(cellsWithCloudIndexes),1)*coordinatesCells2D(tmpMigrateToCell,:);
                                tmpDist=tmpDiff(:,1).^2+tmpDiff(:,2).^2;
                                [~, tmpDistMinIndex]=min(tmpDist);
                                tmpMigrateToCell=cellsWithCloudIndexes(tmpDistMinIndex);
                                %****
                            end
                            tmpMigrateToCellAlwaysMigrate(i)=tmpMigrateToCell;




                            %% myopic
                            tmpPrevLocationIndex=serviceLocationMyopic(i);
                            if serviceLocationMyopic(i)==0
                                tmpMigrateToCell=cellOfUsers(timeslot,i);

                                tmpMigrateToCellMyopicOrig(i)=tmpMigrateToCell;
                                %** block for finding closest cell with cloud (initialization) **
                                tmpDiff=coordinatesCells2D(cellsWithCloudIndexes,:)-ones(length(cellsWithCloudIndexes),1)*coordinatesCells2D(tmpMigrateToCell,:);
                                tmpDist=tmpDiff(:,1).^2+tmpDiff(:,2).^2;
                                [~, tmpDistMinIndex]=min(tmpDist);
                                tmpMigrateToCell=cellsWithCloudIndexes(tmpDistMinIndex);
                                %****                        


                            else
                                tmpOffset=dsearchn(coordinatesCells2D,centerCoordinate+coordinatesCells2D(cellOfUsers(timeslot,i),:)-coordinatesCells2D(serviceLocationMyopic(i),:)); 

                                if tmpOffset==1  %state zero
                                    tmpNewLocationCoordinate=coordinatesCells2D(cellOfUsers(timeslot,i),:);
                                else
                                    tmpNewLocationCoordinate=coordinatesCells2D(cellOfUsers(timeslot,i),:)-(coordinatesCells2D(actionsMyopicEachTimeslot(timeslot,tmpOffset),:)-centerCoordinate);
                                
                                end

                                tmpMigrateToCell = dsearchn(coordinatesCells2D,tmpNewLocationCoordinate);                

                                tmpMigrateToCellMyopicOrig(i)=tmpMigrateToCell;

                                
                                
                                
                                
                                %** block for finding closest cell with cloud (migration, SPECIFIC FOR MYOPIC) **
                                tmpDiff1=coordinatesCells2D(cellsWithCloudIndexes,:)-ones(length(cellsWithCloudIndexes),1)*coordinatesCells2D(cellOfUsers(timeslot,i),:);
                                tmpDist1=sqrt(tmpDiff1(:,1).^2+tmpDiff1(:,2).^2)/cellDist;

                                tmpAbs=abs(tmpDist1-round(tmpDist1));
                                tmpDist1=(tmpAbs<0.00000001).*round(tmpDist1)+(tmpAbs>=0.00000001).*ceil(tmpDist1);
                                
                                cellsWithCloudIndexesSub=cellsWithCloudIndexes(tmpDist1<numStates2D); %only consider less than numStates2D hops from the user
                                
                                tmpDiff1=coordinatesCells2D(cellsWithCloudIndexesSub,:)-ones(length(cellsWithCloudIndexesSub),1)*coordinatesCells2D(tmpPrevLocationIndex,:);
                                tmpDiff2=coordinatesCells2D(cellsWithCloudIndexesSub,:)-ones(length(cellsWithCloudIndexesSub),1)*coordinatesCells2D(cellOfUsers(timeslot,i),:);
                                tmpDist1=sqrt(tmpDiff1(:,1).^2+tmpDiff1(:,2).^2)/cellDist;
                                tmpDist2=sqrt(tmpDiff2(:,1).^2+tmpDiff2(:,2).^2)/cellDist;

                                tmpAbs=abs(tmpDist1-round(tmpDist1));
                                tmpDist1=(tmpAbs<0.00000001).*round(tmpDist1)+(tmpAbs>=0.00000001).*ceil(tmpDist1);

                                tmpAbs=abs(tmpDist2-round(tmpDist2));
                                tmpDist2=(tmpAbs<0.00000001).*round(tmpDist2)+(tmpAbs>=0.00000001).*ceil(tmpDist2);

                                tmpCost1=constFactorMigrateEachTimeslot(timeslot)+proportionalFactorMigrateEachTimeslot(timeslot)*powerFactor.^tmpDist1;
                                tmpCost2=constFactorTransEachTimeslot(timeslot)+proportionalFactorTransEachTimeslot(timeslot)*powerFactor.^tmpDist2;
                                tmpCostForDistMin=tmpCost1.*(cellsWithCloudIndexesSub~=tmpPrevLocationIndex)...
                                    + tmpCost2.*(cellOfUsers(timeslot,i)~=cellsWithCloudIndexesSub);

                                [~,tmpDistMinIndex]=min(tmpCostForDistMin);

                                tmpMigrateToCell=cellsWithCloudIndexesSub(tmpDistMinIndex);
                                %****
                                

                            end  
                            tmpMigrateToCellMyopic(i)=tmpMigrateToCell;





                            %% modified policy iteration
                            tmpPrevLocationIndex=serviceLocationThPolicy(i);
                            if serviceLocationThPolicy(i)==0
                                tmpMigrateToCell=cellOfUsers(timeslot,i);

                                tmpMigrateToCellThPolicyOrig(i)=tmpMigrateToCell;
                                %** block for finding closest cell with cloud (initialization) **
                                tmpDiff=coordinatesCells2D(cellsWithCloudIndexes,:)-ones(length(cellsWithCloudIndexes),1)*coordinatesCells2D(tmpMigrateToCell,:);
                                tmpDist=tmpDiff(:,1).^2+tmpDiff(:,2).^2;
                                [~, tmpDistMinIndex]=min(tmpDist);
                                tmpMigrateToCell=cellsWithCloudIndexes(tmpDistMinIndex);
                                %****                        


                            else
                                tmpNorm1=norm(coordinatesCells2D(cellOfUsers(timeslot,i),:)-coordinatesCells2D(serviceLocationThPolicy(i),:))/cellDist;
                                if abs(tmpNorm1-round(tmpNorm1))<0.00000001
                                    tmpNorm1=round(tmpNorm1);
                                else
                                    tmpNorm1=ceil(tmpNorm1);
                                end

                                if tmpNorm1==0
                                    tmpNewLocationCoordinate=coordinatesCells2D(cellOfUsers(timeslot,i),:);
                                else
                                    tmpNewLocationCoordinate=(coordinatesCells2D(cellOfUsers(timeslot,i),:)-coordinatesCells2D(serviceLocationThPolicy(i),:))*...
                                        (1-(actionsThPolicyEachTimeslot(timeslot,min(tmpNorm1,numStates2D)+1)-1)/tmpNorm1)+coordinatesCells2D(serviceLocationThPolicy(i),:);
                                end

                                tmpMigrateToCell = dsearchn(coordinatesCells2D,tmpNewLocationCoordinate);

                                tmpMigrateToCellThPolicyOrig(i)=tmpMigrateToCell;
                                
                                %** block for finding closest cell with cloud (migration, SPECIFIC FOR PROPOSED) **
                                tmpDiff1=coordinatesCells2D(cellsWithCloudIndexes,:)-ones(length(cellsWithCloudIndexes),1)*coordinatesCells2D(tmpPrevLocationIndex,:);
                                tmpDiff2=coordinatesCells2D(cellsWithCloudIndexes,:)-ones(length(cellsWithCloudIndexes),1)*coordinatesCells2D(cellOfUsers(timeslot,i),:);
                                tmpDist1=sqrt(tmpDiff1(:,1).^2+tmpDiff1(:,2).^2)/cellDist;
                                tmpDist2=sqrt(tmpDiff2(:,1).^2+tmpDiff2(:,2).^2)/cellDist;

                                tmpAbs=abs(tmpDist1-round(tmpDist1));
                                tmpDist1=(tmpAbs<0.00000001).*round(tmpDist1)+(tmpAbs>=0.00000001).*ceil(tmpDist1);

                                tmpAbs=abs(tmpDist2-round(tmpDist2));
                                tmpDist2=(tmpAbs<0.00000001).*round(tmpDist2)+(tmpAbs>=0.00000001).*ceil(tmpDist2);

                                tmpDist2=min(tmpDist2,numStates2D); %neglect distances larger than numStates2D, may need to CHANGE later
                                tmpCostForDistMin=(constFactorMigrateEachTimeslot(timeslot)+proportionalFactorMigrateEachTimeslot(timeslot)*powerFactor.^tmpDist1).*(cellsWithCloudIndexes~=tmpPrevLocationIndex)...
                                    + (constFactorTransEachTimeslot(timeslot)+proportionalFactorTransEachTimeslot(timeslot)*powerFactor.^tmpDist2).*(cellOfUsers(timeslot,i)~=cellsWithCloudIndexes)...
                                    +gamma*PEachSlot(tmpDist2+1,:,timeslot)*valuesThPolicyEachSlot(:,timeslot);

                                [~,tmpDistMinIndex]=min(tmpCostForDistMin);

                                tmpMigrateToCell=cellsWithCloudIndexes(tmpDistMinIndex);
                                %****
                                
                                
                                
                                
                            end            
                            tmpMigrateToCellThPolicy(i)=tmpMigrateToCell;



                        end

                    end


                    %% find number of users in each cloud
                    for i=1:length(coordinatesCells2D(:,1))
                        for j=1:4
                            if j==1
                                numUsersEachCloudNeverMigrate(timeslot,i)=length(find(tmpMigrateToCellNeverMigrate==i));
                            elseif j==2
                                numUsersEachCloudAlwaysMigrate(timeslot,i)=length(find(tmpMigrateToCellAlwaysMigrate==i));
                            elseif j==3
                                numUsersEachCloudMyopic(timeslot,i)=length(find(tmpMigrateToCellMyopic==i));
                            elseif j==4
                                numUsersEachCloudThPolicy(timeslot,i)=length(find(tmpMigrateToCellThPolicy==i));
                            end
                        end
                    end

                    %% find and reassign users exceeding cloud capacity
                    for j=1:4
                        if j==1 %never migrate
                            tmpExceedIndex=find(numUsersEachCloudNeverMigrate(timeslot,:)>maxUserEachCloud);
                            for i=tmpExceedIndex
                                tmpUsersInExceedingCloud=find(tmpMigrateToCellNeverMigrate==i);

                                %plug in cost function here
                                tmpDiff=coordinatesCells2D(tmpMigrateToCellNeverMigrateOrig(tmpUsersInExceedingCloud),:)-coordinatesCells2D(tmpMigrateToCellNeverMigrate(tmpUsersInExceedingCloud),:);
                                tmpDist=tmpDiff(:,1).^2+tmpDiff(:,2).^2;

                                [~,tmpSortedIndexes]=sort(tmpDist,'descend'); %try to remove those with largest cost

                                tmpNumUsersToGoOut=numUsersEachCloudNeverMigrate(timeslot,i)-maxUserEachCloud;

                                for k=1:tmpNumUsersToGoOut
                                    tmpUserIndex=tmpUsersInExceedingCloud(tmpSortedIndexes(k));

                                    %plug in cost function here
                                    tmpDiff=coordinatesCells2D(cellsWithCloudIndexes,:)-ones(length(cellsWithCloudIndexes),1)*coordinatesCells2D(tmpMigrateToCellNeverMigrateOrig(tmpUserIndex),:);
                                    tmpDist=tmpDiff(:,1).^2+tmpDiff(:,2).^2;
                                    [~, tmpDistSortIndex]=sort(tmpDist); %try to place in first cloud (with lowest cost) that still has space

                                    tmpFound=0;
                                    for l=1:length(tmpDistSortIndex);
                                        if numUsersEachCloudNeverMigrate(timeslot,cellsWithCloudIndexes(l))<maxUserEachCloud
                                            tmpFound=1;
                                            tmpMigrateToCellNeverMigrate(tmpUserIndex)=cellsWithCloudIndexes(l);
                                            numUsersEachCloudNeverMigrate(timeslot,cellsWithCloudIndexes(l))=numUsersEachCloudNeverMigrate(timeslot,cellsWithCloudIndexes(l))+1;
                                            numUsersEachCloudNeverMigrate(timeslot,i)=numUsersEachCloudNeverMigrate(timeslot,i)-1;
                                            break
                                        end
                                    end
                                    if tmpFound==0
                                        disp('ERROR: infeasible assignment: never migrate')
                                        return
                                    end


                                end

                            end

                        elseif j==2 %always migrate
                            tmpExceedIndex=find(numUsersEachCloudAlwaysMigrate(timeslot,:)>maxUserEachCloud);
                            for i=tmpExceedIndex
                                tmpUsersInExceedingCloud=find(tmpMigrateToCellAlwaysMigrate==i);

                                %plug in cost function here
                                tmpDiff=coordinatesCells2D(tmpMigrateToCellAlwaysMigrateOrig(tmpUsersInExceedingCloud),:)-coordinatesCells2D(tmpMigrateToCellAlwaysMigrate(tmpUsersInExceedingCloud),:);
                                tmpDist=tmpDiff(:,1).^2+tmpDiff(:,2).^2;

                                [~,tmpSortedIndexes]=sort(tmpDist,'descend'); %try to remove those with largest cost

                                tmpNumUsersToGoOut=numUsersEachCloudAlwaysMigrate(timeslot,i)-maxUserEachCloud;

                                for k=1:tmpNumUsersToGoOut
                                    tmpUserIndex=tmpUsersInExceedingCloud(tmpSortedIndexes(k));

                                    %plug in cost function here
                                    tmpDiff=coordinatesCells2D(cellsWithCloudIndexes,:)-ones(length(cellsWithCloudIndexes),1)*coordinatesCells2D(tmpMigrateToCellAlwaysMigrateOrig(tmpUserIndex),:);
                                    tmpDist=tmpDiff(:,1).^2+tmpDiff(:,2).^2;
                                    [~, tmpDistSortIndex]=sort(tmpDist); %try to place in first cloud (with lowest cost) that still has space

                                    tmpFound=0;
                                    for l=1:length(tmpDistSortIndex);
                                        if numUsersEachCloudAlwaysMigrate(timeslot,cellsWithCloudIndexes(l))<maxUserEachCloud
                                            tmpFound=1;
                                            tmpMigrateToCellAlwaysMigrate(tmpUserIndex)=cellsWithCloudIndexes(l);
                                            numUsersEachCloudAlwaysMigrate(timeslot,cellsWithCloudIndexes(l))=numUsersEachCloudAlwaysMigrate(timeslot,cellsWithCloudIndexes(l))+1;
                                            numUsersEachCloudAlwaysMigrate(timeslot,i)=numUsersEachCloudAlwaysMigrate(timeslot,i)-1;
                                            break
                                        end
                                    end
                                    if tmpFound==0
                                        disp('ERROR: infeasible assignment: always migrate')
                                        return
                                    end


                                end

                            end


                        elseif j==3 %myopic
                            tmpExceedIndex=find(numUsersEachCloudMyopic(timeslot,:)>maxUserEachCloud);
                            for i=tmpExceedIndex
                                tmpUsersInExceedingCloud=find(tmpMigrateToCellMyopic==i);

                                %plug in cost function here
                                
                                
                                tmpDiff1=zeros(length(tmpUsersInExceedingCloud),2);
                                tmpDiff1NonZeroIndexes=find(serviceLocationMyopic(tmpUsersInExceedingCloud)~=0);
                                tmpDiff1(tmpDiff1NonZeroIndexes,:)=coordinatesCells2D(tmpMigrateToCellMyopic(tmpUsersInExceedingCloud(tmpDiff1NonZeroIndexes)),:)-coordinatesCells2D(serviceLocationMyopic(tmpUsersInExceedingCloud(tmpDiff1NonZeroIndexes)),:);
                                tmpDiff2=coordinatesCells2D(tmpMigrateToCellMyopic(tmpUsersInExceedingCloud),:)-coordinatesCells2D(cellOfUsers(timeslot,tmpUsersInExceedingCloud),:);
                                tmpDist1=sqrt(tmpDiff1(:,1).^2+tmpDiff1(:,2).^2)/cellDist;
                                tmpDist2=sqrt(tmpDiff2(:,1).^2+tmpDiff2(:,2).^2)/cellDist;

                                tmpAbs=abs(tmpDist1-round(tmpDist1));
                                tmpDist1=(tmpAbs<0.00000001).*round(tmpDist1)+(tmpAbs>=0.00000001).*ceil(tmpDist1);

                                tmpAbs=abs(tmpDist2-round(tmpDist2));
                                tmpDist2=(tmpAbs<0.00000001).*round(tmpDist2)+(tmpAbs>=0.00000001).*ceil(tmpDist2);

                                tmpCost1=constFactorMigrateEachTimeslot(timeslot)+proportionalFactorMigrateEachTimeslot(timeslot)*powerFactor.^tmpDist1;
                                tmpCost2=constFactorTransEachTimeslot(timeslot)+proportionalFactorTransEachTimeslot(timeslot)*powerFactor.^tmpDist2;
                                tmpCostForDistMin=tmpCost1.*(tmpMigrateToCellMyopic(tmpUsersInExceedingCloud)~=serviceLocationMyopic(tmpUsersInExceedingCloud)')...
                                    + tmpCost2.*(cellOfUsers(timeslot,tmpUsersInExceedingCloud)'~=tmpMigrateToCellMyopic(tmpUsersInExceedingCloud));

                                [~,tmpSortedIndexes]=sort(tmpCostForDistMin,'descend'); %try to remove those with largest cost
                                
                                %end cost function
                                
                                

                                tmpNumUsersToGoOut=numUsersEachCloudMyopic(timeslot,i)-maxUserEachCloud;

                                for k=1:tmpNumUsersToGoOut
                                    tmpUserIndex=tmpUsersInExceedingCloud(tmpSortedIndexes(k));

                                    %plug in cost function here

                                    
                                    tmpDiff1=coordinatesCells2D(cellsWithCloudIndexes,:)-ones(length(cellsWithCloudIndexes),1)*coordinatesCells2D(cellOfUsers(timeslot,tmpUserIndex),:);
                                    tmpDist1=sqrt(tmpDiff1(:,1).^2+tmpDiff1(:,2).^2)/cellDist;

                                    tmpAbs=abs(tmpDist1-round(tmpDist1));
                                    tmpDist1=(tmpAbs<0.00000001).*round(tmpDist1)+(tmpAbs>=0.00000001).*ceil(tmpDist1);

                                    cellsWithCloudIndexesSub=cellsWithCloudIndexes(tmpDist1<numStates2D); %only consider less than numStates2D hops from the user

                                    if serviceLocationMyopic(tmpUserIndex)==0
                                        tmpDiff1=ones(length(cellsWithCloudIndexesSub),1)*[0,0];
                                    else
                                        tmpDiff1=coordinatesCells2D(cellsWithCloudIndexesSub,:)-ones(length(cellsWithCloudIndexesSub),1)*coordinatesCells2D(serviceLocationMyopic(tmpUserIndex),:);
                                    end
                                    tmpDiff2=coordinatesCells2D(cellsWithCloudIndexesSub,:)-ones(length(cellsWithCloudIndexesSub),1)*coordinatesCells2D(cellOfUsers(timeslot,tmpUserIndex),:);
                                    tmpDist1=sqrt(tmpDiff1(:,1).^2+tmpDiff1(:,2).^2)/cellDist;
                                    tmpDist2=sqrt(tmpDiff2(:,1).^2+tmpDiff2(:,2).^2)/cellDist;

                                    tmpAbs=abs(tmpDist1-round(tmpDist1));
                                    tmpDist1=(tmpAbs<0.00000001).*round(tmpDist1)+(tmpAbs>=0.00000001).*ceil(tmpDist1);

                                    tmpAbs=abs(tmpDist2-round(tmpDist2));
                                    tmpDist2=(tmpAbs<0.00000001).*round(tmpDist2)+(tmpAbs>=0.00000001).*ceil(tmpDist2);

                                    tmpCost1=constFactorMigrateEachTimeslot(timeslot)+proportionalFactorMigrateEachTimeslot(timeslot)*powerFactor.^tmpDist1;
                                    tmpCost2=constFactorTransEachTimeslot(timeslot)+proportionalFactorTransEachTimeslot(timeslot)*powerFactor.^tmpDist2;
                                    tmpCostForDistMin=tmpCost1.*(cellsWithCloudIndexesSub~=serviceLocationMyopic(tmpUserIndex)')...
                                        + tmpCost2.*(cellOfUsers(timeslot,tmpUserIndex)'~=cellsWithCloudIndexesSub);
                                    

                                    [~, tmpDistSortIndex]=sort(tmpCostForDistMin); %try to place in first cloud (with lowest cost) that still has space
                                    %end cost function
                                    

                                    tmpFound=0;
                                    for l=1:length(tmpDistSortIndex);
                                        if numUsersEachCloudMyopic(timeslot,cellsWithCloudIndexes(l))<maxUserEachCloud
                                            tmpFound=1;
                                            tmpMigrateToCellMyopic(tmpUserIndex)=cellsWithCloudIndexes(l);
                                            numUsersEachCloudMyopic(timeslot,cellsWithCloudIndexes(l))=numUsersEachCloudMyopic(timeslot,cellsWithCloudIndexes(l))+1;
                                            numUsersEachCloudMyopic(timeslot,i)=numUsersEachCloudMyopic(timeslot,i)-1;
                                            break
                                        end
                                    end
                                    if tmpFound==0
                                        disp('ERROR: infeasible assignment: myopic')
                                        return
                                    end


                                end

                            end



                        elseif j==4 %proposed
                            tmpExceedIndex=find(numUsersEachCloudThPolicy(timeslot,:)>maxUserEachCloud);
                            for i=tmpExceedIndex
                                tmpUsersInExceedingCloud=find(tmpMigrateToCellThPolicy==i);

                                %plug in cost function here
                                tmpDiff1=zeros(length(tmpUsersInExceedingCloud),2);
                                tmpDiff1NonZeroIndexes=find(serviceLocationThPolicy(tmpUsersInExceedingCloud)~=0);
                                tmpDiff1(tmpDiff1NonZeroIndexes,:)=coordinatesCells2D(tmpMigrateToCellThPolicy(tmpUsersInExceedingCloud(tmpDiff1NonZeroIndexes)),:)-coordinatesCells2D(serviceLocationThPolicy(tmpUsersInExceedingCloud(tmpDiff1NonZeroIndexes)),:);
                                tmpDiff2=coordinatesCells2D(tmpMigrateToCellThPolicy(tmpUsersInExceedingCloud),:)-coordinatesCells2D(cellOfUsers(timeslot,tmpUsersInExceedingCloud),:);
                                tmpDist1=sqrt(tmpDiff1(:,1).^2+tmpDiff1(:,2).^2)/cellDist;
                                tmpDist2=sqrt(tmpDiff2(:,1).^2+tmpDiff2(:,2).^2)/cellDist;

                                tmpAbs=abs(tmpDist1-round(tmpDist1));
                                tmpDist1=(tmpAbs<0.00000001).*round(tmpDist1)+(tmpAbs>=0.00000001).*ceil(tmpDist1);

                                tmpAbs=abs(tmpDist2-round(tmpDist2));
                                tmpDist2=(tmpAbs<0.00000001).*round(tmpDist2)+(tmpAbs>=0.00000001).*ceil(tmpDist2);

                                tmpDist2=min(tmpDist2,numStates2D); %neglect distances larger than numStates2D, may need to CHANGE later
                                tmpCost1=constFactorMigrateEachTimeslot(timeslot)+proportionalFactorMigrateEachTimeslot(timeslot)*powerFactor.^tmpDist1;
                                tmpCost2=constFactorTransEachTimeslot(timeslot)+proportionalFactorTransEachTimeslot(timeslot)*powerFactor.^tmpDist2;
                                tmpCostForDistMin=tmpCost1.*(tmpMigrateToCellThPolicy(tmpUsersInExceedingCloud)~=serviceLocationThPolicy(tmpUsersInExceedingCloud)')...
                                    + tmpCost2.*(cellOfUsers(timeslot,tmpUsersInExceedingCloud)'~=tmpMigrateToCellThPolicy(tmpUsersInExceedingCloud))...
                                    +gamma*PEachSlot(tmpDist2+1,:,timeslot)*valuesThPolicyEachSlot(:,timeslot);

                                [~,tmpSortedIndexes]=sort(tmpCostForDistMin,'descend'); %try to remove those with largest cost

                                %end cost function
                                
                                
                                tmpNumUsersToGoOut=numUsersEachCloudThPolicy(timeslot,i)-maxUserEachCloud;

                                for k=1:tmpNumUsersToGoOut
                                    tmpUserIndex=tmpUsersInExceedingCloud(tmpSortedIndexes(k));

                                    %plug in cost function here
                                    tmpDiff1=coordinatesCells2D(cellsWithCloudIndexes,:)-ones(length(cellsWithCloudIndexes),1)*coordinatesCells2D(cellOfUsers(timeslot,tmpUserIndex),:);
                                    tmpDist1=sqrt(tmpDiff1(:,1).^2+tmpDiff1(:,2).^2)/cellDist;

                                    tmpAbs=abs(tmpDist1-round(tmpDist1));
                                    tmpDist1=(tmpAbs<0.00000001).*round(tmpDist1)+(tmpAbs>=0.00000001).*ceil(tmpDist1);

                                    cellsWithCloudIndexesSub=cellsWithCloudIndexes(tmpDist1<numStates2D); %only consider less than numStates2D hops from the user

                                    if serviceLocationThPolicy(tmpUserIndex)==0
                                        tmpDiff1=ones(length(cellsWithCloudIndexesSub),1)*[0,0];
                                    else
                                        tmpDiff1=coordinatesCells2D(cellsWithCloudIndexesSub,:)-ones(length(cellsWithCloudIndexesSub),1)*coordinatesCells2D(serviceLocationThPolicy(tmpUserIndex),:);
                                    end
                                    tmpDiff2=coordinatesCells2D(cellsWithCloudIndexesSub,:)-ones(length(cellsWithCloudIndexesSub),1)*coordinatesCells2D(cellOfUsers(timeslot,tmpUserIndex),:);
                                    tmpDist1=sqrt(tmpDiff1(:,1).^2+tmpDiff1(:,2).^2)/cellDist;
                                    tmpDist2=sqrt(tmpDiff2(:,1).^2+tmpDiff2(:,2).^2)/cellDist;

                                    tmpAbs=abs(tmpDist1-round(tmpDist1));
                                    tmpDist1=(tmpAbs<0.00000001).*round(tmpDist1)+(tmpAbs>=0.00000001).*ceil(tmpDist1);

                                    tmpAbs=abs(tmpDist2-round(tmpDist2));
                                    tmpDist2=(tmpAbs<0.00000001).*round(tmpDist2)+(tmpAbs>=0.00000001).*ceil(tmpDist2);

                                    tmpDist2=min(tmpDist2,numStates2D); %neglect distances larger than numStates2D, may need to CHANGE later
                                    tmpCost1=constFactorMigrateEachTimeslot(timeslot)+proportionalFactorMigrateEachTimeslot(timeslot)*powerFactor.^tmpDist1;
                                    tmpCost2=constFactorTransEachTimeslot(timeslot)+proportionalFactorTransEachTimeslot(timeslot)*powerFactor.^tmpDist2;
                                    tmpCostForDistMin=tmpCost1.*(cellsWithCloudIndexesSub~=serviceLocationThPolicy(tmpUserIndex)')...
                                        + tmpCost2.*(cellOfUsers(timeslot,tmpUserIndex)'~=cellsWithCloudIndexesSub)...
                                        +gamma*PEachSlot(tmpDist2+1,:,timeslot)*valuesThPolicyEachSlot(:,timeslot);
                                    

                                    [~, tmpDistSortIndex]=sort(tmpCostForDistMin); %try to place in first cloud (with lowest cost) that still has space
                                    %end cost function

                                    tmpFound=0;
                                    for l=1:length(tmpDistSortIndex);
                                        if numUsersEachCloudThPolicy(timeslot,cellsWithCloudIndexes(l))<maxUserEachCloud
                                            tmpFound=1;
                                            tmpMigrateToCellThPolicy(tmpUserIndex)=cellsWithCloudIndexes(l);
                                            numUsersEachCloudThPolicy(timeslot,cellsWithCloudIndexes(l))=numUsersEachCloudThPolicy(timeslot,cellsWithCloudIndexes(l))+1;
                                            numUsersEachCloudThPolicy(timeslot,i)=numUsersEachCloudThPolicy(timeslot,i)-1;
                                            break
                                        end
                                    end
                                    if (tmpFound==0)
                                        disp('ERROR: infeasible assignment: proposed')
                                        return
                                    end


                                end

                            end



                        end
                    end


                    %% find costs
                    for i=1:length(cellOfUsers(1,:))
                        if cellOfUsers(timeslot,i)==0
                            serviceLocationNeverMigrate(i)=0;
                            serviceLocationAlwaysMigrate(i)=0;
                            serviceLocationMyopic(i)=0;
                            serviceLocationThPolicy(i)=0;
                            traceOneTimeslotCostNeverMigrate(timeslot,i)=0;
                            traceOneTimeslotCostAlwaysMigrate(timeslot,i)=0;
                            traceOneTimeslotCostMyopic(timeslot,i)=0;
                            traceOneTimeslotCostThPolicy(timeslot,i)=0;
                        else
                            for j=1:4
                                if j==1
                                    tmpPrevLocationIndex=serviceLocationNeverMigrate(i);
                                    tmpMigrateToCell=tmpMigrateToCellNeverMigrate(i);
                                elseif j==2
                                    tmpPrevLocationIndex=serviceLocationAlwaysMigrate(i);
                                    tmpMigrateToCell=tmpMigrateToCellAlwaysMigrate(i);
                                elseif j==3
                                    tmpPrevLocationIndex=serviceLocationMyopic(i);
                                    tmpMigrateToCell=tmpMigrateToCellMyopic(i);
                                elseif j==4
                                    tmpPrevLocationIndex=serviceLocationThPolicy(i);
                                    tmpMigrateToCell=tmpMigrateToCellThPolicy(i);
                                end

                                %** block for cost calculation **
                                if tmpPrevLocationIndex==0
                                    tmpNorm1=0;
                                else
                                    tmpNorm1=norm(coordinatesCells2D(tmpMigrateToCell,:)-coordinatesCells2D(tmpPrevLocationIndex,:))/cellDist;
                                    if abs(tmpNorm1-round(tmpNorm1))<0.00000001
                                        tmpNorm1=round(tmpNorm1);
                                    else
                                        tmpNorm1=ceil(tmpNorm1);
                                    end
                                end

                                tmpNorm2=norm(coordinatesCells2D(cellOfUsers(timeslot,i),:)-coordinatesCells2D(tmpMigrateToCell,:))/cellDist;
                                if abs(tmpNorm2-round(tmpNorm2))<0.00000001
                                    tmpNorm2=round(tmpNorm2);
                                else
                                    tmpNorm2=ceil(tmpNorm2);
                                end

                                tmpCost=(constFactorMigrateEachTimeslot(timeslot)+proportionalFactorMigrateEachTimeslot(timeslot)*powerFactor^tmpNorm1)*(tmpMigrateToCell~=tmpPrevLocationIndex)...
                                    + (constFactorTransEachTimeslot(timeslot)+proportionalFactorTransEachTimeslot(timeslot)*powerFactor^tmpNorm2)*(cellOfUsers(timeslot,i)~=tmpMigrateToCell);

                                %****

                                if j==1
                                    traceOneTimeslotCostNeverMigrate(timeslot,i)=tmpCost;
                                    serviceLocationNeverMigrate(i)=tmpMigrateToCell;
                                elseif j==2
                                    traceOneTimeslotCostAlwaysMigrate(timeslot,i)=tmpCost;
                                    serviceLocationAlwaysMigrate(i)=tmpMigrateToCell;
                                elseif j==3
                                    traceOneTimeslotCostMyopic(timeslot,i)=tmpCost;
                                    serviceLocationMyopic(i)=tmpMigrateToCell;
                                elseif j==4
                                    traceOneTimeslotCostThPolicy(timeslot,i)=tmpCost;
                                    serviceLocationThPolicy(i)=tmpMigrateToCell;
                                end


                            end

                        end
                    end

                end

                availResourceTransFactorStore=[availResourceTransFactorStore; availResourceTransFactor];
                availResourceMigrationFactorStore=[availResourceMigrationFactorStore; availResourceMigrationFactor];        
                maxUserEachCloudStore=[maxUserEachCloudStore, maxUserEachCloud];
                numCellsWithCloudStore=[numCellsWithCloudStore, numCellsWithCloud];        

                avrCostNeverMigrateStore=[avrCostNeverMigrateStore; mean(sum(traceOneTimeslotCostNeverMigrate,2)./totalUsers)];
                avrCostAlwaysMigrateStore=[avrCostAlwaysMigrateStore; mean(sum(traceOneTimeslotCostAlwaysMigrate,2)./totalUsers)];
                avrCostMyopicStore=[avrCostMyopicStore; mean(sum(traceOneTimeslotCostMyopic,2)./totalUsers)];
                avrCostThPolicyStore=[avrCostThPolicyStore; mean(sum(traceOneTimeslotCostThPolicy,2)./totalUsers)];
                stdCostNeverMigrateStore=[stdCostNeverMigrateStore; std(sum(traceOneTimeslotCostNeverMigrate,2)./totalUsers)];
                stdCostAlwaysMigrateStore=[stdCostAlwaysMigrateStore; std(sum(traceOneTimeslotCostAlwaysMigrate,2)./totalUsers)];
                stdCostMyopicStore=[stdCostMyopicStore; std(sum(traceOneTimeslotCostMyopic,2)./totalUsers)];
                stdCostThPolicyStore=[stdCostThPolicyStore; std(sum(traceOneTimeslotCostThPolicy,2)./totalUsers)];


                tmpGainOverNeverMigrate=sum(traceOneTimeslotCostThPolicy,2)./sum(traceOneTimeslotCostNeverMigrate,2);
                tmpGainOverAlwaysMigrate=sum(traceOneTimeslotCostThPolicy,2)./sum(traceOneTimeslotCostAlwaysMigrate,2);
                tmpGainOverMyopic=sum(traceOneTimeslotCostThPolicy,2)./sum(traceOneTimeslotCostMyopic,2);
                tmpGainOverNeverMigrate(isnan(tmpGainOverNeverMigrate))=[];
                tmpGainOverAlwaysMigrate(isnan(tmpGainOverAlwaysMigrate))=[];
                tmpGainOverMyopic(isnan(tmpGainOverMyopic))=[];

                avrCostGainOverNeverMigrateStore=[avrCostGainOverNeverMigrateStore; mean(tmpGainOverNeverMigrate)];
                avrCostGainOverAlwaysMigrateStore=[avrCostGainOverAlwaysMigrateStore; mean(tmpGainOverAlwaysMigrate)];
                avrCostGainOverMyopicStore=[avrCostGainOverMyopicStore; mean(tmpGainOverMyopic)];
                stdCostGainOverNeverMigrateStore=[stdCostGainOverNeverMigrateStore; std(tmpGainOverNeverMigrate)];
                stdCostGainOverAlwaysMigrateStore=[stdCostGainOverAlwaysMigrateStore; std(tmpGainOverAlwaysMigrate)];
                stdCostGainOverMyopicStore=[stdCostGainOverMyopicStore; std(tmpGainOverMyopic)];
                maxCostGainOverNeverMigrateStore=[maxCostGainOverNeverMigrateStore; max(tmpGainOverNeverMigrate)];
                maxCostGainOverAlwaysMigrateStore=[maxCostGainOverAlwaysMigrateStore; max(tmpGainOverAlwaysMigrate)];
                maxCostGainOverMyopicStore=[maxCostGainOverMyopicStore; max(tmpGainOverMyopic)];
                minCostGainOverNeverMigrateStore=[minCostGainOverNeverMigrateStore; min(tmpGainOverNeverMigrate)];
                minCostGainOverAlwaysMigrateStore=[minCostGainOverAlwaysMigrateStore; min(tmpGainOverAlwaysMigrate)];
                minCostGainOverMyopicStore=[minCostGainOverMyopicStore; min(tmpGainOverMyopic)];

                tmpFirstMigratesVector=zeros(length(numberOfUsersInCell(:,1)),1);
                for timeslot=1:length(numberOfUsersInCell(:,1))
                    tmpFirstMigratesVector(timeslot)=find(actionsThPolicyEachTimeslot(timeslot,:)<[1:numStates2D+1],1)-1;
                end

                avrFirstMigrationStateThPolicyStore=[avrFirstMigrationStateThPolicyStore; mean(tmpFirstMigratesVector)];
                stdFirstMigrationStateThPolicyStore=[stdFirstMigrationStateThPolicyStore; std(tmpFirstMigratesVector)];
                maxFirstMigrationStateThPolicyStore=[maxFirstMigrationStateThPolicyStore; max(tmpFirstMigratesVector)];
                minFirstMigrationStateThPolicyStore=[minFirstMigrationStateThPolicyStore; min(tmpFirstMigratesVector)];



                disp(['availResourceTransFactorStore ', num2str(availResourceTransFactorStore')])
                disp(['availResourceMigrationFactorStore ', num2str(availResourceMigrationFactorStore')])

                disp(['avrCostGainOverNeverMigrateStore ', num2str(avrCostGainOverNeverMigrateStore')])
                disp(['avrCostGainOverAlwaysMigrateStore ', num2str(avrCostGainOverAlwaysMigrateStore')])
                disp(['avrCostGainOverMyopicStore ', num2str(avrCostGainOverMyopicStore')])
                disp(['stdCostGainOverNeverMigrateStore ', num2str(stdCostGainOverNeverMigrateStore')])
                disp(['stdCostGainOverAlwaysMigrateStore ', num2str(stdCostGainOverAlwaysMigrateStore')])
                disp(['stdCostGainOverMyopicStore ', num2str(stdCostGainOverMyopicStore')])
                disp(['maxCostGainOverNeverMigrateStore ', num2str(maxCostGainOverNeverMigrateStore')])
                disp(['maxCostGainOverAlwaysMigrateStore ', num2str(maxCostGainOverAlwaysMigrateStore')])
                disp(['maxCostGainOverMyopicStore ', num2str(maxCostGainOverMyopicStore')])
                disp(['minCostGainOverNeverMigrateStore ', num2str(minCostGainOverNeverMigrateStore')])
                disp(['minCostGainOverAlwaysMigrateStore ', num2str(minCostGainOverAlwaysMigrateStore')])
                disp(['minCostGainOverMyopicStore ', num2str(minCostGainOverMyopicStore')])

                disp(['avrFirstMigrationStateThPolicyStore ', num2str(avrFirstMigrationStateThPolicyStore')])
                disp(['stdFirstMigrationStateThPolicyStore ', num2str(stdFirstMigrationStateThPolicyStore')])
                disp(['maxFirstMigrationStateThPolicyStore ', num2str(maxFirstMigrationStateThPolicyStore')])
                disp(['minFirstMigrationStateThPolicyStore ', num2str(minFirstMigrationStateThPolicyStore')])

                disp('')
            end
        end
    end
end

xPointsEpochs=timeMax:-updateTimeStep:timeMin;
tmp=datevec((xPointsEpochs-8*3600)./86400 + datenum(1970,1,1));
xPoints=datenum(tmp);


figure(1)
plot(xPoints,sum(traceOneTimeslotCostNeverMigrate,2)./totalUsers,'k',xPoints,sum(traceOneTimeslotCostAlwaysMigrate,2)./totalUsers,'g',...
    xPoints,sum(traceOneTimeslotCostMyopic,2)./totalUsers,'b',xPoints,sum(traceOneTimeslotCostThPolicy,2)./totalUsers,'r',...
    xPoints,ones(1,length(xPoints))*mean(sum(traceOneTimeslotCostNeverMigrate,2)./totalUsers),'--k',xPoints,ones(1,length(xPoints))*mean(sum(traceOneTimeslotCostAlwaysMigrate,2)./totalUsers),'--g',...
    xPoints,ones(1,length(xPoints))*mean(sum(traceOneTimeslotCostMyopic,2)./totalUsers),'--b',xPoints,ones(1,length(xPoints))*mean(sum(traceOneTimeslotCostThPolicy,2)./totalUsers),'--r')
legend('Never Migrate (A)','Always Migrate (B)','Myopic (C)','Proposed (D)','(Average values)')
xlabel('Time','FontSize',12)
ylabel('Instantaneous cost','FontSize',12)
datetick('x','HH:MM');

figure(2)
plot([xPoints(1), xPoints(end)] ,[1,1]*mean(sum(traceOneTimeslotCostNeverMigrate,2)./totalUsers),'--k',[xPoints(1), xPoints(end)],[1,1]*mean(sum(traceOneTimeslotCostAlwaysMigrate,2)./totalUsers),'--g',...
    [xPoints(1), xPoints(end)],[1,1]*mean(sum(traceOneTimeslotCostMyopic,2)./totalUsers),'--b',[xPoints(1), xPoints(end)],[1,1]*mean(sum(traceOneTimeslotCostThPolicy,2)./totalUsers),'--r')
datetick('x','HH:MM');











