clear all
close all

%% init all loops
timeIterativeAvrStore=[];
timeThAvrStore=[];
timePolicyAvrStore=[];
timeThPolicyAvrStore=[];
valueErrorAvrStore=[];
valuePolicyAvrStore=[];
valueValueResultAvrStore=[];
valueValueActualAvrStore=[];
valueNeverMigrateAvrStore=[];
valueAlwaysMigrateAvrStore=[];
valueMyopicAvrStore=[];
differentActionCountPerctStrore=[];
valuesThPolicyActualAvrStore=[];

gammaStore=[];
migrateProportionalStore=[];
%%

Use2D=1;


gammaVector=[0.5, 0.9, 0.99];
migrateProportionalVector=[0:0.25:1, 2, 4, 6, 10, 15, 20];  %weighting factor for migration cost
    
for gamma=gammaVector  %discount factor

    for migrateProportional=migrateProportionalVector

        %to match with the range in plots since value and policy iterations take long to run
        if (gamma==0.5) && (migrateProportional>2)
            notUseValueIteration=1;
            notUsePolicyIteration=1;
        elseif (gamma==0.9) && (migrateProportional>10)
            notUseValueIteration=1;
            notUsePolicyIteration=1;
        else
            notUseValueIteration=0;
            notUsePolicyIteration=0;
        end

        simSeedVector=1:50;   %specifies how many simulations to run, the average results are plotted

        %% init inner loop
        timeIterativeTotal=0;
        timeThTotal=0;
        timePolicyTotal=0;
        timeThPolicyTotal=0;
        valueErrorTotal=0;
        valuePolicyTotal=0;
        valueValueResultTotal=0;
        valueValueActualTotal=0;
        valueNeverMigrateTotal=0;
        valueAlwaysMigrateTotal=0;
        valueMyopicTotal=0;
        differentActionCountTotal=0;
        valuesThPolicyActualTotal=0;
        %%

        for simSeed=simSeedVector
            %RandStream.setGlobalStream(RandStream('mt19937ar','seed',simSeed));   %for fixed random seed, may not work on some MATLAB versions
            
            %% cost function parameters
            powerFactor=0.8;
            proportionalFactorMigrate=-migrateProportional;
            proportionalFactorTrans=-1;
            constFactorMigrate=1-proportionalFactorMigrate;
            constFactorTrans=0-proportionalFactorTrans;
            
            %% Config transition prob.
            if Use2D==0  %random walk parameters for 1D mobility
                %parameters for 1D mobility
                numStatesLeft=0;
                numStatesRight=10;
                
                p_forward=rand();
                p_back=rand()*(1-p_forward);
                p_same=1-p_forward-p_back;

                if numStatesLeft>0
                    p_out_state_first=0;
                else
                    p_out_state_first=rand();   %transition prob. from state zero for one-sided Markov chain
                end
                p_out_state_last=0;
                
                numStates=numStatesLeft+numStatesRight+1;
                zeroStateIndex=numStatesLeft+1;

            else   %random walk parameters for 2D mobility

                %for 2D mobility transition matrix (hexagon cell)
                
                %----- This part same as in dataset collection -----
                %parameters have to match the dataset
                numStates2D=10;   %excluding state 0
                cellDist=0.005;
                centerCoordinate=[37.762 -122.43];
                
                cellDistX=cellDist;
                cellDistY=cellDist*sqrt(3)/2;
                
                numStates2DTotal=1;
                startingIndexEachCircle=zeros(numStates2D+1,1);
                for i=1:numStates2D
                    startingIndexEachCircle(i)=numStates2DTotal+1;
                    numStates2DTotal=numStates2DTotal+i*6;
                end
                startingIndexEachCircle(numStates2D+1)=numStates2DTotal+1;

                coordinatesCells2D=zeros(numStates2DTotal,2);
                
                %find coordinates
                for i=1:numStates2D
                    for j=1:6
                        coordinatesCells2D(startingIndexEachCircle(i)+(j-1)*i,1)=cellDist*i*cos(2*pi/6*(j-1));
                        coordinatesCells2D(startingIndexEachCircle(i)+(j-1)*i,2)=cellDist*i*sin(2*pi/6*(j-1));     
                    end
                    for j=1:6
                        prevCornerPoint=startingIndexEachCircle(i)+(mod(j-1,6))*i;
                        nextCornerPoint=startingIndexEachCircle(i)+(mod(j,6))*i;
                        for k=prevCornerPoint+1:prevCornerPoint+i-1
                            coordinatesCells2D(k,:)=coordinatesCells2D(prevCornerPoint,:)+ ...
                                (coordinatesCells2D(nextCornerPoint,:)-coordinatesCells2D(prevCornerPoint,:))/i*(k-prevCornerPoint);
                        end
                    end
                end
                
                coordinatesCells2D(:,1)=coordinatesCells2D(:,1)+centerCoordinate(1);
                coordinatesCells2D(:,2)=coordinatesCells2D(:,2)+centerCoordinate(2);
                
                %find neighboring cell indices
                neighborCells2D=zeros(numStates2DTotal,7); %index i,j: i - given cell, j - neighbor (j=1 is itself)
                for i=1:numStates2DTotal
                    tmp=[];
                    for j=1:numStates2DTotal
                        if abs(norm(coordinatesCells2D(i,:)-coordinatesCells2D(j,:))-cellDist)<0.000000000001
                            tmp=[tmp, j];
                        end
                    end
                    
                    for j=tmp
                        delta=round((coordinatesCells2D(j,:)-coordinatesCells2D(i,:))*1e8)/1e8;
                        
                        neighborCells2D(i,1)=i;
                        
                        if (delta(1)>0 && delta(2)==0)
                            neighborCells2D(i,2)=j;
                        elseif (delta(1)>0 && delta(2)>0)
                            neighborCells2D(i,3)=j;
                        elseif (delta(1)<0 && delta(2)>0)
                            neighborCells2D(i,4)=j;
                        elseif (delta(1)<0 && delta(2)==0)
                            neighborCells2D(i,5)=j;
                        elseif (delta(1)<0 && delta(2)<0)
                            neighborCells2D(i,6)=j;
                        elseif (delta(1)>0 && delta(2)<0)
                            neighborCells2D(i,7)=j;
                        end                
                    end
                    
                end
                
                %-------------------------------------------------

                
                %*** for symmetric random walk ***
                p_2D=rand()*1/6;
                P2D=zeros(numStates2DTotal,numStates2DTotal);

                %assign for same circle, starting from circle 1
                for i=1:numStates2D
                    for j=startingIndexEachCircle(i):startingIndexEachCircle(i+1)-1
                        if j+1>startingIndexEachCircle(i+1)-1
                            tmp=startingIndexEachCircle(i);
                        else
                            tmp=j+1;
                        end
                        P2D(j,tmp)=p_2D;
                        P2D(tmp,j)=p_2D;            
                    end
                end

                %assign for connections with zero state
                for i=1:6
                    P2D(1,i+1)=p_2D;
                    P2D(i+1,1)=p_2D;
                end
                %assign connections with higher circle
                for i=1:numStates2D-1
                    for j=1:i*6
                        for k=ceil((i+1)/i*(j-1)) : ceil((i+1)/i*j)
                            if (k==0)
                                tmp=(i+1)*6;
                            else
                                tmp=k;
                            end
                            P2D(startingIndexEachCircle(i)-1+j,startingIndexEachCircle(i+1)-1+tmp)=p_2D;
                            P2D(startingIndexEachCircle(i+1)-1+tmp,startingIndexEachCircle(i)-1+j)=p_2D;
                        end
                    end        
                end    

                %assign for same state
                for i=1:numStates2DTotal
                    P2D(i,i)=1-sum(P2D(i,:));
                end
            
                %********************************
                

                
            end
                
            %% run algorithm
            disp(['** simulation seed: ', num2str(simSeed)])
            disp(['gamma: ', num2str(gamma)])
            disp(['migrateProportional: ', num2str(migrateProportional)])

            run([pwd, '\algorithms.m'])
            if algReturn==1
                return
            end
            
            %% Compute total values
            timeIterativeTotal=timeIterativeTotal+timeValue;
            timePolicyTotal=timePolicyTotal+timePolicy;
            timeThPolicyTotal=timeThPolicyTotal+timeThPolicy;    
            valueErrorTotal=valueErrorTotal+valueError;    
            valuePolicyTotal=valuePolicyTotal+sum(valuesPolicy);
            valueValueResultTotal=valueValueResultTotal+sum(valuesValueResult);
            valueValueActualTotal=valueValueActualTotal+sum(valuesValueActual);
            valueNeverMigrateTotal=valueNeverMigrateTotal+sum(valuesNeverMigrate);
            valueAlwaysMigrateTotal=valueAlwaysMigrateTotal+sum(valuesAlwaysMigrate);
            valueMyopicTotal=valueMyopicTotal+sum(valuesMyopic);
            if Use2D==1
                valuesThPolicyActualTotal=valuesThPolicyActualTotal+sum(valuesThPolicyActual);
            end
        end

        %% record results
        timeIterativeAvr=timeIterativeTotal/length(simSeedVector);
        timeThAvr=timeThTotal/length(simSeedVector);
        timePolicyAvr=timePolicyTotal/length(simSeedVector);
        timeThPolicyAvr=timeThPolicyTotal/length(simSeedVector);
        valueErrorAvr=valueErrorTotal/length(simSeedVector)/numStatesForStandardSolution;
        valuePolicyAvr=valuePolicyTotal/length(simSeedVector)/numStatesForStandardSolution;
        valueValueResultAvr=valueValueResultTotal/length(simSeedVector)/numStatesForStandardSolution;
        valueValueActualAvr=valueValueActualTotal/length(simSeedVector)/numStatesForStandardSolution;
        valueNeverMigrateAvr=valueNeverMigrateTotal/length(simSeedVector)/numStatesForStandardSolution;
        valueAlwaysMigrateAvr=valueAlwaysMigrateTotal/length(simSeedVector)/numStatesForStandardSolution;
        valueMyopicAvr=valueMyopicTotal/length(simSeedVector)/numStatesForStandardSolution;
        differentActionCountPerct=differentActionCountTotal/length(simSeedVector);
        if Use2D==1
            valuesThPolicyActualAvr=valuesThPolicyActualTotal/length(simSeedVector)/numStatesForStandardSolution;
        end

        timeIterativeAvrStore=[timeIterativeAvrStore,timeIterativeAvr];
        timeThAvrStore=[timeThAvrStore,timeThAvr];
        timePolicyAvrStore=[timePolicyAvrStore,timePolicyAvr];
        timeThPolicyAvrStore=[timeThPolicyAvrStore,timeThPolicyAvr];
        valueErrorAvrStore=[valueErrorAvrStore,valueErrorAvr];
        valuePolicyAvrStore=[valuePolicyAvrStore,valuePolicyAvr];
        valueValueResultAvrStore=[valueValueResultAvrStore,valueValueResultAvr];
        valueValueActualAvrStore=[valueValueActualAvrStore,valueValueActualAvr];
        valueNeverMigrateAvrStore=[valueNeverMigrateAvrStore,valueNeverMigrateAvr];
        valueAlwaysMigrateAvrStore=[valueAlwaysMigrateAvrStore,valueAlwaysMigrateAvr];
        valueMyopicAvrStore=[valueMyopicAvrStore,valueMyopicAvr];
        differentActionCountPerctStrore=[differentActionCountPerctStrore,differentActionCountPerct];
        if Use2D==1
            valuesThPolicyActualAvrStore=[valuesThPolicyActualAvrStore,valuesThPolicyActualAvr];
        end

        gammaStore=[gammaStore, gamma];
        migrateProportionalStore=[migrateProportionalStore, migrateProportional];

    end
end


%for plotting figure

sizeGamma=length(gammaVector);
sizeProportional=length(migrateProportionalVector);

timeThPolicyAvrMatrix=zeros(sizeGamma,sizeProportional);
timePolicyAvrMatrix=zeros(sizeGamma,sizeProportional);
timeIterativeAvrMatrix=zeros(sizeGamma,sizeProportional);

valuesThPolicyActualAvrMatrix=zeros(sizeGamma,sizeProportional);
valuePolicyAvrMatrix=zeros(sizeGamma,sizeProportional);
valueNeverMigrateAvrMatrix=zeros(sizeGamma,sizeProportional);
valueAlwaysMigrateAvrMatrix=zeros(sizeGamma,sizeProportional);
valueMyopicAvrMatrix=zeros(sizeGamma,sizeProportional);

gammaFactors=zeros(sizeGamma,1); 
proportionalFactors=zeros(sizeProportional,1); 

for g=1:sizeGamma
    for i=1:sizeProportional
        timeThPolicyAvrMatrix(g, i)=timeThPolicyAvrStore((g-1)*sizeProportional+i);
        timePolicyAvrMatrix(g, i)=timePolicyAvrStore((g-1)*sizeProportional+i);
        timeIterativeAvrMatrix(g, i)=timeIterativeAvrStore((g-1)*sizeProportional+i);

        valuesThPolicyActualAvrMatrix(g, i)=valuesThPolicyActualAvrStore((g-1)*sizeProportional+i);
        valuePolicyAvrMatrix(g, i)=valuePolicyAvrStore((g-1)*sizeProportional+i);
        valueNeverMigrateAvrMatrix(g, i)=valueNeverMigrateAvrStore((g-1)*sizeProportional+i);
        valueAlwaysMigrateAvrMatrix(g, i)=valueAlwaysMigrateAvrStore((g-1)*sizeProportional+i);
        valueMyopicAvrMatrix(g, i)=valueMyopicAvrStore((g-1)*sizeProportional+i);

        gammaFactors(g)=gammaStore((g-1)*sizeProportional+i);
        proportionalFactors(i)=migrateProportionalStore((g-1)*sizeProportional+i);

    end
end

simParamVector=proportionalFactors;

g=1;
i=1:sizeProportional;
j=1;


set(0,'DefaultTextFontname','Times New Roman', ...
'DefaultAxesFontname','Times New Roman')


figure(1)

subplot(2,1,1)
semilogy(simParamVector,timeThPolicyAvrMatrix(g,i,j),'-ko',simParamVector,timePolicyAvrMatrix(g,i,j),'-.bx',simParamVector,timeIterativeAvrMatrix(g,i,j),'--r^','LineWidth',2)
legend('Proposed','Policy iteration','Value iteration')
xlabel('-\it\beta_l','FontSize',12)
ylabel('Computation time (s)','FontSize',12)
xlim([0 2])

subplot(2,1,2)
plot(simParamVector,valuesThPolicyActualAvrMatrix(g,i,j),'-ko',simParamVector,valuePolicyAvrMatrix(g,i,j),':mx',simParamVector,valueNeverMigrateAvrMatrix(g,i,j),'--g^',simParamVector,valueAlwaysMigrateAvrMatrix(g,i,j),'-.bv',simParamVector,valueMyopicAvrMatrix(g,i,j),':rs','LineWidth',2)
legend('Proposed','Optimal','Never migrate','Always migrate','Myopic')
xlabel('-\it\beta_l','FontSize',12)
ylabel('Discounted sum cost','FontSize',12)
axis([0 2 0 2.5])

figure(2)

subplot(2,1,1)
semilogy(simParamVector,timeThPolicyAvrMatrix(2,:),'-ko',simParamVector,timePolicyAvrMatrix(2,:),'-.bx',simParamVector,timeIterativeAvrMatrix(2,:),'--r^','LineWidth',2)
legend('Proposed','Policy iteration','Value iteration')
xlabel('-\it\beta_l','FontSize',12)
ylabel('Computation time (s)','FontSize',12)
xlim([0 8])

subplot(2,1,2)
plot(simParamVector,valuesThPolicyActualAvrMatrix(2,:),'-ko',simParamVector,valuePolicyAvrMatrix(2,:),':mx',simParamVector,valueNeverMigrateAvrMatrix(2,:),'--g^',simParamVector,valueAlwaysMigrateAvrMatrix(2,:),'-.bv',simParamVector,valueMyopicAvrMatrix(2,:),':rs','LineWidth',2)
legend('Proposed','Optimal','Never migrate','Always migrate','Myopic')
xlabel('-\it\beta_l','FontSize',12)
ylabel('Discounted sum cost','FontSize',12)
axis([0 8 0 20])


figure(3)

subplot(2,1,1)
semilogy(simParamVector,timeThPolicyAvrMatrix(3,:),'-ko',simParamVector,timePolicyAvrMatrix(3,:),'-.bx',simParamVector,timeIterativeAvrMatrix(3,:),'--r^','LineWidth',2)
legend('Proposed','Policy iteration','Value iteration')
xlabel('-\it\beta_l','FontSize',12)
ylabel('Computation time (s)','FontSize',12)

subplot(2,1,2)
plot(simParamVector,valuesThPolicyActualAvrMatrix(3,:),'-ko',simParamVector,valuePolicyAvrMatrix(3,:),':mx',simParamVector,valueNeverMigrateAvrMatrix(3,:),'--g^',simParamVector,valueAlwaysMigrateAvrMatrix(3,:),'-.bv',simParamVector,valueMyopicAvrMatrix(3,:),':rs','LineWidth',2)
legend('Proposed','Optimal','Never migrate','Always migrate','Myopic')
xlabel('-\it\beta_l','FontSize',12)
ylabel('Discounted sum cost','FontSize',12)



