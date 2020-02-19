%% init for algorithms

epsilonValueIteration=0.1;
numValueIteration=5;

if ~exist('notUseValueIteration')
    notUseValueIteration=0;
end
if ~exist('notUsePolicyIteration')
    notUsePolicyIteration=0;
end
if ~exist('supressOutput')
    supressOutput=0;
end    
if ~exist('notFindActualValues')
    notFindActualValues=0;
end        

algReturn=0;

if Use2D==1
    %Set parameters for 1D approximate Markov chain
    p_forward=(1-P2D(1,1))/6*2.5;
    p_back=(1-P2D(1,1))/6*1.5;
    p_same=1-p_forward-p_back;
    p_out_state_first=1-P2D(1,1);
    p_out_state_last=(1-P2D(1,1))/6;           

    numStates=numStates2D+1;
    zeroStateIndex=1;
end


% Find 1D (1D approximation if using 2D) transition prob. matrix from random walk parameters
P=zeros(numStates,numStates);
for i=1:numStates
    for j=1:numStates
        if i==j
            if i==1
                P(i,j)=1-p_out_state_first;
            elseif i==numStates
                P(i,j)=1-p_out_state_last;
            else
                P(i,j)=p_same;
            end
        elseif j==i+1
            if i==1
                P(i,j)=p_out_state_first;
            else
                P(i,j)=p_forward;
            end
        elseif j==i-1
            if i==numStates
                P(i,j)=p_out_state_last;
            else
                P(i,j)=p_back;
            end
        end
    end

end

%set transition prob. matrix for solution with standard policies
if Use2D==0
    PForStandardSolution=P;
    numStatesForStandardSolution=numStates;
    zeroStateIndexForStandardSolution=zeroStateIndex;
else
    PForStandardSolution=P2D;
    numStatesForStandardSolution=numStates2DTotal;
    zeroStateIndexForStandardSolution=1;
end


%% standard solution method -- value iteration
tic
valuesValueResult=zeros(numStatesForStandardSolution,1); %valuesValueResult when at corresponding mobile state
actionsValue=zeroStateIndexForStandardSolution*ones(numStatesForStandardSolution,1);  %action to take when at i=mobile state
if notUseValueIteration==0
    for timeslot=1:numValueIteration %assume increasing timeslot (reverse indexing)
        valuesPrev=valuesValueResult;

        for S_mobile=1:numStatesForStandardSolution
            %iterative policy
            tmpMinValue=inf;

            if Use2D==0
                if S_mobile==1
                    if zeroStateIndexForStandardSolution>1
                        possibleActionVector=2:numStatesForStandardSolution;
                    else
                        possibleActionVector=1:numStatesForStandardSolution;
                    end
                elseif S_mobile==numStatesForStandardSolution
                    possibleActionVector=1:numStatesForStandardSolution-1;
                else
                    possibleActionVector=1:numStatesForStandardSolution;
                end
            else 
                if S_mobile>=startingIndexEachCircle(numStates2D)
                    possibleActionVector=1:startingIndexEachCircle(numStates2D)-1;
                else
                    possibleActionVector=1:numStatesForStandardSolution;
                end
            end

            for possibleAction=possibleActionVector
                if Use2D==0
                    tmp=(constFactorMigrate+proportionalFactorMigrate*powerFactor^abs(possibleAction-S_mobile))*(S_mobile~=possibleAction)...
                        +(constFactorTrans+proportionalFactorTrans*powerFactor^abs(zeroStateIndexForStandardSolution-possibleAction))*(zeroStateIndexForStandardSolution~=possibleAction)...
                        +gamma*PForStandardSolution(possibleAction,:)*valuesValueResult;
                else
                    tmpNorm1=norm(coordinatesCells2D(possibleAction,:)-coordinatesCells2D(S_mobile,:))/cellDist;
                    if abs(tmpNorm1-round(tmpNorm1))<0.00000001
                        tmpNorm1=round(tmpNorm1);
                    else
                        tmpNorm1=ceil(tmpNorm1);
                    end

                    tmpNorm2=norm(coordinatesCells2D(zeroStateIndexForStandardSolution,:)-coordinatesCells2D(possibleAction,:))/cellDist;
                    if abs(tmpNorm2-round(tmpNorm2))<0.00000001
                        tmpNorm2=round(tmpNorm2);
                    else
                        tmpNorm2=ceil(tmpNorm2);
                    end

                    tmp=(constFactorMigrate+proportionalFactorMigrate*powerFactor^tmpNorm1)*(S_mobile~=possibleAction)...
                        +(constFactorTrans+proportionalFactorTrans*powerFactor^tmpNorm2)*(zeroStateIndexForStandardSolution~=possibleAction)...
                        +gamma*PForStandardSolution(possibleAction,:)*valuesValueResult;
                end

                if tmp<tmpMinValue
                    tmpMinValue=tmp;
                    actionsValue(S_mobile)=possibleAction;
                end

            end            

            valuesValueResult(S_mobile)=tmpMinValue;


        end

        if abs(valuesPrev-valuesValueResult)<epsilonValueIteration*(1-gamma)/(2*gamma)
            break
        end
    end
end

timeValue=toc;


%% standard solution method -- policy iteration
tic
valuesPolicy=zeros(numStatesForStandardSolution,1); %valuesValueResult when at corresponding mobile state
actionsPolicy=zeroStateIndexForStandardSolution*ones(numStatesForStandardSolution,1);  %action to take when at i=mobile state

if notUsePolicyIteration==0

    tmpCountPolicy=0;

    while(1)
        tmpCountPolicy=tmpCountPolicy+1;
        actionsPolicyPrev=actionsPolicy;
        Pmodified=zeros(numStatesForStandardSolution,numStatesForStandardSolution);
        ck=zeros(numStatesForStandardSolution,1);

        for i=1:length(actionsPolicy)
            Pmodified(i,:)=PForStandardSolution(actionsPolicy(i),:);
            if Use2D==0
                ck(i)=(constFactorMigrate+proportionalFactorMigrate*powerFactor^abs(actionsPolicy(i)-i))*(actionsPolicy(i)~=i)... 
                    + (constFactorTrans+proportionalFactorTrans*powerFactor^abs(actionsPolicy(i)-zeroStateIndexForStandardSolution))*(actionsPolicy(i)~=zeroStateIndexForStandardSolution);
            else
                tmpNorm1=norm(coordinatesCells2D(actionsPolicy(i),:)-coordinatesCells2D(i,:))/cellDist;
                if abs(tmpNorm1-round(tmpNorm1))<0.00000001
                    tmpNorm1=round(tmpNorm1);
                else
                    tmpNorm1=ceil(tmpNorm1);
                end

                tmpNorm2=norm(coordinatesCells2D(actionsPolicy(i),:)-coordinatesCells2D(zeroStateIndexForStandardSolution,:))/cellDist;
                if abs(tmpNorm2-round(tmpNorm2))<0.00000001
                    tmpNorm2=round(tmpNorm2);
                else
                    tmpNorm2=ceil(tmpNorm2);
                end

                ck(i)=(constFactorMigrate+proportionalFactorMigrate*powerFactor^tmpNorm1)*(actionsPolicy(i)~=i)...
                    + (constFactorTrans+proportionalFactorTrans*powerFactor^tmpNorm2)*(actionsPolicy(i)~=zeroStateIndexForStandardSolution);
            end
        end

        valuesPolicy=(eye(length(Pmodified))-gamma*Pmodified)\ck;

        for S_mobile=1:numStatesForStandardSolution
            tmpMinValue=inf;

            if Use2D==0
                if S_mobile==1
                    if zeroStateIndexForStandardSolution>1
                        possibleActionVector=2:numStatesForStandardSolution;
                    else
                        possibleActionVector=1:numStatesForStandardSolution;
                    end
                elseif S_mobile==numStatesForStandardSolution
                    possibleActionVector=1:numStatesForStandardSolution-1;
                else
                    possibleActionVector=1:numStatesForStandardSolution;
                end
            else 
                if S_mobile>=startingIndexEachCircle(numStates2D)
                    possibleActionVector=1:startingIndexEachCircle(numStates2D)-1;
                else
                    possibleActionVector=1:numStatesForStandardSolution;
                end
            end

            for possibleAction=possibleActionVector
                if Use2D==0
                    tmp=(constFactorMigrate+proportionalFactorMigrate*powerFactor^abs(possibleAction-S_mobile))*(S_mobile~=possibleAction)...
                        +(constFactorTrans+proportionalFactorTrans*powerFactor^abs(zeroStateIndexForStandardSolution-possibleAction))*(zeroStateIndexForStandardSolution~=possibleAction)...
                        +gamma*PForStandardSolution(possibleAction,:)*valuesPolicy;
                else
                    tmpNorm1=norm(coordinatesCells2D(possibleAction,:)-coordinatesCells2D(S_mobile,:))/cellDist;
                    if abs(tmpNorm1-round(tmpNorm1))<0.00000001
                        tmpNorm1=round(tmpNorm1);
                    else
                        tmpNorm1=ceil(tmpNorm1);
                    end

                    tmpNorm2=norm(coordinatesCells2D(zeroStateIndexForStandardSolution,:)-coordinatesCells2D(possibleAction,:))/cellDist;
                    if abs(tmpNorm2-round(tmpNorm2))<0.00000001
                        tmpNorm2=round(tmpNorm2);
                    else
                        tmpNorm2=ceil(tmpNorm2);
                    end

                    tmp=(constFactorMigrate+proportionalFactorMigrate*powerFactor^tmpNorm1)*(S_mobile~=possibleAction)...
                        +(constFactorTrans+proportionalFactorTrans*powerFactor^tmpNorm2)*(zeroStateIndexForStandardSolution~=possibleAction)...
                        +gamma*PForStandardSolution(possibleAction,:)*valuesPolicy;
                end

                if tmp<tmpMinValue
                    tmpMinValue=tmp;
                    actionsPolicy(S_mobile)=possibleAction;
                end

            end

        end

        if actionsPolicy==actionsPolicyPrev
            break
        end
    end
end

timePolicy=toc;


%% modified policy iteration with difference equations
tic


valuesThPolicy=zeros(numStates,1); %valuesValueResult when at corresponding mobile state
actionsThPolicy=zeroStateIndex*ones(numStates,1);  %action to take when at i=mobile state

if zeroStateIndex==1  
    %% One direction
    %------ Constant variables for difference equation -------
    alpha1=gamma*p_back/(1-gamma*(1-p_forward-p_back));
    alpha2=gamma*p_forward/(1-gamma*(1-p_forward-p_back));
    alpha3=constFactorTrans/(1-gamma*(1-p_forward-p_back));
    alpha4=proportionalFactorTrans/(1-gamma*(1-p_forward-p_back));
    
    alpha0=(gamma*p_out_state_first)/(1-gamma*(1-p_out_state_first));

    m1=(1+sqrt(1-4*alpha1*alpha2))/(2*alpha2);
    m2=(1-sqrt(1-4*alpha1*alpha2))/(2*alpha2);

    D=alpha3/(1-alpha1-alpha2);
    B=alpha4/(1-alpha1/powerFactor-alpha2*powerFactor);
    %------ End Constant variables for difference equation -------

    tmpCountProposed=0;
    while(1)
        tmpCountProposed=tmpCountProposed+1;
        actionsThPolicyPrev=actionsThPolicy;


        tmpHighestIndex=find(actionsThPolicy<[1:numStates]', 1, 'first' );
        rightMostStateMapTo=actionsThPolicy(tmpHighestIndex);

        %-------- Solve difference equation ----------
        %mapping-specific parameters
        abs_N=tmpHighestIndex-1;
        abs_aN=(rightMostStateMapTo-1);

        A=[1-alpha0*m1, 1-alpha0*m2; ...
            m1^(tmpHighestIndex-1)-m1^(rightMostStateMapTo-1), m2^(tmpHighestIndex-1)-m2^(rightMostStateMapTo-1)];

        const=[D*(alpha0-1)+B*(alpha0*powerFactor-1); ...
            constFactorMigrate+proportionalFactorMigrate*powerFactor^(tmpHighestIndex-rightMostStateMapTo)-B*(powerFactor^(tmpHighestIndex-1)-powerFactor^abs(rightMostStateMapTo-1))];

        coeffDet=A(1,1)*A(2,2)-A(1,2)*A(2,1);

        A1=(const(1)*A(2,2)-A(1,2)*const(2))/coeffDet;
        A2=(const(2)*A(1,1)-A(2,1)*const(1))/coeffDet;            
        
        %find all costs
        for i=1:tmpHighestIndex
            valuesThPolicy(i)=A1*m1^(i-1)+A2*m2^(i-1)+D+B*powerFactor^(i-1);
        end     
        %---------- End Solve difference equation ---------


        while(1)
            tmpHighestIndexPrev=tmpHighestIndex;
            tmpHighestIndex=find(actionsThPolicy(tmpHighestIndexPrev+1:numStates)<[tmpHighestIndexPrev+1:numStates]', 1, 'first' )+(tmpHighestIndexPrev);
            if isempty(tmpHighestIndex)
                break
            end

            rightMostStateMapTo=actionsThPolicy(tmpHighestIndex);

            if rightMostStateMapTo<=tmpHighestIndexPrev
                A=[m1^(tmpHighestIndexPrev-1), m2^(tmpHighestIndexPrev-1); ...
                    m1^(tmpHighestIndex-1), m2^(tmpHighestIndex-1)];

                const=[valuesThPolicy(tmpHighestIndexPrev)-D-B*powerFactor^abs(tmpHighestIndexPrev-1); ...
                    constFactorMigrate+proportionalFactorMigrate*powerFactor^(tmpHighestIndex-rightMostStateMapTo)+valuesThPolicy(rightMostStateMapTo)-D-B*powerFactor^(tmpHighestIndex-1)];
            else
                A=[m1^(tmpHighestIndexPrev-1), m2^(tmpHighestIndexPrev-1); ...
                    m1^(tmpHighestIndex-1)-m1^abs(rightMostStateMapTo-1), ...
                    m2^(tmpHighestIndex-1)-m2^abs(rightMostStateMapTo-1)];

                const=[valuesThPolicy(tmpHighestIndexPrev)-D-B*powerFactor^(tmpHighestIndexPrev-1); ...
                    constFactorMigrate+proportionalFactorMigrate*powerFactor^(tmpHighestIndex-rightMostStateMapTo)-B*(powerFactor^(tmpHighestIndex-1)-powerFactor^abs(rightMostStateMapTo-1))];
            end

            coeffDet=A(1,1)*A(2,2)-A(1,2)*A(2,1);

            A1=(const(1)*A(2,2)-A(1,2)*const(2))/coeffDet;
            A2=(const(2)*A(1,1)-A(2,1)*const(1))/coeffDet;

            for i=tmpHighestIndexPrev+1:tmpHighestIndex
                valuesThPolicy(i)=A1*m1^abs(i-1)+A2*m2^abs(i-1)+D+B*powerFactor^abs(i-1);
            end

        end

        for S_mobile=1:numStates
            tmpMinValue=inf;

            if S_mobile==1
                if zeroStateIndexForStandardSolution>1
                    possibleActionVector=2:S_mobile;
                else
                    possibleActionVector=1:S_mobile;
                end
            elseif S_mobile==numStates
                possibleActionVector=1:S_mobile-1;
            else
                possibleActionVector=1:S_mobile;
            end

            for possibleAction=possibleActionVector

                tmp=(constFactorMigrate+proportionalFactorMigrate*powerFactor^abs(possibleAction-S_mobile))*(S_mobile~=possibleAction)...
                    +(constFactorTrans+proportionalFactorTrans*powerFactor^abs(zeroStateIndex-possibleAction))*(zeroStateIndex~=possibleAction)...
                    +gamma*P(possibleAction,:)*valuesThPolicy;

                if tmp<tmpMinValue
                    tmpMinValue=tmp;
                    actionsThPolicy(S_mobile)=possibleAction;
                end
            end
        end

        if actionsThPolicy==actionsThPolicyPrev
            break
        end
    end            
    
    
else
    %% Two directions (there are some problems for this case, probably due to numerical inaccuracy)
    %------ Constant variables for difference equation -------
    alpha1=gamma*p_back/(1-gamma*(1-p_forward-p_back));
    alpha2=gamma*p_forward/(1-gamma*(1-p_forward-p_back));
    alpha3=constFactorTrans/(1-gamma*(1-p_forward-p_back));
    alpha4=proportionalFactorTrans/(1-gamma*(1-p_forward-p_back));

    m1=(1+sqrt(1-4*alpha1*alpha2))/(2*alpha2);
    m2=(1-sqrt(1-4*alpha1*alpha2))/(2*alpha2);
    m1r=(1+sqrt(1-4*alpha1*alpha2))/(2*alpha1);
    m2r=(1-sqrt(1-4*alpha1*alpha2))/(2*alpha1);

    D=alpha3/(1-alpha1-alpha2);
    B=alpha4/(1-alpha1/powerFactor-alpha2*powerFactor);
    Br=alpha4/(1-alpha2/powerFactor-alpha1*powerFactor);

    H=(alpha1*(D*(1-m2r)+Br*(powerFactor-m2r))+alpha2*(D*(1-m2)+B*(powerFactor-m2)))/(sqrt(1-4*alpha1*alpha2));
    %------ End Constant variables for difference equation -------
    tmpCountProposed=0;
    while(1)
        tmpCountProposed=tmpCountProposed+1;
        actionsThPolicyPrev=actionsThPolicy;


        tmpLowestIndex=find(actionsThPolicy>[1:numStates]', 1, 'last' );
        tmpHighestIndex=find(actionsThPolicy<[1:numStates]', 1, 'first' );
        leftMostStateMapTo=actionsThPolicy(tmpLowestIndex);
        rightMostStateMapTo=actionsThPolicy(tmpHighestIndex);

        %-------- Solve difference equation ----------
        %mapping-specific parameters
        abs_M=zeroStateIndex-tmpLowestIndex;
        abs_aM=abs(leftMostStateMapTo-zeroStateIndex);
        if leftMostStateMapTo<=zeroStateIndex
            V0CoeffLeft=m1r^abs_M-m1r^abs_aM;
            A1CoeffLeft=m2r^abs_M-m1r^abs_M+m1r^abs_aM-m2r^abs_aM;
            constTermLeft=constFactorMigrate+proportionalFactorMigrate*powerFactor^(leftMostStateMapTo-tmpLowestIndex)+Br*(powerFactor^abs_aM-powerFactor^abs_M)+(D+Br-H)*(m2r^abs_M-m2r^abs_aM)+H*(m1r^abs_M-m1r^abs_aM);
        else
            V0CoeffLeft=m1r^abs_M-m2^abs_aM;
            A1CoeffLeft=m2r^abs_M-m1r^abs_M+m2^abs_aM-m1^abs_aM;
            constTermLeft=constFactorMigrate+proportionalFactorMigrate*powerFactor^(leftMostStateMapTo-tmpLowestIndex)+B*powerFactor^abs_aM-Br*powerFactor^abs_M+(D+Br-H)*(m2r^abs_M)+H*m1r^abs_M-(D+B)*m2^abs_aM;
        end

        abs_N=tmpHighestIndex-zeroStateIndex;
        abs_aN=abs(rightMostStateMapTo-zeroStateIndex);
        if rightMostStateMapTo<=zeroStateIndex
            V0CoeffRight=m2^abs_N-m1r^abs_aN;
            A1CoeffRight=m1^abs_N-m2^abs_N+m1r^abs_aN-m2r^abs_aN;
            constTermRight=constFactorMigrate+proportionalFactorMigrate*powerFactor^(tmpHighestIndex-rightMostStateMapTo)+Br*powerFactor^abs_aN-B*powerFactor^abs_N+(H-D-Br)*m2r^abs_aN-H*m1r^abs_aN+(D+B)*m2^abs_N;
        else
            V0CoeffRight=m2^abs_N-m2^abs_aN;
            A1CoeffRight=m1^abs_N-m2^abs_N+m2^abs_aN-m1^abs_aN;
            constTermRight=constFactorMigrate+proportionalFactorMigrate*powerFactor^(tmpHighestIndex-rightMostStateMapTo)+B*(powerFactor^abs_aN-powerFactor^abs_N)+(D+B)*(m2^abs_N-m2^abs_aN);
        end

        coeffDet=V0CoeffLeft*A1CoeffRight-A1CoeffLeft*V0CoeffRight;
        V0=(constTermLeft*A1CoeffRight-A1CoeffLeft*constTermRight)/coeffDet;

        A1=(V0CoeffLeft*constTermRight-constTermLeft*V0CoeffRight)/coeffDet;
        A2=V0-A1-D-B;
        A1r=V0-A1-H;
        A2r=V0-A1r-D-Br;

        %find all costs
        for i=tmpLowestIndex:tmpHighestIndex
            if (i<=zeroStateIndex)
                valuesThPolicy(i)=A1r*m1r^abs(i-zeroStateIndex)+A2r*m2r^abs(i-zeroStateIndex)+D+Br*powerFactor^abs(i-zeroStateIndex);
            elseif (i>zeroStateIndex)
                valuesThPolicy(i)=A1*m1^(i-zeroStateIndex)+A2*m2^(i-zeroStateIndex)+D+B*powerFactor^(i-zeroStateIndex);
            end
        end     
        %---------- End Solve difference equation ---------

        while(1)
            tmpLowestIndexPrev=tmpLowestIndex;
            tmpLowestIndex=find(actionsThPolicy(1:tmpLowestIndexPrev-1)>[1:tmpLowestIndexPrev-1]', 1, 'last' );
            if isempty(tmpLowestIndex)
                break
            end

            leftMostStateMapTo=actionsThPolicy(tmpLowestIndex);

            if leftMostStateMapTo>=tmpLowestIndexPrev
                A=[m1r^abs(tmpLowestIndexPrev-zeroStateIndex), m2r^abs(tmpLowestIndexPrev-zeroStateIndex); ...
                    m1r^abs(tmpLowestIndex-zeroStateIndex), m2r^abs(tmpLowestIndex-zeroStateIndex)];

                const=[valuesThPolicy(tmpLowestIndexPrev)-D-Br*powerFactor^abs(tmpLowestIndexPrev-zeroStateIndex); ...
                    constFactorMigrate+proportionalFactorMigrate*powerFactor^(leftMostStateMapTo-tmpLowestIndex)+valuesThPolicy(leftMostStateMapTo)-D-Br*powerFactor^abs(tmpLowestIndex-zeroStateIndex)];
            else
                A=[m1r^abs(tmpLowestIndexPrev-zeroStateIndex), m2r^abs(tmpLowestIndexPrev-zeroStateIndex); ...
                    m1r^abs(tmpLowestIndex-zeroStateIndex)-m1r^abs(leftMostStateMapTo-zeroStateIndex), ...
                    m2r^abs(tmpLowestIndex-zeroStateIndex)-m2r^abs(leftMostStateMapTo-zeroStateIndex)];

                const=[valuesThPolicy(tmpLowestIndexPrev)-D-Br*powerFactor^abs(tmpLowestIndexPrev-zeroStateIndex); ...
                    constFactorMigrate+proportionalFactorMigrate*powerFactor^(leftMostStateMapTo-tmpLowestIndex)-Br*(powerFactor^abs(tmpLowestIndex-zeroStateIndex)-powerFactor^abs(leftMostStateMapTo-zeroStateIndex))];
            end

            coeffDet=A(1,1)*A(2,2)-A(1,2)*A(2,1);

            A1r=(const(1)*A(2,2)-A(1,2)*const(2))/coeffDet;
            A2r=(const(2)*A(1,1)-A(2,1)*const(1))/coeffDet;

            for i=tmpLowestIndex:tmpLowestIndexPrev-1
                valuesThPolicy(i)=A1r*m1r^abs(i-zeroStateIndex)+A2r*m2r^abs(i-zeroStateIndex)+D+Br*powerFactor^abs(i-zeroStateIndex);
            end

        end

        while(1)
            tmpHighestIndexPrev=tmpHighestIndex;
            tmpHighestIndex=find(actionsThPolicy(tmpHighestIndexPrev+1:numStates)<[tmpHighestIndexPrev+1:numStates]', 1, 'first' )+(tmpHighestIndexPrev);
            if isempty(tmpHighestIndex)
                break
            end

            rightMostStateMapTo=actionsThPolicy(tmpHighestIndex);

            if rightMostStateMapTo<=tmpHighestIndexPrev
                A=[m1^(tmpHighestIndexPrev-zeroStateIndex), m2^(tmpHighestIndexPrev-zeroStateIndex); ...
                    m1^(tmpHighestIndex-zeroStateIndex), m2^(tmpHighestIndex-zeroStateIndex)];

                const=[valuesThPolicy(tmpHighestIndexPrev)-D-B*powerFactor^abs(tmpHighestIndexPrev-zeroStateIndex); ...
                    constFactorMigrate+proportionalFactorMigrate*powerFactor^(tmpHighestIndex-rightMostStateMapTo)+valuesThPolicy(rightMostStateMapTo)-D-B*powerFactor^(tmpHighestIndex-zeroStateIndex)];
            else
                A=[m1^(tmpHighestIndexPrev-zeroStateIndex), m2^(tmpHighestIndexPrev-zeroStateIndex); ...
                    m1^(tmpHighestIndex-zeroStateIndex)-m1^abs(rightMostStateMapTo-zeroStateIndex), ...
                    m2^(tmpHighestIndex-zeroStateIndex)-m2^abs(rightMostStateMapTo-zeroStateIndex)];

                const=[valuesThPolicy(tmpHighestIndexPrev)-D-B*powerFactor^(tmpHighestIndexPrev-zeroStateIndex); ...
                    constFactorMigrate+proportionalFactorMigrate*powerFactor^(tmpHighestIndex-rightMostStateMapTo)-B*(powerFactor^(tmpHighestIndex-zeroStateIndex)-powerFactor^abs(rightMostStateMapTo-zeroStateIndex))];
            end

            coeffDet=A(1,1)*A(2,2)-A(1,2)*A(2,1);

            A1=(const(1)*A(2,2)-A(1,2)*const(2))/coeffDet;
            A2=(const(2)*A(1,1)-A(2,1)*const(1))/coeffDet;

            for i=tmpHighestIndexPrev+1:tmpHighestIndex
                valuesThPolicy(i)=A1*m1^abs(i-zeroStateIndex)+A2*m2^abs(i-zeroStateIndex)+D+B*powerFactor^abs(i-zeroStateIndex);
            end

        end

        %--- for comparison with inverse method ---
        Pmodified=zeros(numStates,numStates);
        ck=zeros(numStates,1);

        for i=1:length(actionsThPolicy)
            Pmodified(i,:)=P(actionsThPolicy(i),:);
            ck(i)=(constFactorMigrate+proportionalFactorMigrate*powerFactor^abs(actionsThPolicy(i)-i))*(actionsThPolicy(i)~=i)... 
                + (constFactorTrans+proportionalFactorTrans*powerFactor^abs(actionsThPolicy(i)-zeroStateIndex))*(actionsThPolicy(i)~=zeroStateIndex);
        end

        valuesThPolicyFromInv=(eye(length(Pmodified))-gamma*Pmodified)\ck;
        %------------------------------------------


        for S_mobile=1:numStates
            tmpMinValue=inf;

            if Use2D==0
                if S_mobile==1
                    if zeroStateIndexForStandardSolution>1
                        possibleActionVector=2:numStatesForStandardSolution;
                    else
                        possibleActionVector=1:numStatesForStandardSolution;
                    end
                elseif S_mobile==numStates
                    possibleActionVector=1:numStates-1;
                else
                    possibleActionVector=1:numStates;
                end
            else 
                possibleActionVector=1:numStates;
            end

            for possibleAction=possibleActionVector

                tmp=(constFactorMigrate+proportionalFactorMigrate*powerFactor^abs(possibleAction-S_mobile))*(S_mobile~=possibleAction)...
                    +(constFactorTrans+proportionalFactorTrans*powerFactor^abs(zeroStateIndex-possibleAction))*(zeroStateIndex~=possibleAction)...
                    +gamma*P(possibleAction,:)*valuesThPolicy;

                if tmp<tmpMinValue
                    tmpMinValue=tmp;
                    actionsThPolicy(S_mobile)=possibleAction;
                end

            end

        end

        if actionsThPolicy==actionsThPolicyPrev
            break
        end
    end    
end


timeThPolicy=toc;    

%% myopic policy
Pmodified=zeros(numStatesForStandardSolution,numStatesForStandardSolution);
ck=zeros(numStatesForStandardSolution,1);

actionsMyopic=zeros(1,numStatesForStandardSolution);



for i=1:numStatesForStandardSolution
    if Use2D==0
        tmp=abs(i-zeroStateIndexForStandardSolution);
    else
        tmpNorm1=norm(coordinatesCells2D(i,:)-coordinatesCells2D(zeroStateIndexForStandardSolution,:))/cellDist;
        if abs(tmpNorm1-round(tmpNorm1))<0.00000001
            tmpNorm1=round(tmpNorm1);
        else
            tmpNorm1=ceil(tmpNorm1);
        end
        tmp=tmpNorm1;
    end
    if tmp==0
        actionsMyopic(i)=i;
    else
        if constFactorMigrate+proportionalFactorMigrate*powerFactor^tmp < constFactorTrans+proportionalFactorTrans*powerFactor^tmp
            actionsMyopic(i)=zeroStateIndexForStandardSolution;
        else
            actionsMyopic(i)=i;
        end
    end
end

%always migrate at border
if Use2D==0
    actionsMyopic(1)=zeroStateIndexForStandardSolution;
    actionsMyopic(numStatesForStandardSolution)=zeroStateIndexForStandardSolution;
else
    if exist('startingIndexEachCircle')
        actionsMyopic(startingIndexEachCircle(numStates2D):end)=ones(1,numStatesForStandardSolution-startingIndexEachCircle(numStates2D)+1);
    end
end



%% Other checks, outputs, and cost computations
if supressOutput==0
    if notUseValueIteration==0
        disp(['actions from value iteration: ', num2str(actionsValue')])
    end
    if notUsePolicyIteration==0
        disp(['actions from policy iteration: ', num2str(actionsPolicy')])
    end
    disp(['actions from modified policy iteration: ', num2str(actionsThPolicy')])
end    

if Use2D==0
    if length(find(actionsThPolicy~=actionsPolicy))>0
        disp('result not same')
        algReturn=1;
        return
    end    
end

%% find actual values
if notFindActualValues==0
    %% find actual value from value iteration
    Pmodified=zeros(numStatesForStandardSolution,numStatesForStandardSolution);
    ck=zeros(numStatesForStandardSolution,1);
    for i=1:numStatesForStandardSolution
        Pmodified(i,:)=PForStandardSolution(actionsValue(i),:);

        if Use2D==0
            ck(i)=(constFactorMigrate+proportionalFactorMigrate*powerFactor^abs(actionsValue(i)-i))*(actionsValue(i)~=i)... 
                + (constFactorTrans+proportionalFactorTrans*powerFactor^abs(actionsValue(i)-zeroStateIndexForStandardSolution))*(actionsValue(i)~=zeroStateIndexForStandardSolution);
        else
            tmpNorm1=norm(coordinatesCells2D(actionsValue(i),:)-coordinatesCells2D(i,:))/cellDist;
            if abs(tmpNorm1-round(tmpNorm1))<0.00000001
                tmpNorm1=round(tmpNorm1);
            else
                tmpNorm1=ceil(tmpNorm1);
            end

            tmpNorm2=norm(coordinatesCells2D(actionsValue(i),:)-coordinatesCells2D(zeroStateIndexForStandardSolution,:))/cellDist;
            if abs(tmpNorm2-round(tmpNorm2))<0.00000001
                tmpNorm2=round(tmpNorm2);
            else
                tmpNorm2=ceil(tmpNorm2);
            end

            ck(i)=(constFactorMigrate+proportionalFactorMigrate*powerFactor^tmpNorm1)*(actionsValue(i)~=i)...
                + (constFactorTrans+proportionalFactorTrans*powerFactor^tmpNorm2)*(actionsValue(i)~=zeroStateIndexForStandardSolution);
        end

    end
    valuesValueActual=(eye(length(Pmodified))-gamma*Pmodified)\ck;


    %% find value for never migrate mechanism
    Pmodified=zeros(numStatesForStandardSolution,numStatesForStandardSolution);
    ck=zeros(numStatesForStandardSolution,1);
    if Use2D==0
        actionsNeverMigrate=[zeroStateIndexForStandardSolution, 2:numStatesForStandardSolution-1, zeroStateIndexForStandardSolution];   %Always migrate to zero state at border
    else
        actionsNeverMigrate=[1:startingIndexEachCircle(numStates2D)-1,ones(1,numStatesForStandardSolution-startingIndexEachCircle(numStates2D)+1)];
    end

    for i=1:numStatesForStandardSolution
        Pmodified(i,:)=PForStandardSolution(actionsNeverMigrate(i),:);

        if Use2D==0
            ck(i)=(constFactorMigrate+proportionalFactorMigrate*powerFactor^abs(actionsNeverMigrate(i)-i))*(actionsNeverMigrate(i)~=i)... 
                + (constFactorTrans+proportionalFactorTrans*powerFactor^abs(actionsNeverMigrate(i)-zeroStateIndexForStandardSolution))*(actionsNeverMigrate(i)~=zeroStateIndexForStandardSolution);
        else
            tmpNorm1=norm(coordinatesCells2D(actionsNeverMigrate(i),:)-coordinatesCells2D(i,:))/cellDist;
            if abs(tmpNorm1-round(tmpNorm1))<0.00000001
                tmpNorm1=round(tmpNorm1);
            else
                tmpNorm1=ceil(tmpNorm1);
            end

            tmpNorm2=norm(coordinatesCells2D(actionsNeverMigrate(i),:)-coordinatesCells2D(zeroStateIndexForStandardSolution,:))/cellDist;
            if abs(tmpNorm2-round(tmpNorm2))<0.00000001
                tmpNorm2=round(tmpNorm2);
            else
                tmpNorm2=ceil(tmpNorm2);
            end

            ck(i)=(constFactorMigrate+proportionalFactorMigrate*powerFactor^tmpNorm1)*(actionsNeverMigrate(i)~=i)...
                + (constFactorTrans+proportionalFactorTrans*powerFactor^tmpNorm2)*(actionsNeverMigrate(i)~=zeroStateIndexForStandardSolution);
        end

    end
    valuesNeverMigrate=(eye(length(Pmodified))-gamma*Pmodified)\ck;

    %% find value for always migrate mechanism
    Pmodified=zeros(numStatesForStandardSolution,numStatesForStandardSolution);
    ck=zeros(numStatesForStandardSolution,1);
    actionsAlwaysMigrate=zeroStateIndexForStandardSolution*ones(1,numStatesForStandardSolution);
    for i=1:numStatesForStandardSolution
        Pmodified(i,:)=PForStandardSolution(actionsAlwaysMigrate(i),:);

        if Use2D==0
            ck(i)=(constFactorMigrate+proportionalFactorMigrate*powerFactor^abs(actionsAlwaysMigrate(i)-i))*(actionsAlwaysMigrate(i)~=i)... 
                + (constFactorTrans+proportionalFactorTrans*powerFactor^abs(actionsAlwaysMigrate(i)-zeroStateIndexForStandardSolution))*(actionsAlwaysMigrate(i)~=zeroStateIndexForStandardSolution);
        else
            tmpNorm1=norm(coordinatesCells2D(actionsAlwaysMigrate(i),:)-coordinatesCells2D(i,:))/cellDist;
            if abs(tmpNorm1-round(tmpNorm1))<0.00000001
                tmpNorm1=round(tmpNorm1);
            else
                tmpNorm1=ceil(tmpNorm1);
            end

            tmpNorm2=norm(coordinatesCells2D(actionsAlwaysMigrate(i),:)-coordinatesCells2D(zeroStateIndexForStandardSolution,:))/cellDist;
            if abs(tmpNorm2-round(tmpNorm2))<0.00000001
                tmpNorm2=round(tmpNorm2);
            else
                tmpNorm2=ceil(tmpNorm2);
            end

            ck(i)=(constFactorMigrate+proportionalFactorMigrate*powerFactor^tmpNorm1)*(actionsAlwaysMigrate(i)~=i)...
                + (constFactorTrans+proportionalFactorTrans*powerFactor^tmpNorm2)*(actionsAlwaysMigrate(i)~=zeroStateIndexForStandardSolution);
        end

    end

    valuesAlwaysMigrate=(eye(length(Pmodified))-gamma*Pmodified)\ck;    


    %% find value for myopic policy (migrate to state zero if one-timeslot transmission cost larger than miration cost, otherwise do not migrate)

    for i=1:numStatesForStandardSolution
        Pmodified(i,:)=PForStandardSolution(actionsMyopic(i),:);

        if Use2D==0
            ck(i)=(constFactorMigrate+proportionalFactorMigrate*powerFactor^abs(actionsMyopic(i)-i))*(actionsMyopic(i)~=i)... 
                + (constFactorTrans+proportionalFactorTrans*powerFactor^abs(actionsMyopic(i)-zeroStateIndexForStandardSolution))*(actionsMyopic(i)~=zeroStateIndexForStandardSolution);
        else
            tmpNorm1=norm(coordinatesCells2D(actionsMyopic(i),:)-coordinatesCells2D(i,:))/cellDist;
            if abs(tmpNorm1-round(tmpNorm1))<0.00000001
                tmpNorm1=round(tmpNorm1);
            else
                tmpNorm1=ceil(tmpNorm1);
            end

            tmpNorm2=norm(coordinatesCells2D(actionsMyopic(i),:)-coordinatesCells2D(zeroStateIndexForStandardSolution,:))/cellDist;
            if abs(tmpNorm2-round(tmpNorm2))<0.00000001
                tmpNorm2=round(tmpNorm2);
            else
                tmpNorm2=ceil(tmpNorm2);
            end

            ck(i)=(constFactorMigrate+proportionalFactorMigrate*powerFactor^tmpNorm1)*(actionsMyopic(i)~=i)...
                + (constFactorTrans+proportionalFactorTrans*powerFactor^tmpNorm2)*(actionsMyopic(i)~=zeroStateIndexForStandardSolution);
        end

    end

    valuesMyopic=(eye(length(Pmodified))-gamma*Pmodified)\ck;    


    %% find true value for modified policy iteration
    if Use2D==1
        %update actions
        actionsThPolicyMatchedToStandard=ones(numStatesForStandardSolution,1);

        for i=2:numStates2D
            if actionsThPolicy(i)==i
                for j=startingIndexEachCircle(i-1):startingIndexEachCircle(i)-1
                    actionsThPolicyMatchedToStandard(j)=j;
                end
            else
                for j=startingIndexEachCircle(i-1):startingIndexEachCircle(i)-1
                    if actionsThPolicy(i)==1
                        actionsThPolicyMatchedToStandard(j)=1;
                    else
                        tmpMinDist=inf;
                        tmpMinIndex=0;
                        for k=startingIndexEachCircle(actionsThPolicy(i)-1):startingIndexEachCircle(actionsThPolicy(i))-1
                            tmp=norm(coordinatesCells2D(j,:)-coordinatesCells2D(k,:))/cellDist;
                            if abs(tmp-round(tmp))<0.00000001
                                tmp=round(tmp);
                            else
                                tmp=ceil(tmp);
                            end

                            if tmp<tmpMinDist
                                tmpMinDist=tmp;
                                tmpMinIndex=k;
                            end
                        end

                        actionsThPolicyMatchedToStandard(j)=tmpMinIndex;
                    end
                end
            end
        end


        Pmodified=zeros(numStatesForStandardSolution,numStatesForStandardSolution);
        ck=zeros(numStatesForStandardSolution,1);

        for i=1:length(actionsThPolicyMatchedToStandard)
            Pmodified(i,:)=PForStandardSolution(actionsThPolicyMatchedToStandard(i),:);
            if Use2D==0
                ck(i)=(constFactorMigrate+proportionalFactorMigrate*powerFactor^abs(actionsThPolicyMatchedToStandard(i)-i))*(actionsThPolicyMatchedToStandard(i)~=i)... 
                    + (constFactorTrans+proportionalFactorTrans*powerFactor^abs(actionsThPolicyMatchedToStandard(i)-zeroStateIndexForStandardSolution))*(actionsThPolicyMatchedToStandard(i)~=zeroStateIndexForStandardSolution);
            else
                tmpNorm1=norm(coordinatesCells2D(actionsThPolicyMatchedToStandard(i),:)-coordinatesCells2D(i,:))/cellDist;
                if abs(tmpNorm1-round(tmpNorm1))<0.00000001
                    tmpNorm1=round(tmpNorm1);
                else
                    tmpNorm1=ceil(tmpNorm1);
                end

                tmpNorm2=norm(coordinatesCells2D(actionsThPolicyMatchedToStandard(i),:)-coordinatesCells2D(zeroStateIndexForStandardSolution,:))/cellDist;
                if abs(tmpNorm2-round(tmpNorm2))<0.00000001
                    tmpNorm2=round(tmpNorm2);
                else
                    tmpNorm2=ceil(tmpNorm2);
                end                

                ck(i)=(constFactorMigrate+proportionalFactorMigrate*powerFactor^tmpNorm1)*(actionsThPolicyMatchedToStandard(i)~=i)... 
                    + (constFactorTrans+proportionalFactorTrans*powerFactor^tmpNorm2)*(actionsThPolicyMatchedToStandard(i)~=zeroStateIndexForStandardSolution);
            end
        end        


        valuesThPolicyActual=(eye(length(Pmodified))-gamma*Pmodified)\ck;  
    end


    if supressOutput==0
        if notUseValueIteration==0
            disp(['time of value iteration: ', num2str(timeValue)])
        end
        if notUsePolicyIteration==0
            disp(['time of policy iteration: ', num2str(timePolicy)])
        end
        disp(['time of modified policy iteration: ', num2str(timeThPolicy)])
        disp(['average cost of never-migrate policy: ', num2str(sum(valuesNeverMigrate)/numStatesForStandardSolution)])
        disp(['average cost of always-migrate policy: ', num2str(sum(valuesAlwaysMigrate)/numStatesForStandardSolution)])
        disp(['average cost of myopic policy: ', num2str(sum(valuesMyopic)/numStatesForStandardSolution)])

        if notUseValueIteration==0
            disp(['average cost of value iteration: ', num2str(sum(valuesValueActual)/numStatesForStandardSolution)])
        end
        if notUsePolicyIteration==0
            disp(['average cost of policy iteration: ', num2str(sum(valuesPolicy)/numStatesForStandardSolution)])
        end

        if Use2D==1
            disp(['average cost of modified policy iteration: ', num2str(sum(valuesThPolicyActual)/numStatesForStandardSolution)])
        end
    end

    valueError=sum(abs(valuesPolicy-valuesValueResult));   %check value error
end