%% Title: DISTRIBUTED SENSOR NETWORK LOCALIZATION USING SOCP RELAXATION
% We have 500 nodes, m of which are anchors and the rest are sensors, the
% percentage of anchor nodes being p.
 function error = localization_socp(n,p,radioRange,nfd,nfa)

clc
%n = 300; p = .6; 
m = (1-p)*n; 
%radioRange = 0.5;
%% generate nodes with true coordinates xt and yt and plot, calculate dt(n,n) matrix, setA matrix, cardA
%Generate random true node positions, using n samples from a uniform distribution: a is lower limit and b is upper limit of region in which 
% nodes are generated.
a = -0.5; b = 0.5;
xt = a + rand(n,1)*(b-a);
yt = a + rand(n,1)*(b-a);
xt=xt.*(max(0, 1+randn(1,1)*nfa)) ;
yt=yt.*(max(0, 1+randn(1,1)*nfa));
scatter(xt(1:m),yt(1:m),'MarkerFaceColor','c');                    % first m nodes are the sensors, to be estimated
axis([-0.5 0.5 -0.5 0.5]); grid on; hold on;
%scatter(xt(m+1:n),yt(m+1:n),'d');

% Calculate distances
dt = zeros(n,n);
cardA = 0; % Cardinality of set A;
% 'setA': set of all neighbor pairs (i,j) which are
% within radio range of each other
for i = 1:n
    for j = 1:n
        if (i~=j) %&& (j>i)
            dt(i,j) = sqrt((xt(i)-xt(j))^2 + (yt(i)-yt(j))^2);
            if (dt(i,j) < radioRange)
                setA(cardA+1,:) = [i,j];                         
                cardA = cardA + 1;
            end
        end
    end
end

%setA = [ xi1 yi1 ;  coordinates of nodes for whom dt < radioRange
%         xi2 yi2 ;
%         .
%         .
%         .       ;]

%% Noisy Measured or estimated distances between neighbors
%nfd = 0.0; % Noise factor
noise = max(0, 1+randn(n)*nfd);
d = dt.*noise;

%% Initialize Node and Anchor positions
x = zeros(2*n,1);
x(2*m+1:2:2*n-1,1) = xt(m+1:n,1);
x(2*m+2:2:2*n,1) = yt(m+1:n,1);

MAXITER = 100; ABSTOL = 0.01;
for niter = 1:MAXITER
    prevx = x(1:2*m);
    for i = 1:m % Loop over all sensor nodes to be positioned
        % find the neighbors of Node 'i'
        colA = setA(:,1);  %1st column of setA, x coordinates of feasible points
        
        neighbors = colA == i; %find indices of nodes with x-coordinate i among feasible points
        NA = setA(neighbors,:);      %coordinates of neighbors of node i
        j = NA(:,2); % y coordinates of Neighbors of Node 'i'
        ni = length(j); % Number of Neighbors of Node 'i'
        %ni_max = max(ni,ni_max);
        ti = zeros(n*(n-1)/2,1);
        for na = 1:length(j)
            ij = sum(n:-1:(n-i+2)) + (j(na)-i);
            ti(ij) = 1;
        end
        b = -[0; zeros(2*(ni+1),1); ones(ni,1); zeros(ni,1)];
        S = []; r = []; U = []; t = [];
        for na = 1:length(j)
            tij = zeros(1,ni);
            tij(na) = 1;
            Sj = [zeros(1,2*ni+3) tij zeros(1,ni);
                zeros(1,3*ni+3) tij];
            rj = [0;d(i,j(na))^2];
            S = [S; Sj]; r = [r; rj];
            ei1 = zeros(1,2*ni+2); ei1(1) = 1;
            ei2 = zeros(1,2*ni+2); ei2(2) = 1;
            Uj = [0.5 zeros(1,3*ni+2) 0.5*tij;
                -0.5 zeros(1,3*ni+2) 0.5*tij;
                0  ei1 zeros(1,2*ni);
                0  ei2 zeros(1,2*ni)];
            tj = [0; 0; x(2*j(na)-1); x(2*j(na))];
            U = [U; Uj]; t = [t; tj];
        end
        Li = [1 zeros(1,4*ni+2)];
        % Formulating the problem in dual form
        A = -[Li; S; U]'; c = -[1; r; t];
        K.f = 1; pars.fid = 0;
        K.q = [2*ones(1,length(j)) 4*ones(1,length(j))];
        if ni>1
            % Optimization performed only if number of neighbors
            % is greater than 1
            [X,Y,info] = sedumi(A,b,c,K,pars);
            xx(:,i) = Y(2:3);
        else
            xx(:,i) = x(2*i-1:2*i);
        end
    end

    % Update node positions based on the computation in this phase
    for i = 1:m
        x(2*i-1:2*i) = xx(:,i);
    end
    % Stopping Criterion
    change = prevx - x(1:2*m);
    if max(abs(change)) <= ABSTOL
          % Estimated Sensor Positions
        xet = x(1:2:2*m-1); yet = x(2:2:2*m);
        scatter(xet, yet, 'r+');
            % Add error lines from the actual to the estimated positions
        for i=1:m
        line([xt(i), xet(i)],[yt(i), yet(i)]);
        end
    error=[xt(1:m) - xet, yt(1:m) - yet]; 
    error=error.^2 ; 
    error=error(:,1)+error(:,2);
    error=mean(error);
    error=error.^.5;  
    
        break;
    end   
end