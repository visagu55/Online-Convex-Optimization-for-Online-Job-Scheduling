t = 150;          % the length of time this program estimate
k = 1;           % the flow time and fractional flow time is in Lk-norm
numofjob = 100;   % how many job in this program
lamda = 4000;      % initial lamda
time = 1:t;

job = zeros(numofjob,8); % id, arrive time, size , cpu, finish rate, whether arrive, finish time of oco, finish time of SJF 
schedule = zeros(numofjob,t); % the scheduling result of this oco method
sjbschedule = zeros(numofjob,t); % the scheduling of SJF(shortest job first)
fft = zeros(1,numofjob); % fractional flow time
ocofractional = zeros(t,1);  % the fractional flow time of each time
ocoflow = zeros(t,1);        % the flow time of each time
sjbfractional = zeros(t,1);  % the fractional flow time of shortest job first
sjbflow = zeros(t,1);        % the flow time of shortest job first

%initial the information of each job
for i = 1:numofjob
    job(i,1) = i;
    if i == 1
        job(i,2) = 0;
    end
    if i > 1
        job(i,2) = job(i-1,2) + round(rand(1));
    end
    job(i,3) = (round(rand(1)*5)+1)/2;
    job(i,4) = round(rand(1),1);
end

% this part is for online convex optimization
for j = 1:t
    for i = 1:numofjob
        if job(i,2) >= j
            fft(1,i) = 0;      
            job(i,6) = 0;            %not arrive
        else
            fft(1,i) = (j-job(i,2))^k/job(i,3) + job(i,3)^(k-1);
            job(i,6) = 1;            %arrive
        end
    end
    
    % solve the offline convex optimization
    cvx_begin
        variable x(numofjob)
        minimize( fft * x + lamda * job(:,6)' * ((job(:,3)-job(:,5)-x ).* (t*ones(numofjob,1)-job(:,2))) + x'*x/100)
        subject to
            job(:,4)'* x <= 1;
            x <= ones(numofjob,1);
            x >= zeros(numofjob,1);
            x <= job(:,3)-job(:,5);
    cvx_end
    
    % update the finish rate of each job
    x = round(x,1);
    
    job(:,5) = job(:,5) + x;
    
    %update lamda
    lamda = lamda +  k*5 * job(:,6)' * (job(:,3)-job(:,5)-x );
    schedule(:,j) = x;
    
    %compute the fractional flow time of each job
    if j >1
        ocofractional(j,1) = ocofractional(j-1,1) + fft * x;
    else
        ocofractional(j,1) =  fft * x;
    end
    
    % compute the finish time of each job
    for i = 1:numofjob
        if x(i) > 0
            if job(i,5) >= job(i,3)
                job(i,7) = j;
            end
        end
    end
    
    % compute the flow time
    for i = 1:numofjob
        if job(i,7) > 0 && job(i,7) <= j
            ocoflow(j,1) = ocoflow(j,1) + (job(i,7) - job(i,2))^k;
        end
    end                                
end

%this part is for shortest job first
sjbfinish = zeros(numofjob,2); % finsih rate, remaining rate
sjbfinish(:,2) = job(:,3);

for j = 1:t
    rcpu = 1;    %remaining cpu can be used
    scpu = 1;    % the shortest cpu of all the job
    
    % find the shortest memory job
    for i = 1:numofjob
        if job(i,2) < j && sjbfinish(i,2) > 0      %arrive and need to work
            if scpu == 1
                scpu = job(i,4);
            else
                if scpu > job(i,4)
                    scpu = job(i,4);
                end
            end
        end
    end

    numofrunning=0;
    tmpflag = zeros(numofjob,1);
    % find which job is the shortest and can run
    while(true)
        numofrunning = numofrunning +1;
        flag = 0;
        
        for i = 1:numofjob
            if job(i,2) < j && sjbfinish(i,2) > 0 && job(i,4) <= rcpu && tmpflag(i,1) ~= 1     %arrive and need to£¨can) work
                if flag == 0
                    flag = i;
                else
                    if sjbfinish(i,2) < sjbfinish(flag,2)
                        flag = i;
                    end
                end    
            end
        end
    
        if flag > 0
            sjbschedule(flag,j) = 1;
            rcpu = rcpu - job(flag,4);
            
            tmpflag(flag,1) = 1;
            
        end
                
        %compute the shortest job
        if flag > 0
            sjbfinish(flag,1) = sjbfinish(flag,1) + 1;
            
            if sjbfinish(flag,1) >= job(flag,3)
                sjbfinish(flag,1) = job(flag,3);
                job(flag,8) = j;
            end
        
            sjbfinish(flag,2) = job(flag,3) - sjbfinish(flag,1);
        end
        
        
        if rcpu < scpu
           break
        end
        
        if numofrunning >=10
            break
        end
    end
     
    
     %compute the fractional flow time of shortest job first
     if j ==1
         for i = 1:numofjob
             sjbfractional(j,1) = sjbfractional(j,1) + tmpflag(i,1) * ((j -job(i,2))^k/job(i,3) + job(i,3)^(k-1));
         end   
     end
     
     if j > 1
         sjbfractional(j,1) = sjbfractional(j-1,1);
         for i = 1:numofjob
             sjbfractional(j,1) = sjbfractional(j,1) + tmpflag(i,1) * tmpflag(i,1) * ((j -job(i,2))^k/job(i,3) + job(i,3)^(k-1));
         end
     end
     
     
     %compute the flow time of the shortest job first
     if j == 1        
         for i = 1:numofjob
             if sjbfinish(i,2) == 0
                 sjbflow(j,1) = sjbflow(j,1) + 1;
             end
         end   
     end
     
     if j > 1  
         sjbflow(j,1) = sjbflow(j-1,1);
         for i = 1:numofjob
             if sjbfinish(i,2) == 0
                 sjbflow(j,1) = sjbflow(j,1) + tmpflag(i,1) * (j - job(i,2))^k;       
             end    
         end 
     end
end

            
% plot the flow time and fractional flow time of our method and shortest job first    
plot(time, ocoflow,'r-', 'LineWidth', 2)
hold on
plot(time, ocofractional,'b-', 'LineWidth', 2)
plot(time,sjbflow,'r--', 'LineWidth', 2)
plot(time,sjbfractional,'b--', 'LineWidth', 2)
xlabel('Time slot', 'FontSize', 16, 'FontWeight','bold')
ylabel('Overall flow time', 'FontSize', 16, 'FontWeight','bold')
% legend({'OCO flowtime', 'SJB flowtime'}, 'FontSize', 12)
set(gca,'fontsize',12)
legend('OCO flowtime', 'OCO fractional flowtime', 'SJF flowtime', 'SJF fractional flowtime', 'FontSize', 24)