
matrix = genmatrix('rb_bravo.txt');
[G,H] = H2G(matrix);
%H=matrix;%
matrix=H;


num_v = size(H,2);
num_c = size(H,1);
num_iterations = 10000;
num_runs=1;
count_s_correct = 0;
count_s_wrong = 0;
count_f_correct = 0;
count_f_wrong = 0;
count_m_correct = 0;
count_m_wrong = 0;
count = 0;

% set scaling
minSNR=1;
SNRstepsize=0.05;
num_points=1;
maxSNR=minSNR+SNRstepsize*(num_points-1);

Ber = zeros(1,num_points);
Ber2 = zeros(1,num_points);
Ber3 = Ber;
Raw  = Ber;

[r,k] = size(G);

for Noise_lvls=1:num_points
    ber = zeros(1,num_runs);
    ber2 = zeros(1,num_runs);
    raw = ber;
    ber3=ber;
    for l=1:num_runs
%% Generate Codeword using G Matrix (Note:13/1/2017 unsure if G is correct)

 X = randi([0 1],1,r); 
 codeword = mod(X*G,2);
    %% Generate random noise variance in accordance with Channel model and modulation
    
SNR = minSNR+((Noise_lvls-1)*SNRstepsize);
     N0 = 1/(10^(SNR/10));

     a_tx = -2*(codeword-0.5);
     a_rx= a_tx + sqrt(N0/2)*(randn(size(a_tx))+1i*randn(size(a_tx)));
     Channel = (abs(a_rx+1).^2-abs(a_rx-1).^2)/N0;
    

    %% Load LLRs into Parity Check
    H = zeros(num_c,num_v);
    H2 = H;
    H3 = H;
    % H3v = H;
  
    %% start iterative decoding
   tic 
    for n = 1:num_iterations
%         tic
        %% variable node calculations
        
        mat1 = [Channel;H];
        mat2 = repmat(sum(mat1),num_c,1);
        mat2(matrix==0) = 0;
        H = mat2 - H;
        %disp(['end of iteration ',num2str(n)]);
%         toc
        
        %% check node calculations
        for c = 1:num_c
            int_LLRs = [];
            for v = 1:num_v
                if matrix(c,v)==1
                    int_LLRs = [int_LLRs H(c,v)];
                end
            end
           % int_LLRs
            pri_LLRs = check_node(int_LLRs);
            
           % pri_LLRs
           
            for v = 1:num_v
          
                if matrix(c,v)==1
                    H(c,v) = pri_LLRs(1);
                    pri_LLRs(1) = [];
                  
                end
            end
        end
    % 
    % HplusC = [Channel;H];
    % sums = sum(HplusC);
    % 
    % 
    % 
    % 
    % 
    % %% decision on codeword based on LLRS
    % xhat = zeros(1,num_v);
    % xhat(sums<0) = 1;
    % s = (mod(xhat*matrix.',2));
    % syndrome_output = sum(s);
      % if syndrome_output==0
      % break
      % end

    end
    toc
    tic
    for n = 1:num_iterations

        for c = 1:num_c
            int_LLrs = [];
            for v = 1:num_v
                if matrix(c,v)==1
                    calc = Channel(v) + sum(H2(:,v))-H2(c,v);
                    int_LLrs = [int_LLrs calc];
                end
            end                
            pri_LLRs = check_node(int_LLrs);
            for v = 1:num_v
                if matrix(c,v)==1
                    H2(c,v) = pri_LLRs(1);
                    pri_LLRs(1) = [];
                end
            end
        
        end
        % HplusC = [Channel;H2];
        % sums = sum(HplusC);
        % 
        % 
        % 
        % 
        % 
        % %% decision on codeword based on LLRS
        % xhat = zeros(1,num_v);
        % xhat(sums<0) = 1;
        % s = (mod(xhat*matrix.',2));
        % syndrome_output = sum(s);
      % if syndrome_output==0
      % break
      % end
    end
    toc
    HplusC = [Channel;H];
    sums = sum(HplusC);% was HplusC
    xhat = zeros(1,num_v);
    xhat(sums<0) = 1;
    err = sum(abs(xhat-codeword));
    ber(l) = err;
    if sum(mod(xhat*matrix.',2)) == 0
        if err == 0
            count_f_correct = count_f_correct + 1;
        else
            count_f_wrong = count_f_wrong + 1;
        end
    end

    HplusC = [Channel;H2];
    sums = sum(HplusC);% was HplusC
    xhat = zeros(1,num_v);
    xhat(sums<0) = 1;
    err = sum(abs(xhat-codeword));
    ber2(l) = err;
    if sum(mod(xhat*matrix.',2)) == 0
        if err == 0
            count_s_correct = count_s_correct + 1;
        else
            count_s_wrong = count_s_wrong + 1;
        end
    end

    xhat = zeros(1,num_v);
    xhat(Channel<0) = 1;
    err = sum(abs(xhat-codeword));
    raw(l) = err;

    for n = 1:num_iterations/2

        mat1 = [Channel;H3];
        mat2 = repmat(sum(mat1),num_c,1);
        mat2(matrix==0) = 0;
        H3 = mat2 - H3;
        %disp(['end of iteration ',num2str(n)]);
%         toc
            
            %% check node calculations
        for c = 1:num_c
            int_LLRs = [];
            for v = 1:num_v     
                if matrix(c,v)==1
                    int_LLRs = [int_LLRs H3(c,v)];
                    
                end
            end
           % int_LLRs
            pri_LLRs = check_node(int_LLRs);
            
           % pri_LLRs
           
            for v = 1:num_v
          
                if matrix(c,v)==1
                    H3(c,v) = pri_LLRs(1);
                    pri_LLRs(1) = [];
                  
                end
            end
        end
        H3plusC = [Channel;H3];
        sums3 = sum(H3plusC);
        xhat3 = zeros(1,num_v);
        xhat3(sums3<0) = 1;
    
        s = (mod(xhat3*matrix.',2));
        syndrome_output = sum(s);   
        if syndrome_output==0
          break
        end

        for c = 1:num_c
            int_LLrs = [];
            for v = 1:num_v
                if matrix(c,v)==1
                    calc = Channel(v) + sum(H3(:,v))-H3(c,v);
                    int_LLrs = [int_LLrs calc];                   
                end
            end                
            pri_LLRs = check_node(int_LLrs);
            for v = 1:num_v
                if matrix(c,v)==1
                    H3(c,v) = pri_LLRs(1);
                    pri_LLRs(1) = [];
                end
            end
        
        end
        
        H3plusC = [Channel;H3];
        sums3 = sum(H3plusC);
        xhat3 = zeros(1,num_v);
        xhat3(sums3<0) = 1;
    
        s = (mod(xhat3*matrix.',2));
        syndrome_output = sum(s);   
        if syndrome_output==0
          break
        end
    end
    H3plusC = [Channel;H3];
    sums = sum(H3plusC);% was HplusC
    xhat = zeros(1,num_v);
    xhat(sums<0) = 1;
    err = sum(abs(xhat-codeword));
    ber3(l) = err;
    if sum(mod(xhat*matrix.',2)) == 0
        if err == 0
            count_m_correct = count_m_correct + 1;
        else
            count_m_wrong = count_m_wrong + 1;
        end
    end

    count = count + 1;
    count
    end
    Ber(Noise_lvls) = sum(ber)/(num_runs*k);
    Ber2(Noise_lvls) = sum(ber2)/(num_runs*k);
    Raw(Noise_lvls) = sum(raw)/(num_runs*k);
    Ber3(Noise_lvls) = sum(ber3)/(num_runs*k);
end


zerolocs=(find(Ber==0));
Ber(zerolocs)=1e-10;


zerolocs=(find(Ber2==0));
Ber2(zerolocs)=1e-10;
% 
 zerolocs=(find(Ber3==0));
 Ber3(zerolocs)=1e-10;

disp(["Correct serial count:" count_s_correct]);
disp(["Wrong serial count:" count_s_wrong]);
disp(["Correct flooding count:" count_f_correct]);
disp(["Wrong flooding count:" count_f_wrong]);
disp(["Wrong mix count:" count_m_wrong]);
disp(["Correct mix count:" count_m_correct]);

figure(2);

semilogy(minSNR:SNRstepsize:maxSNR,Ber,"-square");
title(sprintf("BER after %d iteration",num_iterations))
hold on
semilogy(minSNR:SNRstepsize:maxSNR,Ber2,"-diamond");
semilogy(minSNR:SNRstepsize:maxSNR,Raw,"-o");
semilogy(minSNR:SNRstepsize:maxSNR,Ber3,"-x");
legend("Flooding","Serial","Raw","Mix");
hold off
