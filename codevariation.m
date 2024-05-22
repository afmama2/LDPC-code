
matrix = genmatrix('rb_bravo.txt');
[G,H] = H2G(matrix);
matrix = genmatrix('rb_alpha.txt');
[G1,H1] = H2G(matrix);
matrix = genmatrix('r_alpha.txt');
[G2,H2] = H2G(matrix);
matrix = genmatrix('r_bravo.txt');
[G3,H3] = H2G(matrix);
%H=matrix;%
matrix=H;
matrix1=H1;
matrix2=H2;
matrix3=H3;

size(H)
size(H1)
size(H2)
size(H3)
num_v = size(H,2);
num_c = size(H,1);
num_v1 = size(H1,2);
num_c1 = size(H1,1);
num_v2 = size(H2,2);
num_c2 = size(H2,1);
num_v3 = size(H3,2);
num_c3 = size(H3,1);
num_iterations = 100;
num_runs=1000;


count = 0;

% set scaling
minSNR=-4;
SNRstepsize=0.5;
num_points=20;
maxSNR=minSNR+SNRstepsize*(num_points-1);

Ber = zeros(1,num_points);
Ber2 = zeros(1,num_points);
Ber3 = Ber;
Ber1 = Ber;
Ber0 = Ber;
Ber10 = Ber;
Ber20 = Ber;
Ber30 = Ber;
[r,k] = size(G);
[r1,k1] = size(G1);
[r2,k2] = size(G2);
[r3,k3] = size(G3);


parfor Noise_lvls=1:num_points
    ber = zeros(1,num_runs);
    ber2 = zeros(1,num_runs);
    ber1 = ber;
    ber3 = ber;
    ber0 = ber;
    ber10 = ber;
    ber20 = ber;
    ber30 = ber;
    for l=1:num_runs
%% Generate Codeword using G Matrix (Note:13/1/2017 unsure if G is correct)

 X = randi([0 1],1,r); 
 codeword = mod(X*G,2);
  X = randi([0 1],1,r1); 
 codeword1 = mod(X*G1,2);
  X = randi([0 1],1,r2); 
 codeword2 = mod(X*G2,2);
  X = randi([0 1],1,r3); 
 codeword3 = mod(X*G3,2);
    %% Generate random noise variance in accordance with Channel model and modulation
    
SNR = minSNR+((Noise_lvls-1)*SNRstepsize);
     N0 = 1/(10^(SNR/10));

     a_tx = -2*(codeword-0.5);
     a_rx= a_tx + sqrt(N0/2)*(randn(size(a_tx))+1i*randn(size(a_tx)));
     Channel = (abs(a_rx+1).^2-abs(a_rx-1).^2)/N0;

     a_tx = -2*(codeword1-0.5);
     a_rx= a_tx + sqrt(N0/2)*(randn(size(a_tx))+1i*randn(size(a_tx)));
     Channel1 = (abs(a_rx+1).^2-abs(a_rx-1).^2)/N0;

     a_tx = -2*(codeword2-0.5);
     a_rx= a_tx + sqrt(N0/2)*(randn(size(a_tx))+1i*randn(size(a_tx)));
     Channel2 = (abs(a_rx+1).^2-abs(a_rx-1).^2)/N0;

     a_tx = -2*(codeword3-0.5);
     a_rx= a_tx + sqrt(N0/2)*(randn(size(a_tx))+1i*randn(size(a_tx)));
     Channel3 = (abs(a_rx+1).^2-abs(a_rx-1).^2)/N0;
    

    %% Load LLRs into Parity Check
    H = zeros(num_c,num_v);
    H0 = H;
    H1 = zeros(num_c1,num_v1);
    H10 = H1;
    H2 = zeros(num_c2,num_v2);
    H20 = H2;
    H3 = zeros(num_c3,num_v3);
    H30 = H3;
  
    %% start iterative decoding
   
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

    HplusC = [Channel;H];
    sums = sum(HplusC);





    %% decision on codeword based on LLRS
    xhat = zeros(1,num_v);
    xhat(sums<0) = 1;
    s = (mod(xhat*matrix.',2));
    syndrome_output = sum(s);
      if syndrome_output==0
      break
      end

    end
    for n = 1:num_iterations

        for c = 1:num_c
            int_LLrs = [];
            for v = 1:num_v
                if matrix(c,v)==1
                    calc = Channel(v) + sum(H0(:,v))-H0(c,v);
                    int_LLrs = [int_LLrs calc];
                end
            end                
            pri_LLRs = check_node(int_LLrs);
            for v = 1:num_v
                if matrix(c,v)==1
                    H0(c,v) = pri_LLRs(1);
                    pri_LLRs(1) = [];
                end
            end
        
        end
        HplusC = [Channel;H0];
        sums = sum(HplusC);





        %% decision on codeword based on LLRS
        xhat = zeros(1,num_v);
        xhat(sums<0) = 1;
        s = (mod(xhat*matrix.',2));
        syndrome_output2 = sum(s);
      if syndrome_output2==0
      break
      end
    end
    for n = 1:num_iterations
%         tic
        %% variable node calculations
        
        mat1 = [Channel1;H1];
        mat2 = repmat(sum(mat1),num_c1,1);
        mat2(matrix1==0) = 0;
        H1 = mat2 - H1;
        %disp(['end of iteration ',num2str(n)]);
%         toc
        
        %% check node calculations
        for c = 1:num_c1
            int_LLRs = [];
            for v = 1:num_v1
                if matrix1(c,v)==1
                    int_LLRs = [int_LLRs H1(c,v)];
                end
            end
           % int_LLRs
            pri_LLRs = check_node(int_LLRs);
            
           % pri_LLRs
           
            for v = 1:num_v1
          
                if matrix1(c,v)==1
                    H1(c,v) = pri_LLRs(1);
                    pri_LLRs(1) = [];
                  
                end
            end
        end

    HplusC = [Channel1;H1];
    sums = sum(HplusC);





    %% decision on codeword based on LLRS
    xhat = zeros(1,num_v1);
    xhat(sums<0) = 1;
    s = (mod(xhat*matrix1.',2));
    syndrome_output = sum(s);
      if syndrome_output==0
      break
      end

    end
    for n = 1:num_iterations

        for c = 1:num_c1
            int_LLrs = [];
            for v = 1:num_v1
                if matrix1(c,v)==1
                    calc = Channel1(v) + sum(H10(:,v))-H10(c,v);
                    int_LLrs = [int_LLrs calc];
                end
            end                
            pri_LLRs = check_node(int_LLrs);
            for v = 1:num_v1
                if matrix1(c,v)==1
                    H10(c,v) = pri_LLRs(1);
                    pri_LLRs(1) = [];
                end
            end
        
        end
        HplusC = [Channel1;H10];
        sums = sum(HplusC);





        %% decision on codeword based on LLRS
        xhat = zeros(1,num_v1);
        xhat(sums<0) = 1;
        s = (mod(xhat*matrix1.',2));
        syndrome_output2 = sum(s);
      if syndrome_output2==0
      break
      end
    end
    for n = 1:num_iterations
%         tic
        %% variable node calculations
        
        mat1 = [Channel2;H2];
        mat2 = repmat(sum(mat1),num_c2,1);
        mat2(matrix2==0) = 0;
        H2 = mat2 - H2;
        %disp(['end of iteration ',num2str(n)]);
%         toc
        
        %% check node calculations
        for c = 1:num_c2
            int_LLRs = [];
            for v = 1:num_v2
                if matrix2(c,v)==1
                    int_LLRs = [int_LLRs H2(c,v)];
                end
            end
           % int_LLRs
            pri_LLRs = check_node(int_LLRs);
            
           % pri_LLRs
           
            for v = 1:num_v2
          
                if matrix2(c,v)==1
                    H2(c,v) = pri_LLRs(1);
                    pri_LLRs(1) = [];
                  
                end
            end
        end

    HplusC = [Channel2;H2];
    sums = sum(HplusC);





    %% decision on codeword based on LLRS
    xhat = zeros(1,num_v2);
    xhat(sums<0) = 1;
    s = (mod(xhat*matrix2.',2));
    syndrome_output = sum(s);
      if syndrome_output==0
      break
      end

    end
    for n = 1:num_iterations

        for c = 1:num_c2
            int_LLrs = [];
            for v = 1:num_v2
                if matrix2(c,v)==1
                    calc = Channel2(v) + sum(H20(:,v))-H20(c,v);
                    int_LLrs = [int_LLrs calc];
                end
            end                
            pri_LLRs = check_node(int_LLrs);
            for v = 1:num_v2
                if matrix2(c,v)==1
                    H20(c,v) = pri_LLRs(1);
                    pri_LLRs(1) = [];
                end
            end
        
        end
        HplusC = [Channel2;H20];
        sums = sum(HplusC);





        %% decision on codeword based on LLRS
        xhat = zeros(1,num_v2);
        xhat(sums<0) = 1;
        s = (mod(xhat*matrix2.',2));
        syndrome_output2 = sum(s);
      if syndrome_output2==0
      break
      end
    end
    for n = 1:num_iterations
%         tic
        %% variable node calculations
        
        mat1 = [Channel3;H3];
        mat2 = repmat(sum(mat1),num_c3,1);
        mat2(matrix3==0) = 0;
        H3 = mat2 - H3;
        %disp(['end of iteration ',num2str(n)]);
%         toc
        
        %% check node calculations
        for c = 1:num_c3
            int_LLRs = [];
            for v = 1:num_v3
                if matrix3(c,v)==1
                    int_LLRs = [int_LLRs H3(c,v)];
                end
            end
           % int_LLRs
            pri_LLRs = check_node(int_LLRs);
            
           % pri_LLRs
           
            for v = 1:num_v3
          
                if matrix3(c,v)==1
                    H3(c,v) = pri_LLRs(1);
                    pri_LLRs(1) = [];
                  
                end
            end
        end
        
    HplusC = [Channel3;H3];
    sums = sum(HplusC);





    %% decision on codeword based on LLRS
    xhat = zeros(1,num_v3);
    xhat(sums<0) = 1;
    s = (mod(xhat*matrix3.',2));
    syndrome_output = sum(s);
      if syndrome_output==0
      break
      end

    end
    for n = 1:num_iterations

        for c = 1:num_c3
            int_LLrs = [];
            for v = 1:num_v3
                if matrix3(c,v)==1
                    calc = Channel3(v) + sum(H30(:,v))-H30(c,v);
                    int_LLrs = [int_LLrs calc];
                end
            end                
            pri_LLRs = check_node(int_LLrs);
            for v = 1:num_v3
                if matrix3(c,v)==1
                    H30(c,v) = pri_LLRs(1);
                    pri_LLRs(1) = [];
                end
            end
        
        end
        HplusC = [Channel3;H30];
        sums = sum(HplusC);





        %% decision on codeword based on LLRS
        xhat = zeros(1,num_v3);
        xhat(sums<0) = 1;
        s = (mod(xhat*matrix3.',2));
        syndrome_output2 = sum(s);
      if syndrome_output2==0
      break
      end
    end
    HplusC = [Channel;H];
    sums = sum(HplusC);% was HplusC
    xhat = zeros(1,num_v);
    xhat(sums<0) = 1;
    err = sum(abs(xhat-codeword));
    ber(l) = err;

    HplusC = [Channel;H0];
    sums = sum(HplusC);% was HplusC
    xhat = zeros(1,num_v);
    xhat(sums<0) = 1;
    err = sum(abs(xhat-codeword));
    ber0(l) = err;

    HplusC = [Channel1;H1];
    sums = sum(HplusC);% was HplusC
    xhat1 = zeros(1,num_v1);
    xhat1(sums<0) = 1;
    err1 = sum(abs(xhat1-codeword1));
    ber1(l) = err1;

    HplusC = [Channel1;H10];
    sums = sum(HplusC);% was HplusC
    xhat1 = zeros(1,num_v1);
    xhat1(sums<0) = 1;
    err1 = sum(abs(xhat1-codeword1));
    ber10(l) = err1;

    HplusC = [Channel2;H2];
    sums = sum(HplusC);% was HplusC
    xhat2 = zeros(1,num_v2);
    xhat2(sums<0) = 1;
    err2 = sum(abs(xhat2-codeword2));
    ber2(l) = err2;

    HplusC = [Channel2;H20];
    sums = sum(HplusC);% was HplusC
    xhat2 = zeros(1,num_v2);
    xhat2(sums<0) = 1;
    err2 = sum(abs(xhat2-codeword2));
    ber20(l) = err2;

    H3plusC = [Channel3;H3];
    sums = sum(H3plusC);% was HplusC
    xhat3 = zeros(1,num_v3);
    xhat3(sums<0) = 1;
    err3 = sum(abs(xhat3-codeword3));
    ber3(l) = err3;

    H3plusC = [Channel3;H30];
    sums = sum(H3plusC);% was HplusC
    xhat3 = zeros(1,num_v3);
    xhat3(sums<0) = 1;
    err3 = sum(abs(xhat3-codeword3));
    ber30(l) = err3;

    end
    Ber(Noise_lvls) = sum(ber)/(num_runs*k);
    Ber0(Noise_lvls) = sum(ber0)/(num_runs*k);
    Ber2(Noise_lvls) = sum(ber2)/(num_runs*k2);
    Ber20(Noise_lvls) = sum(ber20)/(num_runs*k2);
    Ber1(Noise_lvls) = sum(ber1)/(num_runs*k1);
    Ber3(Noise_lvls) = sum(ber3)/(num_runs*k3);
    Ber30(Noise_lvls) = sum(ber30)/(num_runs*k3);
    Ber10(Noise_lvls) = sum(ber10)/(num_runs*k1);
    disp(count)
end


zerolocs=(find(Ber==0));
Ber(zerolocs)=1e-10;


zerolocs=(find(Ber2==0));
Ber2(zerolocs)=1e-10;

zerolocs=(find(Ber3==0));
Ber3(zerolocs)=1e-10;

zerolocs=(find(Ber1==0));
Ber1(zerolocs)=1e-10;

zerolocs=(find(Ber10==0));
Ber10(zerolocs)=1e-10;

zerolocs=(find(Ber20==0));
Ber20(zerolocs)=1e-10;

zerolocs=(find(Ber30==0));
Ber30(zerolocs)=1e-10;

zerolocs=(find(Ber0==0));
Ber0(zerolocs)=1e-10;

figure(1);

semilogy(minSNR:SNRstepsize:maxSNR,Ber,"-square");
title(sprintf("BER after %d iteration",num_iterations))
hold on
semilogy(minSNR:SNRstepsize:maxSNR,Ber0);
semilogy(minSNR:SNRstepsize:maxSNR,Ber1,"-o");
semilogy(minSNR:SNRstepsize:maxSNR,Ber10);
semilogy(minSNR:SNRstepsize:maxSNR,Ber2,"-diamond");
semilogy(minSNR:SNRstepsize:maxSNR,Ber20);
semilogy(minSNR:SNRstepsize:maxSNR,Ber3,"-x");
semilogy(minSNR:SNRstepsize:maxSNR,Ber30);
legend("code (980,717) flooding", "code (980,717) serial", "code (980,490) flooding", "code (980,490) serial", "code (847,484) flooding", "code (847,484) serial", "code (847,681) flooding", "code (847,681) serial");
hold off
