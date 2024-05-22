
matrix = genmatrix('rb_bravo.txt');
[G,H] = H2G(matrix);
%H=matrix;%
matrix=H;


num_v = size(H,2);
num_c = size(H,1);
num_iterations = 100;
num_runs=1000;

count = 0;

% set scaling
minSNR=-4;
SNRstepsize=0.5;
num_points=20;
maxSNR=minSNR+SNRstepsize*(num_points-1);

blers = zeros(1,num_points);
blers1= blers;
blers2= blers;
blers3= blers;
blers4= blers;
blers5= blers;
[r,k] = size(G);
parfor Noise_lvls=1:num_points
    error_counts = 0; error_countf=0; error_countm=0; error_countboft=0; error_countbb=0; error_countboff = 0;
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
    for n=1:num_iterations
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
    H2plusC = [Channel;H2];
    sums2 = sum(H2plusC);



          
    
    %% decision on codeword based on LLRS
    xhat2 = zeros(1,num_v);
    xhat2(sums2<0) = 1;
    s2 = (mod(xhat2*matrix.',2));
    syndrome_output2 = sum(s2);
      if syndrome_output2==0
      break
      end

    end
    HplusC = [Channel;H];
    sums = sum(HplusC);% was HplusC
       xhat = zeros(1,num_v);
    xhat(sums<0) = 1;
    err = sum(abs(xhat-codeword));
      if err~=0
        error_countf= error_countf+1;
      end

    HplusC = [Channel;H2];
    sums = sum(HplusC);% was HplusC
    xhat2 = zeros(1,num_v);
    xhat2(sums<0) = 1;
    err2 = sum(abs(xhat2-codeword));
      if err2~=0
        error_counts= error_counts+1;
      end

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
        syndrome_output3 = sum(s);   
        if syndrome_output3==0
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
        syndrome_output3 = sum(s);   
        if syndrome_output3==0
          break
        end
    end
    HplusC = [Channel;H3];
    sums = sum(HplusC);% was HplusC
    xhat3 = zeros(1,num_v);
    xhat3(sums<0) = 1;
    err3 = sum(abs(xhat3-codeword));
      if err3~=0
        error_countm= error_countm+1;
      end



    if err ~= 0 && err2 ~= 0 && err3 ~= 0
        error_countboft = error_countboft + 1;
    end
    
    xhat4 = xhat + xhat2 + xhat3;
    xhat4(xhat4<2) = 0;
    xhat4(xhat4>1) = 1;
    err4 = sum(abs(xhat4-codeword));
      if err4~=0
        error_countbb= error_countbb+1;
      end

    if err ~= 0 && err2 ~= 0 && err3 ~= 0 && err4 ~= 0
        error_countboff = error_countboff + 1;
    end
    end
    blers(Noise_lvls)= error_countf/num_runs;
    blers1(Noise_lvls)= error_counts/num_runs;
    blers2(Noise_lvls)= error_countm/num_runs;
    blers3(Noise_lvls)= error_countboft/num_runs;
    blers4(Noise_lvls)= error_countbb/num_runs;
    blers5(Noise_lvls)= error_countboff/num_runs;
    disp(count)
end


figure(1);

plot(minSNR:SNRstepsize:maxSNR,blers,"-o");
title(sprintf('BLER for %d blocks',num_runs));
hold on
plot(minSNR:SNRstepsize:maxSNR,blers1,'-x');
plot(minSNR:SNRstepsize:maxSNR,blers2,'-^');
plot(minSNR:SNRstepsize:maxSNR,blers3,'-+');
plot(minSNR:SNRstepsize:maxSNR,blers4,'-square');
plot(minSNR:SNRstepsize:maxSNR,blers5,'-diamond');
legend("Flooding","Serial","Mix", "Best of three", "Best bit", "Best of foour");
hold off
