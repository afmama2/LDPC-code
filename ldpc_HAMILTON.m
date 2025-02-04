
matrix = genmatrix('rb_alpha.txt');
[G,H] = H2G(matrix);
%H=matrix;%
matrix=H;
%
int_LLRs=zeros(size(H));

correctcount = 0;
falsecodeword=0;
LDPCuseless=0;
num_v = size(H,2);
num_c = size(H,1);
num_iterations = 100;
syndrome_output=10;
num_runs=10;

% set scaling
minSNR=1;
SNRstepsize=0.5;
num_points=2;
maxSNR=minSNR+SNRstepsize*(num_points-1);

[r,k] = size(G);
for Noise_lvls=1:num_points
    correctcount=0;
    for l=1:num_runs
%% Generate Codeword using G Matrix (Note:13/1/2017 unsure if G is correct)

%To be copied
 X = randi([0 1],1,r); 
 codeword = mod(X*G,2);
 %codeword = DFENCRB;
 
    %% Generate random noise variance in accordance with Channel model and modulation
    
SNR = minSNR+((Noise_lvls-1)*SNRstepsize);
     N0 = 1/(10^(SNR/10));

     a_tx = -2*(codeword-0.5);
     a_rx= a_tx + sqrt(N0)*randn(size(a_tx));
     Channel = 2*a_rx/N0;
    
    disp(['codeword: ' mat2str(codeword)]);
    disp(['received LLRs: ' mat2str(Channel,2)]);
    
    %% decision on unprocessed codeword based on LLRS 
    
    rawLLRstore(Noise_lvls,l)=sum(abs(Channel));
    zhat = zeros(1,num_v);
    zhat(Channel<0) = 1;

    %% Load LLRs into Parity Check
    H(H==0) = NaN;
    for c = 1:num_c
        for v = 1:num_v
            if matrix(c,v)==1
                H(c,v) = Channel(v);
            end
        end
    end
    
  
    %% start iterative decoding
    H(isnan(H))=0;
    clear MI_v
    clear mutual_information_c
    
    for n = 1:num_iterations
%         tic
            HplusC = [Channel;H];
    sums = sum(HplusC);% was HplusC
    
    %% decision on codeword based on LLRS
    xhat = zeros(1,num_v);
    xhat(sums<0) = 1;
    s = (mod(xhat*matrix.',2));
    syndrome_output = sum(s);
      if syndrome_output==0
      break
  end
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
        
        P0_c = exp(H)./(1+exp(H));
            P1_c = 1-P0_c;
            entropies_c = -P0_c.*log2(P0_c)-P1_c.*log2(P1_c);
            mutual_information_c(n) = 1-sum(entropies_c(find (matrix > 0)))/nnz(matrix);
        %% variable node calculations
       
        for v = 1:num_v
            int_LLRs = Channel(v);
            for c = 1:num_c
                if matrix(c,v)==1
                    int_LLRs = [int_LLRs H(c,v)];
                end
            end
           % int_LLRs
            pri_LLRs = variable_node(int_LLRs);
            
           % pri_LLRs
            pri_LLRs(1) = []; %remove channel value before re-writing to H
           
            for c = 1:num_c
                
                if matrix(c,v)==1
                    H(c,v) = pri_LLRs(1);
                    pri_LLRs(1) = [];
                  
                end
            end
        end
        %disp(['end of iteration ',num2str(n)]);
%         toc
        H(isnan(H)) = 0;
        
        P0_v = exp(nonzeros(H))./(1+exp(nonzeros(H)));
            P1_v = 1-P0_v;
            entropies_v = -P0_v.*log2(P0_v)-P1_v.*log2(P1_v);
            MI_v(n) = 1-sum(sum(entropies_v(~isnan(entropies_v)),1))/numel(entropies_v);
  HplusC = [Channel;H];
    sums = sum(HplusC);
    end
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
    sums = sum(HplusC);% was HplusC
       xhat = zeros(1,num_v);
    xhat(sums<0) = 1;
    s = (mod(xhat*matrix.',2));
    syndrome_output = sum(s);
   
    %% what is the Mutual information at the end of the iteration
    
    
    P0 = exp(Channel)./(1+exp(Channel));
    P1 = 1-P0;
    entropies = -P0.*log2(P0)-P1.*log2(P1);
    MI_Channel(Noise_lvls,l) = 1-sum(entropies(~isnan(entropies)))/numel(entropies);
      
    %% display LLRS
      disp(sprintf('LLRs after %d iterations: %s',n,mat2str(sums,2)));
    processedLLRstore(Noise_lvls,l)=sum(abs(sums));
    
    %% display EXIT
    if n>=1
    if (mod(xhat*matrix.',2)) == 0
    MI_EXIT_V = [MI_Channel(Noise_lvls,l) MI_v 1];
    MI_EXIT_C = [0 mutual_information_c 1];
    else
        MI_EXIT_V = [MI_Channel(Noise_lvls,l) MI_v];
    MI_EXIT_C = [0 mutual_information_c];
    end
    figure (1)
    cla
    hold off
    axis equal
    grid minor
    xlim([0 1])
    ylim([0 1])
    title('Full EXIT Function')
    xlabel('I_A')
    ylabel('I_E')
    hold on
    plot (MI_EXIT_C,MI_EXIT_V,'d--'),stairs(MI_EXIT_C,MI_EXIT_V,'-')
    
    %find equivaelnt CND EXIT function
    [xCND,yCND]=stairs(MI_EXIT_C,MI_EXIT_V);
    MI_CNDx=xCND(2:2:end);
     MI_CNDy=yCND(2:2:end);
     plot(MI_CNDx,MI_CNDy,'+-.')
     
     drawnow
    
    end
    %% Calculate Syndrome
     
    s = (mod(xhat*matrix.',2));
    if (sum(s) == 0)
    %[~,indx] = ismember(xhat,codewords,'rows');
    %if (indx > 0)
        disp(['Valid codeword found: ' mat2str(xhat) sprintf('\n')]);
        correctcount = correctcount + 1;
        if sum(xhat-codeword)~=0
            falsecodeword=falsecodeword+1;
        end
         if sum(zhat-xhat)==0
            LDPCuseless=LDPCuseless+1;
        end
    else
        disp(['WRONG CODEWORD: ' mat2str(xhat) sprintf('\n')]);
    


    end
    correctcount
    
    %% Gather Data at end of iterations
    rawErrors(Noise_lvls,l)=sum(abs(zhat-codeword));
    Errors(Noise_lvls,l)=sum(abs(xhat-codeword));
    Syndrome(Noise_lvls,l)=sum(s);
    
    end
    BLER(Noise_lvls)=1-(correctcount/100);
end
figure

scale = minSNR:SNRstepsize:maxSNR;
rerr=sum(rawErrors,2);
rBER=rerr/(r*l);
semilogy(scale,rBER)
title('Raw Bit Error Rate in AWGN Channel')
xlabel('SNR')
ylabel('Raw Bit Error Rate')



figure

err=sum(Errors,2);
BER=err/(r*l);
semilogy(scale,BER)
title('Bit Error Rate for AWGN Channel after LDPC')
xlabel('SNR')
ylabel('Bit Error Rate')

figure

semilogy(scale,BER)
title('Bit Error Rate for AWGN Channel after LDPC')
xlabel('SNR')
ylabel('Bit Error Rate')


figure

semilogy(scale,BLER)
title('Block Error Rate for AWGN Channel')
xlabel('SNR')
ylabel('Syndrome Error Rate')


figure

plot(Errors,Syndrome,'.')
title('Correlation between Bit Errors and Syndrome Errors in AWGN channel w/BPSK')
xlabel('Bit Error Rate')
ylabel('Syndrome Error Rate')

figure

plot(BER,SER,'.')
title('Correlation between Bit Errors and Syndrome Error Rates in AWGN channel w/BPSK (averaged)')
xlabel('Bit Error Rate')
ylabel('Syndrome Error Rate')

figure

plot(scale,MI_c)
title('Mutual Information of LLRs post LDPC')
xlabel('SNR')
ylabel('Final MI')


figure

plot(scale,MI_Channel)
title('Mutual Information of LLRs at output of Modulator')
xlabel('SNR')
ylabel('Intial MI')


figure

plot(MI_Channel,MI_c,'.')
title('Relationship between MIs of Channel LLRs and post-LDPC LLRs')
xlabel('Channel MI')
ylabel('Output MI')
