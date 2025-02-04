
matrix = genmatrix('rb_bravo.txt');
[G,H] = H2G(matrix);
%H=matrix;%
G1 = Systematic(G);
H1 = [G1(:,718:end).' eye(263)];
matrix=H;
H2= H1;
%
int_LLRs=zeros(size(H));

correctcount = 0;
falsecodeword=0;
LDPCuseless=0;
num_v = size(H,2);
num_c = size(H,1);
num_iterations = 50;
syndrome_output=10;
num_runs=1000;

% set scaling
minSNR=2;
SNRstepsize=0.05;
num_points=20;
maxSNR=minSNR+SNRstepsize*(num_points-1);
rawLLRstore = zeros([num_points,num_runs]);
rawErrors = zeros([num_points,num_runs]);
bfpp = zeros([num_points,num_runs]);
Errors = zeros([num_points,num_runs]);
Syndrome = zeros([num_points,num_runs]);
MI = zeros([num_points,num_runs]);
MI_Channel = zeros([num_points,num_runs]);

[r,k] = size(G);
for Noise_lvls=1:num_points
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
    
    disp(['codeword: ' mat2str(codeword)]);
    disp(['received LLRs: ' mat2str(Channel,2)]);
    
    %% decision on unprocessed codeword based on LLRS 
    
    rawLLRstore(Noise_lvls,l)=sum(abs(Channel));
    zhat = zeros(1,num_v);
    zhat(Channel<0) = 1;

    %% Load LLRs into Parity Check

    H = zeros(num_c,num_v);
  
    %% start iterative decoding

    for n = 1:num_iterations
%         tic
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
        for c = 1:num_c
            int_LLrs = [];
            for v = 1:num_v
                if matrix(c,v)==1
                    calc = Channel(v) + sum(H(:,v))-H(c,v);
                    int_LLrs = [int_LLrs calc];                
                end
            end                
            pri_LLRs = check_node(int_LLrs);
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
        P0_c = exp(H)./(1+exp(H));
            P1_c = 1-P0_c;
            entropies_c = -P0_c.*log2(P0_c)-P1_c.*log2(P1_c);
            mutual_information_c(n) = 1-sum(entropies_c(~isnan(entropies_c)))/numel(entropies_c);
        
       
        
        P0_v = exp(nonzeros(H))./(1+exp(nonzeros(H)));
            P1_v = 1-P0_v;
            entropies_v = -P0_v.*log2(P0_v)-P1_v.*log2(P1_v);
            mutual_information_v(n) = 1-sum(sum(entropies_v(~isnan(entropies_v)),1))/numel(entropies_v);

    end

    
    HplusC = [Channel;H];
    sums = sum(HplusC);% was HplusC
    xhat = zeros(1,num_v);
    xhat(sums<0) = 1;
    s = (mod(xhat*matrix.',2));
    syndrome_output = sum(s);
    %% display LLRS
    disp(sprintf('LLRs after %d iterations: %s',n,mat2str(sums,2)));
    processedLLRstore(Noise_lvls,l)=sum(abs(sums));
    
    
    
    %% Calculate Syndrome
     
    
    if (sum(s) == 0)

        disp(['Valid codeword found: ' mat2str(xhat) sprintf('\n')]);
        
        if sum(xhat-codeword)==0
            correctcount = correctcount + 1;
        
        else 
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
    
    %% what is the Mutual information at the end of the iteration & list decoding
    P0 = exp(sums)./(1+exp(sums));
    P1 = 1-P0;
    entropies = -P0.*log2(P0)-P1.*log2(P1);
    MI(Noise_lvls,l) = 1-sum(entropies(~isnan(entropies)))/numel(entropies);
    
    P0 = exp(Channel)./(1+exp(Channel));
    P1 = 1-P0;
    entropies = -P0.*log2(P0)-P1.*log2(P1);
    MI_Channel(Noise_lvls,l) = 1-sum(entropies(~isnan(entropies)))/numel(entropies);
    
    end
end
figure

scale = minSNR:SNRstepsize:maxSNR;
rerr=sum(rawErrors,2);
rBER=rerr/(r*l);
zerolocs=(find(rBER==0));
rBER(zerolocs)=1e-10;
semilogy(scale,rBER)
title('Raw Bit Error Rate in AWGN Channel')
xlabel('SNR')
ylabel('Raw Bit Error Rate')



figure

err=sum(Errors,2);
BERflood_RBAlpha=err/(r*l);
zerolocs=(find(BERflood_RBAlpha==0));
BERflood_RBAlpha(zerolocs)=1e-10;
semilogy(scale,BERflood_RBAlpha)
title('Bit Error Rate for AWGN Channel after LDPC')
xlabel('SNR')
ylabel('Bit Error Rate')

figure

semilogy(scale,BERflood_RBAlpha)
hold on
semilogy(scale,rBER)
title('Bit Error Rate for AWGN Channel comparison')
xlabel('SNR')
ylabel('Bit Error Rate')
legend("Mix","raw")
hold off

figure

syn=sum(Syndrome,2);
SER=syn/(r*l);
zerolocs=(find(SER==0));
SER(zerolocs)=1e-10;
semilogy(scale,SER)
title('Syndrome Error Rate for AWGN Channel')
xlabel('SNR')
ylabel('Syndrome Error Rate')


figure

plot(Errors,Syndrome,'.')
title('Correlation between Bit Errors and Syndrome Errors in AWGN channel w/BPSK')
xlabel('Bit Error Rate')
ylabel('Syndrome Error Rate')

figure

plot(BERflood_RBAlpha,SER,'.')
title('Correlation between Bit Errors and Syndrome Error Rates in AWGN channel w/BPSK (averaged)')
xlabel('Bit Error Rate')
ylabel('Syndrome Error Rate')

figure

plot(scale,MI)
title('Mutual Information of LLRs post LDPC')
xlabel('SNR')
ylabel('Final MI')


figure

plot(scale,MI_Channel)
title('Mutual Information of LLRs at output of Modulator')
xlabel('SNR')
ylabel('Intial MI')


figure

plot(MI_Channel,MI,'.')
title('Relationship between MIs of Channel LLRs and post-LDPC LLRs')
xlabel('Channel MI')
ylabel('Output MI')
