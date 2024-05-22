
matrix = genmatrix('rb_bravo.txt');
[G,H] = H2G(matrix);
%H=matrix;%
matrix=H;



num_v = size(H,2);
num_c = size(H,1);
num_iterations = 100;
num_iterations1 = 1000;
num_iterations2 = 2000;
num_runs = 10;

count_s_correct = 0;
count_s_wrong = 0;
count_f_correct = 0;
count_f_wrong = 0;

Rows = cell(num_c);
Columns = cell(num_v);

for n=1:num_c
    row = [];
    for m=1:num_v
        if matrix(n,m)==1
            row = [row m];
        end
    end
    Rows{n} = row;
end

for n=1:num_v
    column = [];
    for m=1:num_c
        if matrix(m,n)==1
            column = [column m];
        end
    end
    Columns{n} = column;
end

% set scaling
minSNR=1;
SNRstepsize=0.05;
num_points=20;
maxSNR=minSNR+SNRstepsize*(num_points-1);

Ber2 = zeros(1,num_points);
Ber = Ber2;
Raw  = Ber2;
[r,k] = size(G);


for Noise_lvls=1:num_points
    ber2 = zeros(1,num_runs);
    ber = ber2;
    raw = ber2;
    parfor l=1:num_runs
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
    H2 = zeros(num_c,num_v);
    H = H2;
    % H3v = H;
    matrix1 = matrix;
    %% start iterative decoding
    for n = 1:num_iterations
            for c = 1:num_c
                for i=Rows{c}
                    H(c,i) = Channel(i) + sum(H(Columns{i},i)) - H(c,i);
                end
                H(c,Rows{c}) = check_node(H(c,Rows{c}));
            end
    HplusC = [Channel;H];
    sums = sum(HplusC);





    %% decision on codeword based on LLRS
    xhat = zeros(1,num_v);
    xhat(sums<0) = 1;
    s = (mod(xhat*matrix1.',2));
    syndrome_output = sum(s);
      if syndrome_output==0
      break
      end

    end
    

    HplusC = [Channel;H];
    sums = sum(HplusC);% was HplusC
    xhat = zeros(1,num_v);
    xhat(sums<0) = 1;
    err = sum(abs(xhat-codeword));
    ber(l) = err;
    syndrome_output = sum(mod(xhat*matrix1.',2));
    if syndrome_output == 0
        if err == 0
            count_f_correct = count_f_correct + 1;
        else
            count_f_wrong = count_f_wrong + 1;
        end
    end
    for n = 1:num_iterations1
            for c = 1:num_c
                for i=Rows{c}
                    H2(c,i) = Channel(i) + sum(H2(Columns{i},i));
                end
                H2(c,Rows{c}) = check_node(H2(c,Rows{c}));
            end
    HplusC = [Channel;H2];
    sums = sum(HplusC);





    %% decision on codeword based on LLRS
    xhat = zeros(1,num_v);
    xhat(sums<0) = 1;
    s = (mod(xhat*matrix1.',2));
    syndrome_output = sum(s);
      if syndrome_output==0
      break
      end

    end
    

    HplusC = [Channel;H2];
    sums = sum(HplusC);% was HplusC
    xhat2 = zeros(1,num_v);
    xhat2(sums<0) = 1;
    err2 = sum(abs(xhat2-codeword));
    ber2(l) = err2;
    syndrome_output2 = sum(mod(xhat2*matrix1.',2));
    if syndrome_output2 == 0
        if err2 == 0
            count_s_correct = count_s_correct + 1;
        else
            count_s_wrong = count_s_wrong + 1;
        end
    end

    xhatr = zeros(1,num_v);
    xhatr(Channel<0) = 1;
    errr = sum(abs(xhatr-codeword));

    end
    Ber2(Noise_lvls) = sum(ber2)/(num_runs*k);
    Ber(Noise_lvls) = sum(ber)/(num_runs*k);
    Raw(Noise_lvls) = sum(raw)/(num_runs*k);
    ber2(ber2==0) = 500;
    disp(min(ber2))
end

zerolocs=(find(Ber==0));
Ber(zerolocs)=1e-10;


zerolocs=(find(Ber2==0));
Ber2(zerolocs)=1e-10;

disp(["Correct count after 100:" count_s_correct]);
disp(["Wrong count after 100:" count_s_wrong]);
disp(["Correct count serial 100:" count_f_correct]);
disp(["Wrong count serial 100:" count_f_wrong]);

figure(1);
semilogy(minSNR:SNRstepsize:maxSNR,Ber2,"-diamond");
title(sprintf("BER after %d iteration",num_iterations))
hold on
semilogy(minSNR:SNRstepsize:maxSNR,Ber,"-o");
semilogy(minSNR:SNRstepsize:maxSNR,Raw,"-|");
legend("Serial_c SFA","Serial_c","Raw");
hold off
