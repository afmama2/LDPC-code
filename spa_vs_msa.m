
matrix = genmatrix('rb_bravo.txt');
[G,H] = H2G(matrix);
%H=matrix;%
matrix=H;



num_v = size(H,2);
num_c = size(H,1);
num_iterations = 100;
num_runs = 10;

count_s_correct = 0;
count_s_wrong = 0;
count_m_correct = 0;
count_m_wrong = 0;
count_w_correct = 0;
count_w_wrong = 0;
count_fm_correct = 0;
count_fm_wrong = 0;
count_scm_correct = 0;
count_scm_wrong = 0;
count_svm_correct = 0;
count_svm_wrong = 0;

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
Ber3 = Ber2;
Raw  = Ber2;
Ber4 = Ber;
Ber5 = Ber;
Ber6 = Ber;
[r,k] = size(G);

for Noise_lvls=1:num_points
    ber2 = zeros(1,num_runs);
    raw = ber2;
    ber3 = ber2;
    ber = ber2;
    ber4 = ber;
    ber5 = ber;
    ber6 = ber;
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
    H2 = zeros(num_c,num_v);
    Hcv = H2;
    Hvc = H2;
    H3 = H2;
    H4 = H2;
    H5 = H2;
    H6cv = H2;
    H6vc = H2;
    % H3v = H;
    matrix1 = matrix;
    %% start iterative decoding
    
    for n = 1:num_iterations
%         tic
        %% variable node calculations
        
        mat1 = [Channel;H2];
        mat2 = repmat(sum(mat1),num_c,1);
        mat2(matrix==0) = 0;
        H2 = mat2 - H2;
        %% check node calculations
        for c = 1:num_c
            H2(c,Rows{c}) = check_node(H2(c,Rows{c}));
        end

    HplusC = [Channel;H2];
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
    

    HplusC = [Channel;H2];
    sums = sum(HplusC);% was HplusC
    xhat2 = zeros(1,num_v);
    xhat2(sums<0) = 1;
    err2 = sum(abs(xhat2-codeword));
    ber2(l) = err2;
    syndrome_output2 = sum(mod(xhat2*matrix.',2));
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
    raw(l) = errr;
    
    
    for n = 1:num_iterations
%         tic
        for c = 1:num_c
            for i=Rows{c}
                H3(c,i) = Channel(i) + sum(H3(Columns{i},i)) - H3(c,i);
            end
            H3(c,Rows{c}) = check_node(H3(c,Rows{c}));
        end
    HplusC = [Channel;H3];
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
    
    
    H3plusC = [Channel;H3];
    sums = sum(H3plusC);% was HplusC
    xhat3 = zeros(1,num_v);
    xhat3(sums<0) = 1;
    err3 = sum(abs(xhat3-codeword));
    ber3(l) = err3;
    
    syndrome_output3 = sum(mod(xhat3*matrix.',2));
    
    if syndrome_output3 == 0
        if err3 == 0
            count_m_correct = count_m_correct + 1;
        else
            count_m_wrong = count_m_wrong + 1;
        end
    end
    mat2 = repmat(Channel,num_c,1);
    mat2(matrix==0) = 0;
    Hvc = mat2;
    
    for n = 1:num_iterations
%         tic
        for v = 1:num_v
            for i=Columns{v}
                Hcv(i,v) = 2*atanh(exp(sum(log(tanh(0.5*Hvc(i,Rows{i})))) - log(tanh(0.5*Hvc(i,v)))));
            end
            Hvc(Columns{v},v) = Channel(v) + sum(Hcv(Columns{v},v)) - Hcv(Columns{v},v);
        end
    HplusC = [Channel;Hcv];
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
    
    HplusC = [Channel;Hcv];
    sums = sum(HplusC);% was HplusC
    xhat = zeros(1,num_v);
    xhat(sums<0) = 1;
    err = sum(abs(xhat-codeword));
    ber(l) = err;
    syndrome_output = sum(mod(xhat*matrix.',2));
    if syndrome_output == 0
        if err == 0
            count_w_correct = count_w_correct + 1;
        else
            count_w_wrong = count_w_wrong + 1;
        end
    end
    


    for n = 1:num_iterations
%         tic
        %% variable node calculations
        
        mat1 = [Channel;H4];
        mat2 = repmat(sum(mat1),num_c,1);
        mat2(matrix==0) = 0;
        H4 = mat2 - H4;
        %% check node calculations
        for c = 1:num_c
            H4(c,Rows{c}) = check_node1(H4(c,Rows{c}));
        end

    HplusC = [Channel;H4];
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
    

    HplusC = [Channel;H4];
    sums = sum(HplusC);% was HplusC
    xhat4 = zeros(1,num_v);
    xhat4(sums<0) = 1;
    err4 = sum(abs(xhat4-codeword));
    ber4(l) = err4;
    syndrome_output4 = sum(mod(xhat4*matrix.',2));
    if syndrome_output4 == 0
        if err4 == 0
            count_fm_correct = count_fm_correct + 1;
        else
            count_fm_wrong = count_fm_wrong + 1;
        end
    end
    
    for n = 1:num_iterations
%         tic
        for c = 1:num_c
            for i=Rows{c}
                H5(c,i) = Channel(i) + sum(H5(Columns{i},i)) - H5(c,i);
            end
            H5(c,Rows{c}) = check_node1(H5(c,Rows{c}));
        end
    HplusC = [Channel;H5];
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
    
    H3plusC = [Channel;H5];
    sums = sum(H3plusC);% was HplusC
    xhat5 = zeros(1,num_v);
    xhat5(sums<0) = 1;
    err5 = sum(abs(xhat5-codeword));
    ber5(l) = err5;
    
    syndrome_output5 = sum(mod(xhat5*matrix.',2));
    
    if syndrome_output5 == 0
        if err5 == 0
            count_scm_correct = count_scm_correct + 1;
        else
            count_scm_wrong = count_scm_wrong + 1;
        end
    end
    mat2 = repmat(Channel,num_c,1);
    mat2(matrix==0) = 0;
    H6vc = mat2;
    
    for n = 1:num_iterations
%         tic
        for v = 1:num_v
            for i=Columns{v}
                brrr = check_node1(H6vc(i,Rows{i}));
                H6cv(i,v) = brrr(Rows{i}==v);
            end
            H6vc(Columns{v},v) = Channel(v) + sum(H6cv(Columns{v},v)) - H6cv(Columns{v},v);
        end
    HplusC = [Channel;H6cv];
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

    HplusC = [Channel;H6cv];
    sums = sum(HplusC);% was HplusC
    xhat6 = zeros(1,num_v);
    xhat6(sums<0) = 1;
    err6 = sum(abs(xhat6-codeword));
    ber6(l) = err6;
    syndrome_output6 = sum(mod(xhat6*matrix.',2));
    if syndrome_output6 == 0
        if err6 == 0
            count_svm_correct = count_svm_correct + 1;
        else
            count_svm_wrong = count_svm_wrong + 1;
        end
    end

    end
    Ber2(Noise_lvls) = sum(ber2)/(num_runs*k);
    Ber3(Noise_lvls) = sum(ber3)/(num_runs*k);
    Ber(Noise_lvls) = sum(ber)/(num_runs*k);
    Raw(Noise_lvls) = sum(raw)/(num_runs*k);
    Ber4(Noise_lvls) = sum(ber4)/(num_runs*k);
    Ber5(Noise_lvls) = sum(ber5)/(num_runs*k);
    Ber6(Noise_lvls) = sum(ber6)/(num_runs*k);
    disp(Noise_lvls)
end


zerolocs=(find(Ber3==0));
Ber3(zerolocs)=1e-10;


zerolocs=(find(Ber2==0));
Ber2(zerolocs)=1e-10;

zerolocs=(find(Ber==0));
Ber(zerolocs)=1e-10;

zerolocs=(find(Ber4==0));
Ber4(zerolocs)=1e-10;

zerolocs=(find(Ber5==0));
Ber5(zerolocs)=1e-10;

zerolocs=(find(Ber6==0));
Ber6(zerolocs)=1e-10;

disp(["Correct count flooding SPA:" count_s_correct]);
disp(["Wrong count flooding SPA:" count_s_wrong]);
disp(["Wrong count serial_c SPA:" count_m_wrong]);
disp(["Correct count serial_c SPA:" count_m_correct]);
disp(["Wrong count serial_v SPA:" count_w_wrong]);
disp(["Correct count serial_v SPA:" count_w_correct]);
disp(["Correct count flooding MSA:" count_fm_correct]);
disp(["Wrong count flooding MSA:" count_fm_wrong]);
disp(["Correct count serial_c MSA:" count_scm_correct]);
disp(["Wrong count serial_c MSA:" count_scm_wrong]);
disp(["Correct count serial_v MSA:" count_svm_correct]);
disp(["Wrong count serial_v MSA:" count_svm_wrong]);

figure(1);
semilogy(minSNR:SNRstepsize:maxSNR,Ber2,"-diamond");
title(sprintf("BER after %d iteration",num_iterations))
hold on

semilogy(minSNR:SNRstepsize:maxSNR,Ber3,"-x");
semilogy(minSNR:SNRstepsize:maxSNR,Ber,"-+");
semilogy(minSNR:SNRstepsize:maxSNR,Ber4,"--");
semilogy(minSNR:SNRstepsize:maxSNR,Ber5,"--");
semilogy(minSNR:SNRstepsize:maxSNR,Ber6,"--");
semilogy(minSNR:SNRstepsize:maxSNR,Raw,"-o");
legend("Flooding SPA", "Serial_c SPA", "Serial_v SPA","Flooding MSA","Serial_c MSA","Serial_v MSA","Raw");
hold off
