
matrix = genmatrix('rb_bravo.txt');
[G,H] = H2G(matrix);
%H=matrix;%
matrix=H;


num_v = size(H,2);
num_c = size(H,1);
num_iterations = 1000;
num_runs=1000;

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
num_points=20;
maxSNR=minSNR+SNRstepsize*(num_points-1);

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


Ber = zeros(1,num_points);
Ber2 = zeros(1,num_points);
Ber3 = Ber;
Raw  = Ber;
Ber4 = Ber;
Ber6 = Ber;
Ber7 = Ber;
[r,k] = size(G);

parfor Noise_lvls=1:num_points
    ber = zeros(1,num_runs);
    ber2 = zeros(1,num_runs);
    raw = ber;
    ber3 = ber;
    ber4 = ber;
    ber6 = ber;
    ber7 = ber;
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
    H4 = H;
    % H3v = H;
  
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
            H(c,Rows{c}) = check_node(H(c,Rows{c}));
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
            for i=Rows{c}
                H2(c,i) = Channel(i) + sum(H2(Columns{i},i)) - H2(c,i);
            end
            H2(c,Rows{c}) = check_node(H2(c,Rows{c}));
        end
        HplusC = [Channel;H2];
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
    
    HplusC = [Channel;H];
    sums = sum(HplusC);% was HplusC
    xhat = zeros(1,num_v);
    xhat(sums<0) = 1;
    err = sum(abs(xhat-codeword));
    ber(l) = err;
    syndrome_output = sum(mod(xhat*matrix.',2));
    if syndrome_output == 0
        if err == 0
            count_f_correct = count_f_correct + 1;
        else
            count_f_wrong = count_f_wrong + 1;
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

    for n = 1:num_iterations/2

        mat1 = [Channel;H3];
        mat2 = repmat(sum(mat1),num_c,1);
        mat2(matrix==0) = 0;
        H3 = mat2 - H3;
        %disp(['end of iteration ',num2str(n)]);
%         toc
            
            %% check node calculations
        for c = 1:num_c
            H3(c,Rows{c}) = check_node(H3(c,Rows{c}));
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
            for i=Rows{c}
                H3(c,i) = Channel(i) + sum(H3(Columns{i},i)) - H3(c,i);
            end
            H3(c,Rows{c}) = check_node(H3(c,Rows{c}));
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
    
    xhat4 = xhat + xhat2 + xhat3;
    xhat4(xhat4<2) = 0;
    xhat4(xhat4>1) = 1;
    ber6(l) = sum(abs(xhat4-codeword));
    
    if syndrome_output == 0
        ber7(l) = err;
    elseif syndrome_output3 == 0
        ber7(l) = err3;
    elseif syndrome_output2 == 0
        ber7(l) = err2;
    else
        ber7(l) = sum(abs(xhat4-codeword));
    end
    

    for n = 1:(3*num_iterations)

        for c = 1:num_c
            for i=Rows{c}
                H4(c,i) = Channel(i) + sum(H4(Columns{i},i)) - H4(c,i);
            end
            H4(c,Rows{c}) = check_node(H4(c,Rows{c}));
        end
        HplusC = [Channel;H4];
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
    HplusC = [Channel;H4];
    sums = sum(HplusC);% was HplusC
    xhat5 = zeros(1,num_v);
    xhat5(sums<0) = 1;
    err5 = sum(abs(xhat5-codeword));
    ber4(l) = err5;

    end
    Ber(Noise_lvls) = sum(ber)/(num_runs*k);
    Ber2(Noise_lvls) = sum(ber2)/(num_runs*k);
    Raw(Noise_lvls) = sum(raw)/(num_runs*k);
    Ber3(Noise_lvls) = sum(ber3)/(num_runs*k);
    Ber4(Noise_lvls) = sum(ber4)/(num_runs*k);
    Ber6(Noise_lvls) = sum(ber6)/(num_runs*k);
    Ber7(Noise_lvls) = sum(ber7)/(num_runs*k);
    disp(count)
end


zerolocs=(find(Ber==0));
Ber(zerolocs)=1e-10;


zerolocs=(find(Ber2==0));
Ber2(zerolocs)=1e-10;

zerolocs=(find(Ber3==0));
Ber3(zerolocs)=1e-10;

zerolocs=(find(Ber4==0));
Ber4(zerolocs)=1e-10;


zerolocs=(find(Ber6==0));
Ber6(zerolocs)=1e-10;

zerolocs=(find(Ber7==0));
Ber7(zerolocs)=1e-10;

disp(["Correct serial count:" count_s_correct]);
disp(["Wrong serial count:" count_s_wrong]);
disp(["Correct flooding count:" count_f_correct]);
disp(["Wrong flooding count:" count_f_wrong]);
disp(["Wrong mix count:" count_m_wrong]);
disp(["Correct mix count:" count_m_correct]);

figure(1);

semilogy(minSNR:SNRstepsize:maxSNR,Ber,"-square");
title(sprintf("BER after %d iteration",num_iterations))
hold on
semilogy(minSNR:SNRstepsize:maxSNR,Ber2,"-diamond");
semilogy(minSNR:SNRstepsize:maxSNR,Raw,"-o");
semilogy(minSNR:SNRstepsize:maxSNR,Ber3,"-x");
semilogy(minSNR:SNRstepsize:maxSNR,Ber4,"-^");
semilogy(minSNR:SNRstepsize:maxSNR,Ber6,"-+");
semilogy(minSNR:SNRstepsize:maxSNR,Ber7,"-|");
legend("Flooding","Serial","Raw","Mix","3x Serial","Most likely bit","Best of four");
hold off
