
matrix = genmatrix('rb_bravo.txt');
[G,H] = H2G(matrix);
%H=matrix;%
matrix=H;


num_v = size(H,2);
num_c = size(H,1);
num_iterations = 300;
num_runs=10000;

sdiff = 3;
fdiff = 3;

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

Ber = zeros(num_points,num_iterations/3);
Ber2 = zeros(num_points,num_iterations/3);
Ber3 = zeros(num_points,num_iterations);
Raw  = Ber3;
Ber6 = Ber3;
[r,k] = size(G);

parfor Noise_lvls=1:num_points
    ber = zeros(num_runs,num_iterations/3);
    ber2 = zeros(num_runs,num_iterations/3);
    ber3 = zeros(num_runs,num_iterations);
    raw = ber3;
    ber6 = ber3;
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
   
    count = 0;
    syndrome_output = 1;
    syndrome_output2 = 1;
    syndrome_output3 = 1;
    syndrome_output4 = 1;


    xhatr = zeros(1,num_v);
    xhatr(Channel<0) = 1;
    errr = sum(abs(xhatr-codeword));
    raw(l,:) = errr*ones(1,num_iterations);
    for n = 1:num_iterations
%         tic
        count = count + 1;


        %% Flooding_________________________________

        if syndrome_output~=0 && rem(n,3)==0
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
        end
        %% Serial________________________________________________

        if syndrome_output2~=0 && rem(n,3) == 0
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
   
        end
        if rem(n,3)==0
    HplusC = [Channel;H];
    sums = sum(HplusC);% was HplusC
    xhat = zeros(1,num_v);
    xhat(sums<0) = 1;
    err = sum(abs(xhat-codeword));
    ber(l,n/3) = err;

    HplusC = [Channel;H2];
    sums = sum(HplusC);% was HplusC
    xhat2 = zeros(1,num_v);
    xhat2(sums<0) = 1;
    err2 = sum(abs(xhat2-codeword));
    ber2(l,n/3) = err2;
        end

   
        %% Mix______________________________________
        if syndrome_output3~=0 
        mat1 = [Channel;H3];
        mat2 = repmat(sum(mat1),num_c,1);
        mat2(matrix==0) = 0;
        H3 = mat2 - H3;
        %disp(['end of iteration ',num2str(n)]);
%         toc
            
            %% check node calculations
        for c = 1:num_c
            H3(c,Rows{c}) = check_node1(H3(c,Rows{c}));
        end
        H3plusC = [Channel;H3];
        sums3 = sum(H3plusC);
        xhat3 = zeros(1,num_v);
        xhat3(sums3<0) = 1;

        s = (mod(xhat3*matrix.',2));
        syndrome_output3 = sum(s);   
       
        end
   
   
        if syndrome_output4~=0
        for c = 1:num_c
            for i=Rows{c}
                H4(c,i) = Channel(i) + sum(H4(Columns{i},i)) - H4(c,i);
            end
            H4(c,Rows{c}) = check_node1(H4(c,Rows{c}));
        end
        
        H4plusC = [Channel;H4];
        sums4 = sum(H4plusC);
        xhat4 = zeros(1,num_v);
        xhat4(sums4<0) = 1;

        s = (mod(xhat4*matrix.',2));
        syndrome_output4 = sum(s);   
        end
   
    H3plusC = [Channel;H3];
    sums = sum(H3plusC);% was HplusC
    xhat3 = zeros(1,num_v);
    xhat3(sums<0) = 1;
    err3 = sum(abs(xhat3-codeword));
    ber3(l,n) = err3;


    
    H4plusC = [Channel;H4];
    sums = sum(H4plusC);% was HplusC
    xhat4 = zeros(1,num_v);
    xhat4(sums<0) = 1;
    err4 = sum(abs(xhat4-codeword));
    ber6(l,n) = sum(abs(xhat4-codeword));
    
    
    end
    end
    Ber(Noise_lvls,:) = sum(ber)/(num_runs*k);
    Ber2(Noise_lvls,:) = sum(ber2)/(num_runs*k);
    Raw(Noise_lvls,:) = sum(raw)/(num_runs*k);
    Ber3(Noise_lvls,:) = sum(ber3)/(num_runs*k);
    Ber6(Noise_lvls,:) = sum(ber6)/(num_runs*k);
    disp(Noise_lvls)
end


zerolocs=(find(Ber==0));
Ber(zerolocs)=1e-10;


zerolocs=(find(Ber2==0));
Ber2(zerolocs)=1e-10;

zerolocs=(find(Ber3==0));
Ber3(zerolocs)=1e-10;

zerolocs=(find(Ber6==0));
Ber6(zerolocs)=1e-10;



[X, Y] = meshgrid(1:num_iterations,minSNR:SNRstepsize:maxSNR);
[X1, Y1] = meshgrid(3:3:num_iterations,minSNR:SNRstepsize:maxSNR);
figure(1);

surf(X1,Y1,Ber,'FaceColor',[1 0 0]);
set(gca,'zscale','log')
title("BER for different schedulings")
hold on
surf(X1,Y1,Ber2,'FaceColor',[0 1 0]);
set(gca,'zscale','log')
surf(X,Y,Raw,'FaceColor',[0 0 1]);
set(gca,'zscale','log')
surf(X,Y,Ber3,'FaceColor',[0.5 0.5 0.5]);
set(gca,'zscale','log')
surf(X,Y,Ber6,'FaceColor',[1 0 1]);
set(gca,'zscale','log')
legend("Flooding SPA","Serial_v SPA","Raw","Flooding MSA","Serial_v MSA");
hold off
