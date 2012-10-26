% Assumes the .mat file Ca2SignalsUdatedAfterViva.mat 
% exists in the current directory;
% The file Ca2SignalsUdatedAfterViva.mat is
% Produced by S. Iftikhar and must contain the
% following fields:
%   <TO DO>
%
% 
% AAB Comment - this script is now trusted; 
% this in preference to Scripts1-4.  
% Current limitation: does not study neighbourhood
% statistics
%
% TO DO: Include studies using Poisson stats
% for the final P() calcluation
% OR replace Poisson stats with Gamma stats for
% the first plot.
%
% AAB Updated 13/07/2011 - Added lines for "AFTER FLOW"
% AAB Updated 12/07/2011

% Close all open windows
close all

load Ca2SignalsUdatedAfterViva.mat;
fps = 1/5;

timeaxis = [0:1:143]/fps;

% Get Centroids for ease of processing
 Cr = cat(1,Cells(:).Centroid1);
 Cc = cat(1,Cells(:).Centroid2);
 
 figure(1); % Place the slide cartoon into Figure 1
 plot(Cr,Cc,'o','MarkerFaceColor','b');
 axis ij; % Mode compat with images
 axis off;
 for n = 1:length(Cr) % number of cells
     text(Cr(n)-5,Cc(n)+10,num2str(n));
 end
 
hf = figure(2); % Place traces in Figure 2
set(hf,'Name','During Flow');
MatrixOfPeaksDuring = zeros(50,144); 
for c = 1:50
    t = SpikeResultsDur(c).locs+1; % t: Array
    MatrixOfPeaksDuring(c,t) = 1; 
end
imagesc(MatrixOfPeaksDuring);colormap(gray(2));
xlabel('Time (s)');ylabel('Cell Label');
title('During Flow');
OriginalMatrixOfPeaksDuring = MatrixOfPeaksDuring;

 hf = figure(3); % Place traces in Figure 3
 set(hf,'Name','After Flow Cessation');
 
MatrixOfPeaksAfter = zeros(50,144); 
for c = 1:50
    t = SpikeResultsAft(c).locs+1; % t: Array
    MatrixOfPeaksAfter(c,t) = 1; 
end
imagesc(MatrixOfPeaksAfter);colormap(gray(2));
xlabel('Time (s)');ylabel('Cell Label');
title('After Flow Cessation');

OriginalMatrixOfPeaksAfter = MatrixOfPeaksAfter;

%% Simple rate stats; during and after
% Both of these rates have units s^{-1}
RateOfPeaksDuring = sum(MatrixOfPeaksDuring,2)/max(timeaxis);
RateOfPeaksAfter = sum(MatrixOfPeaksAfter,2)/max(timeaxis);
hPR = figure;

% Find range to create a single axis for the histogram
RateAxisMin = min([RateOfPeaksDuring,RateOfPeaksAfter]);
RateAxisMax = max([RateOfPeaksDuring,RateOfPeaksAfter]);

RateAxis = RateAxisMin:(RateAxisMax-RateAxisMin)/9:RateAxisMax;

subplot(1,2,1);
[N,X]=hist(RateOfPeaksDuring,RateAxis);title('During Flow');
bar(X,N);
xlabel('Rate (Peaks per Second)');
ylabel('Count');

subplot(1,2,2);
[N,X]=hist(RateOfPeaksAfter,RateAxis);title('After Flow Cessation');
bar(X,N);
xlabel('Rate (Peaks per Second)');
ylabel('Count');

%% Bootstrap by sampling the rows to get an estimate 
% of the "lambda" distribution, given Poisson rates.  
% Chose to use 30 samples of cells each time, in 
% order to be conservative, but these
% are selected, with replacement, from all 50 cells
% for i = 1:10000
%      RowSamples = randi([1,50],30,1);
%      SamplingOfPeaksDuring = MatrixOfPeaksDuring(RowSamples,:);
%      BERatesDuring(i)=mean(sum(SamplingOfPeaksDuring,2)/max(timeaxis));
% end
%      
% for i = 1:10000
%      RowSamples = randi([1,50],30,1);
%      SamplingOfPeaksAfter = MatrixOfPeaksAfter(RowSamples,:);
%      BERatesAfter(i)=mean(sum(SamplingOfPeaksAfter,2)/max(timeaxis));
% end

for i = 1:10000
     RowSamples = randi([1,50],30,1);
     SamplingOfPeaksDuring = MatrixOfPeaksDuring(RowSamples,:);
     SamplingOfPeaksAfter = MatrixOfPeaksAfter(RowSamples,:);
     BERatesDuring(i)=mean(sum(SamplingOfPeaksDuring,2)/max(timeaxis));
     BERatesAfter(i)=mean(sum(SamplingOfPeaksAfter,2)/max(timeaxis));
end




AllBERates = [BERatesDuring,BERatesAfter];
MinBERates = min(AllBERates);MaxBERates=max(AllBERates);
Rates = linspace(MinBERates,MaxBERates,80);

[CD,RD]=hist(BERatesDuring,Rates);
[CA,RA]=hist(BERatesAfter,Rates);

figure;
bar(Rates,CD/sum(CD),'b'); hold on
B=bar(Rates,CA/sum(CA),'r'); hold off
ch = get(B,'child');
set(ch,'facea',.3);
title('Rates During Flow and After Flow Sessation');
legend('During','After');
xlabel('Parameter \lambda, s^{-1}');
ylabel('P(\lambda)');


%% The code from here onwards ONLY considers the 
% Coincident statistics DURING flow, NOT after !
% To see if there is any difference in the 
% Coincidence matrix permute each row indendently
for i = 1:50
    ColOrder = randperm(length(timeaxis));
    ReorderedMatrixOfPeaksDuring(i,:) = OriginalMatrixOfPeaksDuring(i,ColOrder);
end

% Let's now mine this matrix 
% Everything First During Flow; 
% original order

CoincidenceDuring = zeros(50,50);
for cA = 1:50
    for cB = 1:50
        CoincidenceDuring(cA,cB) = sum(MatrixOfPeaksDuring(cA,:) &...
                            MatrixOfPeaksDuring(cB,:));
    end
end

% Same Calculation, Permuted Indices

CoincidenceDuringReordered = zeros(50,50);
for cA = 1:50
    for cB = 1:50
        CoincidenceDuringReordered(cA,cB) = sum(ReorderedMatrixOfPeaksDuring(cA,:) &...
                            ReorderedMatrixOfPeaksDuring(cB,:));
    end
end
disp('');

% Calculate number of non-zero elements for the original
% co-incidence matrix

% The following would find only cells that fire
% together, and does not consider the *number* of times
% they would fire together; more conservative
% NCD=length(find(CoincidenceDuring(:)));
% RNCD=length(find(ReorderedCoincidenceDuring(:)));

% The following (uncommented), calculates the 
% number of same-interval peak co-incidences
% This includes wach cell with itself
NCD=sum(sum(triu(CoincidenceDuring))); % Upper diagonal
RNCD=sum(sum(triu(CoincidenceDuringReordered))); % Upper diag 


% Find the total number of cell peaks
% from the trace of the coincidence
TotalNPeaksDuring = trace(CoincidenceDuring);
TotalNPeaksDuringReordered = trace(CoincidenceDuringReordered);

% Sanity check - these should be equal in ordered and re-ordered sequence
if TotalNPeaksDuring ~= TotalNPeaksDuringReordered
    error('Something went wrong during permute...');
end

% The following calculation corrects the co-incidences
% by accounting for each cell peaking with itself
NCD = NCD - TotalNPeaksDuring; % Original sequence
RNCD = RNCD - TotalNPeaksDuringReordered; % Re-ordered sequence

% Finally, we need to correct fot he following effect.  Given a cel firing
% pattern like this:
% Cell 1: 1 1 0 1 0 1 1
% Cell 2: 1 1 0 1 0 1 1 
% Cell 3: 1 1 0 1 0 1 1 
% Cell 4: 1 1 0 1 0 1 1 

% The coincidence matrix will look ike this
% == C1   C2   C3   C4
% C1  5    5    5    5  
% C2  5    5    5    5
% C3  5    5    5    5
% C4  5    5    5    5

% The upper trianular part will be
% 
% The coincidence matrix will look like this
% == C1   C2   C3   C4
% C1  *    5    5    5  
% C2  0    *    5    5
% C3  0    0    *    5
% C4  0    0    0    *

% Giving a total score of 30; the max score for a given pattern
% containing K spikes will have K entries, all in the positions as shown
% above.  For a N cells (e.g. 4), the max number of entries containing K
% will be (N-1)+(N-2)+(N-3) = 3N - 6, or the sum of the first N-1 terms of
% an arithmetic series or (N-1)(1 + N-1)/2 = N(N-1)/2.  Thus, the
% normalising factor should be (K/2).N.(N-1)
%%
% This cell is a demonstration of the principle 
% of Monte Carlo resampling applied to the
% array of peak times; in particular, this shows
% how co-incidences change with a single randomisation
%  By studying how the distribution of any measures of co
% -incidence firing change over many pseudo-radom experiments
% it is possible to get robust estimates of the probability
% of accidentally obtaining such correspondences in
% independently peaking cells.
hCAF = figure; % Coincidence array figure
colormap(jet(18)); % Based on max number of known spikes
subplot(1,2,1);imagesc(CoincidenceDuring);colorbar
xlabel('Cell Label');ylabel('Cell Label');
title(['#(coincident peaks) in original: ',num2str(NCD)]);
subplot(1,2,2);imagesc(CoincidenceDuringReordered);colorbar
xlabel('Cell Label');ylabel('Cell Label');
title(['#(conincident peaks) in reordered: ',num2str(RNCD)])
set(hCAF,'Name','Demonstration of permutation statistics');
disp('');

%%
% Proper Monte-Carlo study begins now...
% To see if there is any difference in the Coincidence matrix
% Permute each row indendently
disp('Preparing for true M-C simulations...')
pause;

for trial = 1:1002 % Should be at least 10000..
    for i = 1:50
        ColOrder = randperm(length(timeaxis));
        MatrixOfPeaksDuringReordered(i,:) = OriginalMatrixOfPeaksDuring(i,ColOrder);
    end

    for cA = 1:50
        for cB = 1:50
            CoincidenceDuringReordered(cA,cB) = sum(MatrixOfPeaksDuringReordered(cA,:) &...
                            MatrixOfPeaksDuringReordered(cB,:));
        end
    end
   % RNCD = length(find(ReorderedCoincidenceDuring(:)));
    RNCD = sum(sum(triu(CoincidenceDuringReordered)));
    RNCD = RNCD - TotalNPeaksDuring;
    
    ThisTrial(trial) = RNCD;
end
disp('');

figdist = figure;
subplot(2,2,1);hist(ThisTrial,1000);
xlabel('#(Coincident Peaks)');
ylabel('#(Trials)');

% Now, for fitting and plotting the best fits
subplot(2,2,2);[N,X] = hist(ThisTrial,1000);
% Calculate DX
DX = X(2) - X(1); % Assumed equal for all X values

% Presumed empirical probabilities...
Phat1 = N/(sum(N));

% First bar plot of empirical probabilities
bar(X,Phat1); hold on
xlabel('#(Coincident Peaks), C_{obs}');
ylabel('P(C_{obs})');
%
[LHat, Lci] = poissfit(ThisTrial,0.01);
RawPDF = poisspdf(round(X),LHat);
plot(round(X),RawPDF,'--');

% Optimise Graphical fit...
% This should be analytically do-able
% but since different parametric curves are being tried
% this approach is flexible.
Aopt = RawPDF * pinv(Phat1);
plot(round(X),RawPDF/Aopt,'-');hold off
title('Best Poisson model');
[lh,oh]=legend('','PDF fit','Fit, optimised for display');
delete(oh(4));

% Next, try a Gamma Fit
subplot(2,2,3);
bar(X,Phat1); hold on
xlabel('#(Coincident Peaks), C_{obs}');
ylabel('P(C_{obs})');

[PARMHAT,PARMCI] = gamfit(ThisTrial,0.01);
RawPDF = gampdf(X,PARMHAT(1),PARMHAT(2));
h1=plot(X,RawPDF,'r--'); 
set(h1,'LineWidth',2);

% Optimise Graphical fit...
Aopt = RawPDF * pinv(Phat1);
h2=plot(X,RawPDF/Aopt,'r-');hold off
set(h2,'LineWidth',2);
title('Best Gamma model');
[lh,oh]=legend('','PDF fit','Fit, optimised for display');
delete(oh(4));

subplot(2,2,4);
% Need to re-do X axis sensibly....
Xmin = 0; Xmax = max([max(X),ThisTrial,NCD]);
NewXaxis = linspace(Xmin,Xmax*1.1,100); % 1.1 is a fudge factor


% Plot PDF curve for these parameters, extending out to the
% 10pc beyond the biggest coincident value we have found...
ScaledPDF = gampdf(NewXaxis,PARMHAT(1),PARMHAT(2))/Aopt;
plot(NewXaxis,ScaledPDF,'b');hold on
hL=line([NCD,NCD],[0 max(ScaledPDF)]);
set(hL,'Color','r','LineStyle','--'); hold off;
set(hL,'LineWidth',2);
title('Gamma fit and P(C_{obs}|Gamma;a,b)')
PCEqNCD = gampdf(NCD,PARMHAT(1),PARMHAT(2));
format shorte;
text(NCD,max(ScaledPDF),['P=',num2str(PCEqNCD)]); 
xlabel('#(Coincident Peaks), C_{obs}');
ylabel('P(C_{obs})');
hold off;


%%  Added lines on 13/07  AFTER FLOW
% The code from here onwards ONLY considers the 
% Coincident statistics AFTER flow, NOT BEFORE !
% To see if there is any difference in the 
% Coincidence matrix permute each row indendently
for i = 1:50
    ColOrder = randperm(length(timeaxis));
    ReorderedMatrixOfPeaksAfter(i,:) = OriginalMatrixOfPeaksAfter(i,ColOrder);
end

% Let's now mine this matrix 
% Everything First During Flow; 
% original order

for cA = 1:50
    for cB = 1:50
        CoincidenceAfter(cA,cB) = sum(MatrixOfPeaksAfter(cA,:) &...
                            MatrixOfPeaksAfter(cB,:));
    end
end

% Same Calculation, Permuted Indices

for cA = 1:50
    for cB = 1:50
        CoincidenceAfterReordered(cA,cB) = sum(ReorderedMatrixOfPeaksAfter(cA,:) &...
                            ReorderedMatrixOfPeaksAfter(cB,:));
    end
end
disp('');

% Calculate number of non-zero elements for the original
% co-incidence matrix

% The following would find only cells that fire
% together, and does not consider the *number* of times
% they would fire together; more conservative
% NCD=length(find(CoincidenceDuring(:)));
% RNCD=length(find(ReorderedCoincidenceDuring(:)));

% The following (uncommented), calculates the 
% number of same-interval peak co-incidences
% This includes wach cell with itself
NCA=sum(sum(triu(CoincidenceAfter))); % upper diag
RNCA=sum(sum(triu(CoincidenceAfterReordered)));


% Find the total number of cell peaks
% from the trace of the coincidence
TotalNPeaksAfter = trace(CoincidenceAfter);
TotalNPeaksAfterReordered = trace(CoincidenceAfterReordered);

% Sanity check - these should be equal in ordered and re-ordered sequence
if TotalNPeaksAfter ~= TotalNPeaksAfterReordered
    error('Something went wrong during permute...');
end

% The following calculation corrects the co-incidences
% by accounting for each cell peaking with itself
NCA = NCA - TotalNPeaksAfter; % Original sequence
RNCA = RNCA - TotalNPeaksAfterReordered; % Re-ordered sequence

%%
% This cell is a demonstration of the principle 
% of Monte Carlo resampling applied to the
% array of peak times; in particular, this shows
% how co-incidences change with a single randomisation
%  By studying how the distribution of any measures of co
% -incidence firing change over many pseudo-radom experiments
% it is possible to get robust estimates of the probability
% of accidentally obtaining such correspondences in
% independently peaking cells.
hCAF = figure; % Coincidence array figure
colormap(jet(18));
subplot(1,2,1);imagesc(CoincidenceAfter);
xlabel('Cell Number');ylabel('Cell Number');
title(['# Sim. peaks in original (After): ',num2str(NCA)]);
subplot(1,2,2);imagesc(CoincidenceAfterReordered);
xlabel('Cell Number');ylabel('Cell Number');
title(['# Sim. peaks in reordered: ',num2str(RNCA)])
set(hCAF,'Name','Demonstration of permutation statistics');
disp('');

%%
% Proper Monte-Carlo study begins now...
% To see if there is any difference in the Coincidence matrix
% Permute each row indendently
disp('Preparing for true M-C simulations...')
pause;

for trial = 1:1000 % Should be at least 10000..
    for i = 1:50
        ColOrder = randperm(length(timeaxis));
        MatrixOfPeaksAfterReordered(i,:) = OriginalMatrixOfPeaksAfter(i,ColOrder);
    end

    for cA = 1:50
        for cB = 1:50
            CoincidenceAfterReordered(cA,cB) = sum(MatrixOfPeaksAfterReordered(cA,:) &...
                            MatrixOfPeaksAfterReordered(cB,:));
        end
    end
   % RNCD = length(find(ReorderedCoincidenceDuring(:)));
    RNCA = sum(sum(triu(CoincidenceAfterReordered))); % Upper diag
    RNCA = RNCA - TotalNPeaksAfter;
    
    ThisTrial(trial) = RNCA;
end
disp('');

figdist = figure;
subplot(2,2,1);hist(ThisTrial,1000);

% Now, for fitting and plotting the best fits
subplot(2,2,2);[N,X] = hist(ThisTrial,1000);
% Calculate DX
DX = X(2) - X(1); % Assumed equal for all X values

% Presumed empirical probabilities...
Phat1 = N/(sum(N));

% First bar plot of empirical probabilities
bar(X,Phat1); hold on

%
[LHat, Lci] = poissfit(ThisTrial,0.01);
RawPDF = poisspdf(round(X),LHat);
plot(round(X),RawPDF,'--');

% Optimise Graphical fit...
% This should be analytically do-able
% but since different parametric curves are being tried
% this approach is flexible.
Aopt = RawPDF * pinv(Phat1);
plot(round(X),RawPDF/Aopt,'-');hold off
title('Best Poisson model: After');
[lh,oh]=legend('','PDF fit','Fit, optimised for display');
delete(oh(4));

% Next, try a Gamma Fit
subplot(2,2,3);
bar(X,Phat1); hold on

[PARMHAT,PARMCI] = gamfit(ThisTrial,0.01);
RawPDF = gampdf(X,PARMHAT(1),PARMHAT(2));
h1=plot(X,RawPDF,'r--'); 
set(h1,'LineWidth',2);

% Optimise Graphical fit...
Aopt = RawPDF * pinv(Phat1);
h2=plot(X,RawPDF/Aopt,'r-');hold off
set(h2,'LineWidth',2);
title('Best Gamma model: After');
[lh,oh]=legend('','PDF fit','Fit, optimised for display');
delete(oh(4));

subplot(2,2,4);
% Need to re-do X axis sensibly....
Xmin = 0; Xmax = max([max(X),ThisTrial,NCA]);
NewXaxis = linspace(Xmin,Xmax*1.1,100); % 1.1 is a fudge factor


% Plot PDF curve for these parameters, extending out to the
% 10pc beyond the biggest coincident value we have found...
ScaledPDF = gampdf(NewXaxis,PARMHAT(1),PARMHAT(2))/Aopt;
plot(NewXaxis,ScaledPDF,'b');hold on
hL=line([NCA,NCA],[0 max(ScaledPDF)]);
set(hL,'Color','r','LineStyle','--'); hold off;
set(hL,'LineWidth',2);
title('Gamma fit and P(C_{obs}|Gamma;a,b)')
PCEqNCA = gampdf(NCA,PARMHAT(1),PARMHAT(2));
format shorte;
text(NCA,max(ScaledPDF),['P=',num2str(PCEqNCA)]); 
hold off;



