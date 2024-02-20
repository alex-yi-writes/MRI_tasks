clear all;

 subjects={'1_hf','2_rkc','3_ro','4_dn','5_ls','6_co','7_rr','8_ms','9_jd','10_at','11_bg','12_tm','13_bl','14_cc','15_pa','16_jt','17_sc','18_em','19_sd','20_jn','21_jy','22_mk2','23_rd','24_cb','25_sb','26_vc','27_jg','28_ns','29_jm','30_cb','31_eg','32_ps','33_ac','34_dn','35_jh','36_mf','37_ja','38_cc','40_lk','41_js','42_ah','43_sw'};

learning=[2 1 1 2 1 2 1 1 2 2 1 1 2 1 2 1 2 1 1 2 1 2 1 1 1 2 1 1 1 1];
% subjects=subjects(learning==2)

% Perf_GoWin=zeros(length(subjects),60); Perf_GoLose=zeros(length(subjects),60); Perf_NoGoWin=zeros(length(subjects),60); Perf_NoGoLose=zeros(length(subjects),60);
% CPerf_GoWin=zeros(length(subjects),60); CPerf_GoLose=zeros(length(subjects),60); CPerf_NoGoWin=zeros(length(subjects),60); CPerf_NoGoLose=zeros(length(subjects),60);
% RT_GW=zeros(length(subjects),60); RT_GL=zeros(length(subjects),60); RT_NGW=zeros(length(subjects),60); RT_NGL=zeros(length(subjects),60);

for z=1:length(subjects)
    subject=subjects{z};
    
    
    
    dirsubj=['C:\from_Ddrive\RewPunGoN\all_data\' num2str(subject)];
    
    cd(dirsubj);
    file1 = [num2str(subject) '_TaskDataLearning'];
    
    %     cuefile = ['RandCSs_' num2str(subject)];
    
    %     load(cuefile);
    %     cues(z,:)=RandCSs;
    
    load(file1);
    TaskData=TaskDataLearning;
    
    
    TaskDataLearning=TaskData;
    
    cd ..
    
    analysis
    
    Perf_GoWin(z,:)=GW;
    Perf_GoLose(z,:)=GL;
    Perf_NoGoWin(z,:)=NGW;
    Perf_NoGoLose(z,:)=NGL;
    
    CPerf_GoWin(z,:)=CGW;
    CPerf_GoLose(z,:)=CGL;
    CPerf_NoGoWin(z,:)=CNGW;
    CPerf_NoGoLose(z,:)=CNGL;
    
    RT_GW(z,:)=RTGoWin';
    RT_GL(z,:)=RTGOLose';
    RT_NGW(z,:)=RTNoGoWin';
    RT_NGL(z,:)=RTNoGoLose';
    AllData_Behav_MT(:,:,z)=TaskDataLearning;
    
    Action=TaskDataLearning(:,12);
    Action(Action==0)=2;
    AA{z}=Action;
    SS{z}=TaskDataLearning(:,2);
    RR{z}=TaskDataLearning(:,17)*2;
end
save DataModel AA RR SS
A=[sum(CPerf_GoWin,2) sum(CPerf_GoLose,2) sum(CPerf_NoGoWin,2) sum(CPerf_NoGoLose,2)];

TotalCorrect=[sum(CPerf_GoWin,2) sum(CPerf_GoLose,2) sum(CPerf_NoGoWin,2) sum(CPerf_NoGoLose,2)];

CPerf_GoWin_10=[sum(CPerf_GoWin(:,1:10),2) sum(CPerf_GoWin(:,11:20),2) sum(CPerf_GoWin(:,21:30),2) sum(CPerf_GoWin(:,31:40),2) sum(CPerf_GoWin(:,41:50),2) sum(CPerf_GoWin(:,51:60),2)];
CPerf_GoLose_10=[sum(CPerf_GoLose(:,1:10),2) sum(CPerf_GoLose(:,11:20),2) sum(CPerf_GoLose(:,21:30),2) sum(CPerf_GoLose(:,31:40),2) sum(CPerf_GoLose(:,41:50),2) sum(CPerf_GoLose(:,51:60),2)];
CPerf_NoGoWin_10=[sum(CPerf_NoGoWin(:,1:10),2) sum(CPerf_NoGoWin(:,11:20),2) sum(CPerf_NoGoWin(:,21:30),2) sum(CPerf_NoGoWin(:,31:40),2) sum(CPerf_NoGoWin(:,41:50),2) sum(CPerf_NoGoWin(:,51:60),2)];
CPerf_NoGoLose_10=[sum(CPerf_NoGoLose(:,1:10),2) sum(CPerf_NoGoLose(:,11:20),2) sum(CPerf_NoGoLose(:,21:30),2) sum(CPerf_NoGoLose(:,31:40),2) sum(CPerf_NoGoLose(:,41:50),2) sum(CPerf_NoGoLose(:,51:60),2)];

figure; subplot(2,2,1); imagesc(Perf_GoWin);
subplot(2,2,2); imagesc(Perf_GoLose);
subplot(2,2,3); imagesc(Perf_NoGoWin);
subplot(2,2,4); imagesc(Perf_NoGoLose);

figure; subplot(2,2,1); imagesc(CPerf_GoWin);
subplot(2,2,2); imagesc(CPerf_GoLose);
subplot(2,2,3); imagesc(CPerf_NoGoWin);
subplot(2,2,4); imagesc(CPerf_NoGoLose);
%
figure; plot(mean(CPerf_GoWin));
hold on; plot(mean(CPerf_GoLose),'r');
figure; plot(mean(CPerf_NoGoWin));
hold on; plot(mean(CPerf_NoGoLose),'r');
%
figure; plot(mean(CPerf_GoWin_10),'.-b');
hold on; plot(mean(CPerf_GoLose_10),'.-r');
figure; plot(mean(CPerf_NoGoWin_10),'.-b');
hold on; plot(mean(CPerf_NoGoLose_10),'.-r');
%
figure; subplot(2,2,3); hist(sum(CPerf_NoGoWin,2))
subplot(2,2,4); hist(sum(CPerf_NoGoLose,2))
subplot(2,2,2); hist(sum(CPerf_GoLose,2))
subplot(2,2,1); hist(sum(CPerf_GoWin,2))

figure; plot(mean(CPerf_GoWin_10),'s-b');
hold on; plot(mean(CPerf_GoLose_10),'v-b');
plot(mean(CPerf_NoGoWin_10),'s-r');
plot(mean(CPerf_NoGoLose_10),'v-r');

errorbar(1:6,mean(CPerf_GoWin_10),std(CPerf_GoWin_10)/sqrt(length(subjects)),std(CPerf_GoWin_10)/sqrt(length(subjects)),'b');
errorbar(1:6,mean(CPerf_GoLose_10),std(CPerf_GoLose_10)/sqrt(length(subjects)),std(CPerf_GoLose_10)/sqrt(length(subjects)),'b');
errorbar(1:6,mean(CPerf_NoGoWin_10),std(CPerf_NoGoWin_10)/sqrt(length(subjects)),std(CPerf_NoGoWin_10)/sqrt(length(subjects)),'r');
errorbar(1:6,mean(CPerf_NoGoLose_10),std(CPerf_NoGoLose_10)/sqrt(length(subjects)),std(CPerf_NoGoLose_10)/sqrt(length(subjects)),'r');

figure; bar(mean(TotalCorrect/60));
hold on; errorbar(1:4,mean(TotalCorrect/60),std(TotalCorrect/60)/sqrt(length(subjects)),std(TotalCorrect/60)/sqrt(length(subjects)));
