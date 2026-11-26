close all

legends=['A';'B';'C';'D'];

onset_thr=50:200:4*200;

for cc1=1:4

    figure(1)
    subplot(2,2,cc1)

    %function plot_SW_subepidemic(flag1_pass,r_pass,p_pass,a_pass,K_pass,q_pass,npatches_pass,onset_thr_pass,typedecline1_pass,I0_pass,windowsize1_pass)

    curve=plot_SW_subepidemic(1,0.18,0.9,[],1000,0.24,8,onset_thr(cc1),1,5,200)

    title(strcat('C_{thr}=',num2str(onset_thr(cc1))))

    text(5,max(curve)-1,strcat('\fontsize{24}',legends(cc1)))


end
