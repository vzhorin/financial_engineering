time_to_maturity = 1;
delta = .4;   

num_sim = 5000;

num_bonds(1) = 5;
num_bonds(2) = 10;
num_bonds(3) = 50;

def_corr(1) = 0.;
def_corr(2) = 0.2;
def_corr(3) = 0.7;

def_prob(1) = 0.01;
def_prob(2) = 0.1;
def_prob(3) = 0.4;

total_bonds = length(num_bonds);
total_corr= length(def_corr);
total_prob = length(def_prob);

cdo_equity = 0.05; %upper limit
cdo_mezz = 0.15;

x = randn(num_bonds(end),num_sim); 

total_loss_mean = zeros(total_bonds,total_prob*total_corr);
total_loss_std = zeros(total_bonds,total_prob*total_corr);

equity_loss_mean = zeros(total_bonds,total_prob*total_corr);
equity_loss_std = zeros(total_bonds,total_prob*total_corr);

mezz_loss_mean = zeros(total_bonds,total_prob*total_corr);
mezz_loss_std = zeros(total_bonds,total_prob*total_corr);

senior_loss_mean = zeros(total_bonds,total_prob*total_corr);
senior_loss_std = zeros(total_bonds,total_prob*total_corr);


for i = 1:total_bonds
    i_temp = 1;

    for j=1:total_corr
        sigma_mat = diag(ones(num_bonds(i),1));
        sigma_mat(sigma_mat==0) = def_corr(j);

        for k = 1:total_prob
            % see Teaching Note 8, page 519
            C = chol(sigma_mat)'; 
            r = C*x(1:num_bonds(i),:);     
 
            U = normcdf(r);      
            tau = -log(U)/def_prob(k);      

            sample_loss = (1-delta)*sum(tau<time_to_maturity)/num_bonds(i);   %in percents to portfolio value
            if (i == 3 && k ==3) 
            figure(i*100+j*10+k); hist(sample_loss, 100); 
            title(['Loss distibution for ' num2str(num_bonds(i)) ' bonds with pdef=' num2str(def_prob(k)) ' and corr= ' num2str(def_corr(j)) ], ...
                'FontSize',16, 'FontWeight', 'bold');
            % hist(sample_loss)
            end
            total_loss_mean(i,i_temp)=mean(sample_loss);
            total_loss_std(i,i_temp)=std(sample_loss);
            
            % either take percentage loss or full loss of layer
            sample_loss_equity = [sample_loss(sample_loss <= cdo_equity)/cdo_equity ...
               sample_loss(sample_loss > cdo_equity)./sample_loss(sample_loss > cdo_equity)];
           
            % either zero, or take percentage loss or full loss of layer
            sample_loss_mezz = [sample_loss(sample_loss <= cdo_equity)-sample_loss(sample_loss <= cdo_equity) ...
                (sample_loss(sample_loss <= cdo_mezz & sample_loss > cdo_equity)-cdo_equity)/(cdo_mezz-cdo_equity) ...
               sample_loss(sample_loss > cdo_mezz)./sample_loss(sample_loss > cdo_mezz)];
           
            % either zero, or take percentage loss, never larger than 100% 
            sample_loss_senior = [sample_loss(sample_loss <= cdo_mezz)-sample_loss(sample_loss <= cdo_mezz) ...
                (sample_loss(sample_loss > cdo_mezz+cdo_equity)-cdo_equity-cdo_mezz)/(1-cdo_mezz-cdo_equity)];

           
            equity_loss_mean(i,i_temp)=mean(sample_loss_equity);
            equity_loss_std(i,i_temp)=std(sample_loss_equity);
            
            mezz_loss_mean(i,i_temp)=mean(sample_loss_mezz);
            mezz_loss_std(i,i_temp)=std(sample_loss_mezz);
            
            senior_loss_mean(i,i_temp)=mean(sample_loss_senior);
            senior_loss_std(i,i_temp)=std(sample_loss_senior);
            
            i_temp = i_temp+1;

        end 
    end 
end 

for i = 1:total_bonds
 if (i == 1)     
figure(i); 
bar3(reshape(total_loss_mean(i,:)*100,3,3)', 0.25);
ylabel('default correlation','FontSize',16, 'FontWeight', 'bold'); 
xlabel('default probability','FontSize',16, 'FontWeight', 'bold'); 
zlabel('mean loss %','FontSize',16, 'FontWeight', 'bold'); 
title(['Mean loss for ' num2str(num_bonds(i)) ' bonds'],'FontSize',24, 'FontWeight', 'bold')

figure(i+10); 
bar3(reshape(total_loss_std(i,:),3,3)', 0.25);
ylabel('default correlation','FontSize',16, 'FontWeight', 'bold'); 
xlabel('default probability','FontSize',16, 'FontWeight', 'bold'); 
zlabel('\sigma ','FontSize',16, 'FontWeight', 'bold'); 
title(['St.dev. of losses for ' num2str(num_bonds(i)) ' bonds'],'FontSize',24, 'FontWeight', 'bold')

figure(i+20); 
bar3(reshape(equity_loss_mean(i,:)*100,3,3)', 0.25);
ylabel('default correlation','FontSize',16, 'FontWeight', 'bold'); 
xlabel('default probability','FontSize',16, 'FontWeight', 'bold'); 
zlabel('equity loss %','FontSize',16, 'FontWeight', 'bold'); 
title(['Equity loss for ' num2str(num_bonds(i)) ' bonds'],'FontSize',24, 'FontWeight', 'bold')

figure(i+30); 
bar3(reshape(mezz_loss_mean(i,:)*100,3,3)', 0.25);
ylabel('default correlation','FontSize',16, 'FontWeight', 'bold'); 
xlabel('default probability','FontSize',16, 'FontWeight', 'bold'); 
zlabel('mezz loss %','FontSize',16, 'FontWeight', 'bold'); 
title(['Mezzanine loss for ' num2str(num_bonds(i)) ' bonds'],'FontSize',24, 'FontWeight', 'bold')

figure(i+40); 
bar3(reshape(senior_loss_mean(i,:)*100,3,3)', 0.25);
ylabel('default correlation','FontSize',16, 'FontWeight', 'bold'); 
xlabel('default probability','FontSize',16, 'FontWeight', 'bold'); 
zlabel('senior loss %','FontSize',16, 'FontWeight', 'bold'); 
title(['Senior loss for ' num2str(num_bonds(i)) ' bonds'],'FontSize',24, 'FontWeight', 'bold')
 end
end
