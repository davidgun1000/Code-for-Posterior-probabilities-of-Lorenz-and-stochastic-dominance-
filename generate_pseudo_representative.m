function [obs_FPBB]=generate_pseudo_representative(income_data,sampling_weight,n)

    N=600000; %this must be set to a really big number
    l_boot=zeros(n,1);
    Nnn = (N-n)/n;
    weight_selected_sample = sampling_weight.*(N/sum(sampling_weight));
    com = [income_data weight_selected_sample];
    
    for k=1:(N-n)
        
         newweights_num=com(:,2)-1+l_boot.*Nnn;
         newweights_den=(N-n)+(k-1)*Nnn;
         newweights=newweights_num./newweights_den;
         [y_select(k,1),idx]=datasample(com(:,1),1,'replace',true,'weights',newweights);
         lk=zeros(n,1);
         lk(idx,1)=1;
         l_boot=l_boot+lk;
     end
     
     FPBB_population=[com(:,1);y_select];
     obs_FPBB=datasample(FPBB_population,n,'replace',true);


end