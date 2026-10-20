function t_mae_auc=...
    get_Future_MAE_Constant_v1(y,labels,ages,future_timepoints)

btstrp=1;
I=size(ages,1);
K=size(y,1);

t_mae_auc=nan(I,K+btstrp+2,length(future_timepoints));

for j=1:length(future_timepoints)
    mask_future= ~isnan(ages);
    timepoints=nan(I,1);
    ages_aux=ages;
    for i=1:I
        n_samples=sum(mask_future(i,:));
        if((n_samples/2)>=future_timepoints(j))
            mask_future(i,n_samples-future_timepoints(j)+1:n_samples)=0;    
            timepoints(i)=ages(i,n_samples)-ages(i,n_samples-future_timepoints(j));
            ages_aux(i,n_samples-future_timepoints(j)+1:n_samples-1)=nan;
        end
    end
    
    [label_test_fit,ae_test]=get_Future_MAE_Constant(y,labels,ages_aux,mask_future,btstrp);

    mae_subj=nan(I,size(y,1));
    label_subj=nan(I,btstrp+1);
    for i=1:I
        mask_subj_future = mask_future(i,:)==0 & ~isnan(ages_aux(i,:));
        if(sum(mask_subj_future)>0)
            mae_subj(i,:)= squeeze(ae_test(:,i,mask_subj_future,:))';
            label_subj(i,:)=[labels(i,mask_subj_future),squeeze(label_test_fit(i,mask_subj_future,:))'];
        end
    end
    t_mae_auc(:,:,j)=[timepoints,mae_subj,label_subj];
end

end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [label_test_fit,ae_test]=...
    get_Future_MAE_Constant(y_test,labels_test,ages_test,mask_future,btstrp) 
%% Constant Predictive DPM Testing

% Estimating test subject parameters using biomaker values and ages
[K2, I2, J2] = size(y_test); % number of biomarkers, subjects, and visits
label_test_fit = nan(I2, J2, btstrp);
ae_test= nan(K2,I2,J2,btstrp);

[label_t, y_test_fit]=constant_test_future(y_test,labels_test,ages_test, mask_future);
ae_test(:,:,:,1) = abs(y_test - y_test_fit); % absolute errors
label_test_fit(:,:,:,1)=label_t;

end



function [l_fit, y_fit]=constant_test_future(y,l,t,mask_future)
[K, I, J] = size(y);
y_fit = NaN(K, I, J);
l_fit = NaN(I, J);

for i = 1 : I
    y_i = squeeze(y(:, i, :));
    l_i = squeeze(l(i, :));   

    %%%%%%%%%%%%%%%%%%%%%%%%%%% future
    y_i=y_i(:,mask_future(i,:));
    l_i=l_i(mask_future(i,:));   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    mask_s_i = mask_future(i,:)==0 & ~isnan(t(i,:));
    if(sum(mask_s_i)>0)
        y_end = ultimosNoNaN(y_i);
        l_end = ultimosNoNaN(l_i);

        y_fit(:,i,mask_s_i)=y_end;
        l_fit(i,mask_s_i)=l_end;
    end
   
end


end


function ultimos_valores = ultimosNoNaN(matriz)
    [filas, ~] = size(matriz);
    ultimos_valores = NaN(filas, 1);
    
    for i = 1:filas
        fila_actual = matriz(i,:);
        ultimo_indice = find(~isnan(fila_actual), 1, 'last');
        if ~isempty(ultimo_indice)
            ultimos_valores(i) = fila_actual(ultimo_indice);
        end
    end
end

