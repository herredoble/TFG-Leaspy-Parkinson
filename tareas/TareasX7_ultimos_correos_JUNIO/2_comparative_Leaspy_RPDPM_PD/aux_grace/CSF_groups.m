function PPMI=CSF_groups(PPMI)

ID=unique(PPMI.PATNO);

mask_Abeta_mas=[];
mask_Abeta_menos=[];
censuringTimeAbeta=[];
convertedTimeAbeta=[];

mask_Tau_mas=[];
mask_Tau_menos=[];
censuringTimeTau=[];
convertedTimeTau=[];

mask_TauAbeta_mas=[];
mask_TauAbeta_menos=[];
censuringTimeTauAbeta=[];
convertedTimeTauAbeta=[];

thres_Abeta=880; % Hasson18
thres_Tau=290; % Hasson18
thres_TauAbeta=0.33; % Hasson18

for i=1:length(ID)
    index=find(PPMI.PATNO==ID(i));
    subj=PPMI(index,:);
    
    Abeta=subj.abeta;
    Tau=subj.tau;
    TauAbeta=Tau./Abeta;
   
    %% Abeta
    mask_Abeta=~isnan(Abeta);   
    years_Abeta=subj.YEAR(mask_Abeta);    
    Abeta=Abeta(mask_Abeta);
    if(isempty(Abeta)==0)      
        if(sum(Abeta<=thres_Abeta) > 0)
            mask_Abeta_mas=[mask_Abeta_mas;index];
            convertedTimeAbeta=[convertedTimeAbeta;...
                repmat(min(years_Abeta(Abeta<=thres_Abeta)),length(index),1)];
        else
            mask_Abeta_menos=[mask_Abeta_menos;index];
            censuringTimeAbeta=[censuringTimeAbeta;...
                repmat(max(years_Abeta),length(index),1)];
           
        end
    else
        fprintf(2,'Subject %d without ABETA data\n',ID(i));
    end
    
    %% Tau
    mask_Tau=~isnan(Tau);   
    years_Tau=subj.YEAR(mask_Tau);
    Tau=Tau(mask_Tau);
    if(isempty(Tau)==0)      
        if(sum(Tau>=thres_Tau) > 0)
            mask_Tau_mas=[mask_Tau_mas;index];
            convertedTimeTau=[convertedTimeTau;...
                repmat(min(years_Tau(Tau>=thres_Tau)),length(index),1)];
        else
            mask_Tau_menos=[mask_Tau_menos;index];
            censuringTimeTau=[censuringTimeTau;...
                repmat(max(years_Tau),length(index),1)];          
        end
    else
        fprintf(2,'Subject %d without TAU data\n',ID(i));
    end
    
    %% TauAbeta
    mask_TauAbeta=~isnan(TauAbeta);   
    years_TauAbeta=subj.YEAR(mask_TauAbeta);
    TauAbeta=TauAbeta(mask_TauAbeta);
    if(isempty(TauAbeta)==0)      
        if(sum(TauAbeta>=thres_TauAbeta) > 0)
            mask_TauAbeta_mas=[mask_TauAbeta_mas;index];
            convertedTimeTauAbeta=[convertedTimeTauAbeta;...
                repmat(min(years_TauAbeta(TauAbeta>=thres_TauAbeta)),length(index),1)];
        else
            mask_TauAbeta_menos=[mask_TauAbeta_menos;index];
            censuringTimeTauAbeta=[censuringTimeTauAbeta;...
                repmat(max(years_TauAbeta),length(index),1)];          
        end
    else
        fprintf(2,'Subject %d without TAU and ABeta data\n',ID(i));
    end
    
end

PPMI.ABETAmas=zeros(size(PPMI,1),1);
PPMI.ABETAmas(mask_Abeta_mas)=1;
PPMI.ABETAconvertTime=zeros(size(PPMI,1),1);
PPMI.ABETAconvertTime(mask_Abeta_mas)=convertedTimeAbeta;
PPMI.ABETAconvertTime(mask_Abeta_menos)=censuringTimeAbeta;

PPMI.TAUmas=zeros(size(PPMI,1),1);
PPMI.TAUmas(mask_Tau_mas)=1;
PPMI.TAUconvertTime=zeros(size(PPMI,1),1);
PPMI.TAUconvertTime(mask_Tau_mas)=convertedTimeTau;
PPMI.TAUconvertTime(mask_Tau_menos)=censuringTimeTau;

PPMI.TAUABETAmas=zeros(size(PPMI,1),1);
PPMI.TAUABETAmas(mask_TauAbeta_mas)=1;
PPMI.TAUABETAconvertTime=zeros(size(PPMI,1),1);
PPMI.TAUABETAconvertTime(mask_TauAbeta_mas)=convertedTimeTauAbeta;
PPMI.TAUABETAconvertTime(mask_TauAbeta_menos)=censuringTimeTauAbeta;


end