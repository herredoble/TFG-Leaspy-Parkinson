function PPMI=estimateCogState(PPMI,ID_sinCog)

ID_sinCog=ID_sinCog(:,2);
for i=1:length(ID_sinCog)
    index=find(PPMI.PATNO==ID_sinCog(i));
    MCI_test=PPMI.MCI_testscores_end(index(1));
    MCI_test_ConverTime=PPMI.convertTimeMCI_test(index(1));
    np1cog=PPMI.np1cog_end(index(1));
    np1cog_ConvertTime=PPMI.convertTimeNp1cog(index(1));
    if(((MCI_test>0) && (np1cog>=1)) )
        PPMI.cogstate_end(index)=2;
        PPMI.convertTimeCogstate(index)=(MCI_test_ConverTime+...
            np1cog_ConvertTime)/2;
    elseif(((MCI_test==0) && (np1cog<2)) || (np1cog==0))
        PPMI.cogstate_end(index)=1;
        PPMI.convertTimeCogstate(index)=(MCI_test_ConverTime+...
            np1cog_ConvertTime)/2;
        
    else
        fprintf('Discrepance criteria between MCI-Test vs np1cog from %d\n',...
            ID_sinCog(i));
    end
        
        
    
    
end


end