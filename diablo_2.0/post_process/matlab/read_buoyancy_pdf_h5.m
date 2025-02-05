% filename='../buoypdf.h5';
filename=[base_dir '/buoypdf.h5'];


for i=1:nk
    if i<10
        binname=['/bins/000' num2str(i)];
        dbinname=['/dbins/000' num2str(i)];
        countname=['/counts/000' num2str(i)];
    elseif (i>=10)&&(i<100)
        binname=['/bins/00' num2str(i)];
        dbinname=['/dbins/00' num2str(i)];
        countname=['/counts/00' num2str(i)];  
    elseif (i>=100)&&(i<1000)
        binname=['/bins/0' num2str(i)];
        dbinname=['/dbins/0' num2str(i)];
        countname=['/counts/0' num2str(i)]; 
    else
        binname=['/bins/' num2str(i)];
        dbinname=['/dbins/' num2str(i)];
        countname=['/counts/' num2str(i)];        
    end
    

    bins=h5read(filename,binname);
    dbins=h5read(filename,dbinname);
    counts=h5read(filename,countname);
    time(i)=h5readatt(filename,binname,'Time');
    
    binval=bins+dbins/2;
    
    normalization_factor=trapz(binval,counts);
    buoyancy_pdf(:,i)=counts/normalization_factor;
end

return
figure;plot(binval,buoyancy_pdf(:,1),binval,buoyancy_pdf(:,end));
saveas(gcf,[base_dir '/fig_buoypdf.fig'])