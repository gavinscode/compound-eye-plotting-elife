function [ arrayOut ] = loadTIFStack( address, loadSmall)
%Written by Gavin Taylor, 2017. MIT License

    cd(address); dirContents = dir; includedFiles = [];
    for i = 1:length(dirContents)
        if strfind(dirContents(i).name, 'tif')     
            includedFiles = [ includedFiles i ];
        end
    end

    warning('im reader update to _w_by_h version. can effect subscripts and plot3');
                
    counter = 1;
    numberOfImages = length(includedFiles);
    for  i = 1:numberOfImages
        if i == 1
            [temp] = imread_w_by_h(dirContents(includedFiles(i)).name, 'TIFF');
            if loadSmall
                arrayOut = zeros(size(temp,1), size(temp,2), numberOfImages, 'uint8');
                arrayOut(:,:,counter) = uint8(temp);
            else
                arrayOut = zeros(size(temp,1), size(temp,2), numberOfImages);
                arrayOut(:,:,counter) = temp;
            end
        else
            if loadSmall
                arrayOut(:,:,counter) = uint8(imread_w_by_h(dirContents(includedFiles(i)).name, 'TIFF'));
            else
                arrayOut(:,:,counter) = imread_w_by_h(dirContents(includedFiles(i)).name, 'TIFF');
            end
        end
        counter = counter + 1;
    end
end

