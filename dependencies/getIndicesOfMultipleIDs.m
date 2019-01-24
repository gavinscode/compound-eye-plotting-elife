function [ foundIndices ] = getIndicesOfMultipleIDs( image, IDs )
%Written by Gavin Taylor, 2017. MIT License

    foundIndices = [];

    for i = 1:length(IDs)
        foundIndices = [foundIndices find(image == IDs(i))'];
    end

    foundIndices = foundIndices';
end

