function linkCap = getLinkCapacity(linkCapacitiesAvailable)

    % AUTHOR: Afaf
    % LAST MODIFIED: 27/10/2020
    
    % Choose randomly a value for link capacity from the available ones
    
    % OUTPUT:
    % linkCap : the capacity of the link in Mbps

    linkCap = linkCapacitiesAvailable(randperm(length(linkCapacitiesAvailable),1));


end