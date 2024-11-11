function    varargout = num2csl( num )
% num2csl converts a numerical vector into a comma separated list
    assert( isnumeric(num), 'num2csl:WrongClass', 'Input must be numeric' )
    cac = num2cell( num );
    varargout = cac( 1 : nargout );
end