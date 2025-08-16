function dateStrings = create_quarterly_date_strings(startYear,...
    startQuarter,nPeriods)

nYears = ceil(nPeriods/4)+1;
quartersInAYear = [1; 2; 3; 4];
quarterSequence = repmat(quartersInAYear,nYears,1);
firstQuarterToSelect = find(quarterSequence==startQuarter,1);
lastQuarterToSelect = firstQuarterToSelect+nPeriods-1;
quarterSequence = ...
    quarterSequence(firstQuarterToSelect:lastQuarterToSelect,1);
quarterStrings = [repmat('Q',nPeriods,1) num2str(quarterSequence)];
yearAdditions = repmat([1; 0; 0; 0],nYears,1);
years = startYear - 1 + cumsum(yearAdditions);
yearStrings = num2str(years(firstQuarterToSelect:lastQuarterToSelect,1));
dateStrings = [yearStrings quarterStrings];
